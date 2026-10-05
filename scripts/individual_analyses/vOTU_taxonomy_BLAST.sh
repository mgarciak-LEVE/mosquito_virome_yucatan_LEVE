#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script assign taxonomy to viral OTUs and individual ORFs with geNOMAD and BLAST.
# 29/09/2026
# Ver. 1.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

set -uo pipefail

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directories ---
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/vOTU"
# STAR assembly directory (parent indicates mapping tool)
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/STAR"
# Bowtie assembly directory (parent indicates mapping tool)
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/Bowtie"

# Output directories ---
OUTPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_DIR_TAXONOMY="${OUTPUT_DIR_BASE}/taxonomy"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
PRODIGAL_CONTAINER="${CONTAINERS}/prodigal_2.6.3--h577a1d6_11.sif"
GENOMAD_CONTAINER="${CONTAINERS}/genomad_fixed.sif"
DIAMOND_CONTAINER="${CONTAINERS}/diamond.sif" 

# Databases paths
DB_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/databases"
GENOMAD_DB="${DB_BASE}/genomad_db"
REFSEQ_DB="${DB_BASE}/blast_db/viral_refseq/ref_viruses_rep_genomes"
DIAMOND_DB="${DB_BASE}/diamond_db/diamond_db_build"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"
module load blast/2.15.0--pl5321h6f7f691_1 || echo "BLAST module not available"
module load conda || echo "Conda  module not available" 

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genomad

THREADS=36

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Identification Pipeline"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  Prodigal container: ${PRODIGAL_CONTAINER}"
echo "  geNOMAD container: ${GENOMAD_CONTAINER}"
echo "  Input directories (mapping_tool/assembler/read_type/sample):"
echo "    - STAR: ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output for geNOMAD: ${OUTPUT_DIR_TAXONOMY}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting viral identification pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE-ASSEMBLER LIST  ####
####==================================####

# Create arrays to store sample-mapping-assembler combinations
declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a READ_TYPE_ARRAY
declare -a ASSEMBLER_ARRAY
declare -a OTU_FILE_ARRAY

echo "Building sample-assembler list..."

### /lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/vOTU/STAR/paired/MEGAhit/PM3050

# Function to add sample-assembler combination
add_combination() {
    local mapping_tool=$1
    local read_type=$2         
    local assembler=$3
    local sample=$4
    local otu_file=$5

    if [[ -f "$otu_file" ]] && [[ -s "$otu_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        ASSEMBLER_ARRAY+=("$assembler")
        OTU_FILE_ARRAY+=("$otu_file")
        echo "  Added: $mapping_tool/$read_type/$assembler/$sample - $(basename "$otu_file")"
    else
        echo "  WARNING: Assembly file not found or empty for $mapping_tool/$read_type/$assembler/$sample: $otu_file"
    fi
}

# Function to process a mapping tool directory
process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2

    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    # Iterate read_type first
    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        read_type=$(basename "$read_type_dir")
        echo "    Processing read type: $read_type"

        for assembler_dir in "$read_type_dir"/*/; do
            [[ -d "$assembler_dir" ]] || continue
            assembler=$(basename "$assembler_dir")
            echo "      Processing assembler: $assembler"

            for sample_dir in "$assembler_dir"/*/; do
                [[ -d "$sample_dir" ]] || continue
                sample=$(basename "$sample_dir")
                otu_file="${sample_dir}/vOTUs.fa"
                add_combination "$mapping_tool" "$read_type" "$assembler" "$sample" "$otu_file"

                done
            done
        done
}

# Process STAR assemblies
process_mapping_tool "STAR" "$STAR_ASSEMBLY_DIR"

# Process Bowtie assemblies
process_mapping_tool "Bowtie" "$BOWTIE_ASSEMBLY_DIR"

# Total number of jobs
TOTAL_JOBS=${#SAMPLE_ARRAY[@]}

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No valid assembly files found"
    tg_send "ERROR: No valid assembly files found for identification" 2>/dev/null || true
    exit 1
fi

star_count=0; bowtie_count=0
for m in "${MAPPING_ARRAY[@]}"; do
    [[ "$m" == "STAR"   ]] && ((star_count++))   || true
    [[ "$m" == "Bowtie" ]] && ((bowtie_count++)) || true
done

echo ""
echo "========================================="
echo "  Total jobs: ${TOTAL_JOBS}"
echo "  Unique samples: $(printf '%s\n' "${SAMPLE_ARRAY[@]}" | sort -u | wc -l)"
echo "  Summary by mapping tool:"
echo "    - STAR:   ${star_count}"
echo "    - Bowtie: ${bowtie_count}"
echo "========================================="

####==================================####
####          GET ARRAY TASK          ####
####==================================####

# LSB_JOBINDEX is the array task ID (1-based)
JOB_INDEX=$((LSB_JOBINDEX - 1))

if [[ $JOB_INDEX -ge ${#SAMPLE_ARRAY[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#SAMPLE_ARRAY[@]})"
    exit 1
fi

# Get sample, mapping tool, and assembler for this job
SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"
MAPPING_TOOL="${MAPPING_ARRAY[$JOB_INDEX]}"
READ_TYPE="${READ_TYPE_ARRAY[$JOB_INDEX]}"       
ASSEMBLER="${ASSEMBLER_ARRAY[$JOB_INDEX]}"
OTU_FILE="${OTU_FILE_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Sample: ${SAMPLE}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Contig file: ${OTU_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####=======================================####
####      CREATE OUTPUT DIR FOR geNOMAD    ####
####=======================================####

# Output structure: mapping_tool/read_type/assembler/sample
SAMPLE_OUTPUT_GENOMAD="${OUTPUT_DIR_TAXONOMY}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"

####=========================####
####    GENOMAD TAXONOMY     ####
####=========================####

run_genomad() {
    echo ""
    echo "--- Running geNOMAD annotate for taxonomy ---"
    local output_dir="${SAMPLE_OUTPUT_GENOMAD}"
    mkdir -p "$output_dir"
    
    if [[ -f "${output_dir}/vOTUs_annotate/vOTUs_taxonomy.tsv" ]] && \
       [[ -s "${output_dir}/vOTUs_annotate/vOTUs_taxonomy.tsv" ]]; then 
        echo "geNOMAD already completed for ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi

    echo "  Running genomad end-to-end..."
    apptainer exec \
        --bind "$(dirname "$OTU_FILE")":/input:ro \
        --bind "${output_dir}":/output \
        --bind "${GENOMAD_DB}":/database \
        "${GENOMAD_CONTAINER}" \
        genomad end-to-end \
        /input/vOTUs.fa \
        /output \
        /database \
        --threads "${THREADS}" \
        --splits 8 \
        --cleanup \
        2>&1 | tee "${output_dir}/genomad.log"

    local genomad_exit=${PIPESTATUS[0]}
    
    if [[ ${genomad_exit} -ne 0 ]]; then
        echo "  geNOMAD failed (exit ${genomad_exit})"
        return 1
    fi
    
    local taxonomy_file="${output_dir}/vOTUs_annotate/vOTUs_taxonomy.tsv"
    if [[ -f "${taxonomy_file}" ]]; then
        echo "  Top viral families:"
        awk -F'\t' 'NR>1 {print $NF}' "${taxonomy_file}" | \
            grep -oP '(?<=;)[^;]+$' | \
            sort | uniq -c | sort -rn | head -20
    else
        echo "  WARNING: Taxonomy file not found at ${taxonomy_file}"
    fi
}

####=======================####
####    BLAST TAXONOMY     ####
####=======================####

##  STEP 1 - PRODIGAL ORF PREDICTION  ##

# Output structure: mapping_tool/read_type/assembler/sample
SAMPLE_OUTPUT_PRODIGAL="${OUTPUT_DIR_TAXONOMY}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}/ORFs"

run_prodigal() {
    echo ""
    echo "--- Running Prodigal ---"
    local output_dir="${SAMPLE_OUTPUT_PRODIGAL}"
    mkdir -p "$output_dir"
    
    # Check if it has been completed already 
    if [[ -f "${output_dir}/vOTUs_orfs.faa" ]] && [[ -s "${output_dir}/vOTUs_orfs.faa" ]]; then 
        echo "Prodigal already completed for ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi


    # Verify the Checkv container exists
    if [[ ! -f "${PRODIGAL_CONTAINER}" ]]; then
        echo "  ERROR: Container not found at ${PRODIGAL_CONTAINER}"
        return 1
    fi

    # Run cd-hit using apptainer container
    apptainer exec \
        --bind "$(dirname "$OTU_FILE")":/input:ro \
        --bind "${output_dir}":/output \
        "${PRODIGAL_CONTAINER}" \
        prodigal \
        -i /input/vOTUs.fa \
        -a /output/vOTUs_orfs.faa \
        -d /output/vOTUs_orfs.fna \
        -p meta \
        -m -c \
        -f gff \
        -o /output/vOTUs_prodigal.gff \
        2>&1 | tee "${output_dir}/prodigal.log"

    local prodigal_exit=${PIPESTATUS[0]}
    
    if [[ ${prodigal_exit} -ne 0 ]]; then
        echo "  Prodigal failed (exit ${prodigal_exit})"
        return 1
    fi

    #  Check how many ORFs were predicted:
    for prodigal_orfs in $output_dir; do
        echo "  ORFs predicted: $(grep -c '^>' "${output_dir}/vOTUs_orfs.faa")"
    done
}

##  STEP 2 - ORF BLAST  ##

SAMPLE_OUTPUT_BLAST="${OUTPUT_DIR_TAXONOMY}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}/blast"

run_blastp() {
    echo ""
    echo "--- Running BLASTp ---"
    local output_dir="${SAMPLE_OUTPUT_BLAST}"
    mkdir -p "$output_dir"

    local orf_file="${SAMPLE_OUTPUT_PRODIGAL}/vOTUs_orfs.faa"
    local out_file="${output_dir}/vOTUs_orfs_vs_refseq.tsv"

    # Skip if already done
    if [[ -f "$out_file" ]] && [[ -s "$out_file" ]]; then
        echo "BLASTp already completed for ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi

    # Need Prodigal ORFs first
    if [[ ! -s "$orf_file" ]]; then
        echo "  ERROR: ORF file not found: $orf_file"
        return 1
    fi

    # Verify BLAST DB exists
    # makeblastdb produces files like <prefix>.phr, .pin, .psq
    if [[ ! -f "${REFSEQ_DB}.phr" ]] && [[ ! -f "${REFSEQ_DB}.pin" ]]; then
        echo "  ERROR: BLAST DB not found at ${REFSEQ_DB} (expected ${REFSEQ_DB}.phr/.pin/.psq)"
        return 1
    fi

    blastp \
        -query "$orf_file" \
        -db "${REFSEQ_DB}" \
        -out "$out_file" \
        -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames scomnames" \
        -evalue 1e-5 \
        -max_target_seqs 5 \
        -num_threads "${THREADS}" \
        2>&1 | tee "${output_dir}/blastp.log"

    local blastp_exit=${PIPESTATUS[0]}
    if [[ ${blastp_exit} -ne 0 ]]; then
        echo "  BLASTp failed (exit ${blastp_exit})"
        return 1
    fi

    echo "  ORFs with hits: $(cut -f1 "$out_file" | sort -u | wc -l)"
    echo "  Total hits:     $(wc -l < "$out_file")"
}



####==================================####
####          RUN ALL TOOLS          ####
####==================================####

echo ""
echo "========================================="
echo "  Starting viral tsxonomy tools"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
echo "========================================="

# Run each tool
if ! run_genomad; then
    echo "ERROR: geNOMAD failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
    tg_send "ERROR: geNOMAD failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

if ! run_prodigal; then
    echo "ERROR: Prodigal failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
    tg_send "ERROR: Prodigal failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

if ! run_blastp; then
    echo "ERROR: BLASTp failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
    tg_send "ERROR: BLASTp failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

####==================================####
####          COMPLETION              ####
####==================================####

echo ""
echo "========================================="
echo "  Completed job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${SAMPLE_OUTPUT_GENOMAD} & ${SAMPLE_OUTPUT_PRODIGAL} & ${SAMPLE_OUTPUT_BLAST}"
echo "========================================="

tg_send "Completed: ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

