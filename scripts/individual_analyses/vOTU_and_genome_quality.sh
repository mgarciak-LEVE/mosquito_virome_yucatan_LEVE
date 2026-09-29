#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script assign viral OTUs and assess genome quality
# 28/09/2026
# Ver. 1.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directories
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/contig_concatenation"

# STAR assembly directory (parent indicates mapping tool)
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/STAR"

# Bowtie assembly directory (parent indicates mapping tool)
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/Bowtie"

# Output directory
OUTPUT_DIR_OTU="${PROJECT_SCRATCH}/results/vOTU"
OUTPUT_DIR_COMPLETENESS="${PROJECT_SCRATCH}/results/completeness"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
CHECKV_CONTAINER="${CONTAINERS}/checkv.sif"
CDHIT_CONTAINER="${CONTAINERS}/cdhit_4.8.1--c697693.sif"

# Databases paths
DB_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/databases"
CHECKV_DB="${DB_BASE}/checkv_db/checkv-db-v1.5"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=12
CORES=6

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Identification Pipeline"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  CD-HIT container: ${CDHIT_CONTAINER}"
echo "  Input directories (mapping_tool/assembler/read_type/sample):"
echo "    - STAR: ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output for vOTUs: ${OUTPUT_DIR_OTU}"
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
declare -a CONTIG_FILE_ARRAY

echo "Building sample-assembler list..."

# Function to add sample-assembler combination
add_combination() {
    local mapping_tool=$1
    local read_type=$2         
    local assembler=$3
    local sample=$4
    local contig_file=$5

    if [[ -f "$contig_file" ]] && [[ -s "$contig_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        ASSEMBLER_ARRAY+=("$assembler")
        CONTIG_FILE_ARRAY+=("$contig_file")
        echo "  Added: $mapping_tool/$read_type/$assembler/$sample - $(basename "$contig_file")"
    else
        echo "  WARNING: Assembly file not found or empty for $mapping_tool/$read_type/$assembler/$sample: $contig_file"
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
                contig_file="${sample_dir}/all.fa"
                add_combination "$mapping_tool" "$read_type" "$assembler" "$sample" "$contig_file"

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
CONTIG_FILE="${CONTIG_FILE_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Sample: ${SAMPLE}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Contig file: ${CONTIG_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####=======================================####
####      CREATE OUTPUT DIR FOR vOTUS      ####
####=======================================####

# Output structure: mapping_tool/read_type/assembler/sample
SAMPLE_OUTPUT_OTU="${OUTPUT_DIR_OTU}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"

####=========================####
####    vOTU assignment      ####
####=========================####

run_cdhit() {
    echo ""
    echo "--- Running CDHIT ---"
    local output_dir="${SAMPLE_OUTPUT_OTU}/"
    mkdir -p "$output_dir"
    
    # Check if it has been completed already 
    if [[ -f "${output_dir}/vOTUs.fa" ]] && [[ -s "${output_dir}/vOTUs.fa" ]]; then 
        echo "CD-hit already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi


    # Verify the Checkv container exists
    if [[ ! -f "${CDHIT_CONTAINER}" ]]; then
        echo "  ERROR: Container not found at ${CDHIT_CONTAINER}"
        return 1
    fi

    # Run cd-hit using apptainer container
    apptainer exec \
        --bind "$(dirname "$CONTIG_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT_OTU}":/output \
        "${CDHIT_CONTAINER}" \
        cd-hit-est \
        -i /input/all.fa \
        -o /output/vOTUs.fa \
        -c 0.95 \
        -aS 0.85 \
        -T "${CORES}" \
        2>&1 | tee "${output_dir}/cdhit.log"

    local cdhit_exit=${PIPESTATUS[0]}
    
    if [[ ${cdhit_exit} -ne 0 ]]; then
        echo "  CD-HIT failed (exit ${cdhit_exit})"
        return 1
    fi

}


####============================================####
####      CREATE OUTPUT DIR FOR COMPLETENESS    ####
####============================================####

# Output structure: mapping_tool/read_type/assembler/sample
SAMPLE_OUTPUT_COMPLETENESS="${OUTPUT_DIR_COMPLETENESS}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"

####==============================####
####      Genome Completeness     ####
####==============================####

run_checkv() {
    echo ""
    echo "--- Running CheckV---"
    local output_dir="${SAMPLE_OUTPUT_COMPLETENESS}/"
    mkdir -p "$output_dir"
    
    # Check if it has been completed already 
    if [[ -f "${output_dir}/quality_summary.tsv" ]]; then 
        echo "CheckV already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi

    # Check if database exists
    if [[ ! -d "${CHECKV_DB}" ]]; then
        echo "  WARNING: Checkv database not found or incomplete at ${CHECKV_DB}"
        echo "  Skipping CheckV completeness analysis"
        return 0
    fi

    # Verify the Checkv container exists
    if [[ ! -f "${CHECKV_CONTAINER}" ]]; then
        echo "  ERROR: Container not found at ${CHECKV_CONTAINER}"
        return 1
    fi

    # Run checkv using apptainer container
    apptainer exec \
        --bind "${SAMPLE_OUTPUT_OTU}":/input:ro \
        --bind "${SAMPLE_OUTPUT_COMPLETENESS}":/output \
        --bind "${CHECKV_DB}":/database \
        "${CHECKV_CONTAINER}" \
        checkv end_to_end \
        /input/vOTUs.fa \
        /output \
        -d /database \
        -t "${THREADS}" \
        2>&1 | tee "${output_dir}/checkv.log"

    local checkv_exit=${PIPESTATUS[0]}
    
    if [[ ${checkv_exit} -ne 0 ]]; then
        echo "  Checkv failed (exit ${checkv_exit})"
        return 1
    fi

}

####==================================####
####          RUN ALL TOOLS          ####
####==================================####

echo ""
echo "========================================="
echo "  Starting viral OTUs tools"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
echo "========================================="

# Run each tool
if ! run_cdhit; then
    echo "ERROR: CD-HIT failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
    tg_send "ERROR: CD-HIT failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

if ! run_checkv; then
    echo "ERROR: CheckV failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
    tg_send "ERROR: CheckV failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true
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
echo "  Output: ${SAMPLE_OUTPUT_OTU} & ${SAMPLE_OUTPUT_COMPLETENESS}"
echo "========================================="

tg_send "Completed: ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

