#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to identify assembled contigs
# 06/08/2026
# Ver. 1.0.0 (farm-ready)

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
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/assembly"

# STAR assembly directory (parent indicates mapping tool)
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/star_assembly"

# Bowtie assembly directory (parent indicates mapping tool)
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/bowtie_assembly"

# Output directory
OUTPUT_DIR="${PROJECT_SCRATCH}/results/identification"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
VIRSORTER2_CONTAINER="${CONTAINERS}/spades_3.15.5.sif"
DEEPVIRFINDER_CONTAINER="${CONTAINERS}/megahit_1.2.9.sif"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Identification Pipeline"
echo "  Project: ${PROJECT_NAME}"
echo "  Input directories (mapping_tool/assembler/sample):"
echo "    - STAR: ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting viral identification pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE-ASSEMBLER LIST  ####
####==================================####

# Create arrays to store sample-mapping-assembler combinations
declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a ASSEMBLER_ARRAY
declare -a ASSEMBLY_FILE_ARRAY

echo "Building sample-assembler list..."

# Function to add sample-assembler combination
add_combination() {
    local mapping_tool=$1
    local assembler=$2
    local sample=$3
    local assembly_file=$4
    
    if [[ -f "$assembly_file" ]] && [[ -s "$assembly_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        ASSEMBLER_ARRAY+=("$assembler")
        ASSEMBLY_FILE_ARRAY+=("$assembly_file")
        echo "  Added: $mapping_tool/$assembler/$sample - $(basename "$assembly_file")"
    else
        echo "  WARNING: Assembly file not found or empty for $mapping_tool/$assembler/$sample: $assembly_file"
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
    
    # For each assembler under this mapping tool
    for assembler_dir in "$base_dir"/*/; do
        if [[ -d "$assembler_dir" ]]; then
            assembler=$(basename "$assembler_dir")
            echo "    Processing assembler: $assembler"
            
            # For each sample under this assembler
            for sample_dir in "$assembler_dir"/*/; do
                if [[ -d "$sample_dir" ]]; then
                    sample=$(basename "$sample_dir")
                    
                    # Determine which assembly file to use based on assembler
                    case $assembler in
                        "MEGAhit")
                            assembly_file="${sample_dir}/final.contigs.fa"
                            ;;
                        "metaSPAdes"|"metaviralSPAdes")
                            assembly_file="${sample_dir}/contigs.fasta"
                            ;;
                        "rnaSPAdes")
                            if [[ -f "${sample_dir}/hard_filtered_transcripts.fasta" ]]; then
                                assembly_file="${sample_dir}/hard_filtered_transcripts.fasta"
                            else
                                assembly_file="${sample_dir}/transcripts.fasta"
                            fi
                            ;;
                        *)
                            echo "      WARNING: Unknown assembler: $assembler"
                            continue
                            ;;
                    esac
                    
                    add_combination "$mapping_tool" "$assembler" "$sample" "$assembly_file"
                fi
            done
        fi
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

echo ""
echo "========================================="
echo "  Total jobs: ${TOTAL_JOBS}"
echo "  Unique samples: $(printf '%s\n' "${SAMPLE_ARRAY[@]}" | sort -u | wc -l)"
echo "  Summary by mapping tool:"
echo "    - STAR: $(grep -c "STAR" <<< "${MAPPING_ARRAY[@]}" || echo "0")"
echo "    - Bowtie: $(grep -c "Bowtie" <<< "${MAPPING_ARRAY[@]}" || echo "0")"
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
ASSEMBLER="${ASSEMBLER_ARRAY[$JOB_INDEX]}"
ASSEMBLY_FILE="${ASSEMBLY_FILE_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Sample: ${SAMPLE}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Assembly file: ${ASSEMBLY_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####==================================####
####          CREATE OUTPUT DIR       ####
####==================================####

# Output structure: mapping_tool/assembler/sample
SAMPLE_OUTPUT="${OUTPUT_DIR}/${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
mkdir -p "${SAMPLE_OUTPUT}"

# Change to output directory
cd "$SAMPLE_OUTPUT" || exit 1

####==================================####
####      FUNCTION: RUN VIRSORTER2    ####
####==================================####

run_virsorter2() {
    echo ""
    echo "--- Running VirSorter2 ---"
    
    local output_dir="${SAMPLE_OUTPUT}/virsorter2"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/final-viral-combined.fa" ]] && [[ -s "${output_dir}/final-viral-combined.fa" ]]; then
        echo "  VirSorter2 already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running VirSorter2..."
    
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        "${VIRSORTER2_CONTAINER}" \
        virsorter run \
            -i "/input/$(basename "$ASSEMBLY_FILE")" \
            -o "/output/virsorter2" \
            --include-groups dsDNAphage,ssDNA,ssRNA,dsRNA \
            --min-length 1000 \
            --min-score 0.5 \
            --threads ${THREADS} \
            2>&1 | tee "${output_dir}/virsorter2.log"
    
    if [[ $? -eq 0 ]]; then
        echo "VirSorter2 completed successfully"
        if [[ -f "${output_dir}/final-viral-combined.fa" ]]; then
            viral_count=$(grep -c '^>' "${output_dir}/final-viral-combined.fa")
            echo "  Found ${viral_count} viral contigs"
        fi
        return 0
    else
        echo "VirSorter2 failed"
        return 1
    fi
}

####==================================####
####    FUNCTION: RUN DEEPVIRFINDER   ####
####==================================####

run_deepvirfinder() {
    echo ""
    echo "--- Running DeepVirFinder ---"
    
    local output_dir="${SAMPLE_OUTPUT}/deepvirfinder"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/${SAMPLE}_prediction.txt" ]] && [[ -s "${output_dir}/${SAMPLE}_prediction.txt" ]]; then
        echo "  DeepVirFinder already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running DeepVirFinder..."
    
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        "${DEEPVIRFINDER_CONTAINER}" \
        python /DeepVirFinder/dvf.py \
            -i "/input/$(basename "$ASSEMBLY_FILE")" \
            -o "/output/deepvirfinder" \
            -l 1000 \
            -c ${THREADS} \
            2>&1 | tee "${output_dir}/deepvirfinder.log"
    
    if [[ $? -eq 0 ]]; then
        echo "  DeepVirFinder completed successfully"
        if [[ -f "${output_dir}/${SAMPLE}_prediction.txt" ]]; then
            awk '$3 <= 0.05' "${output_dir}/${SAMPLE}_prediction.txt" > \
                "${output_dir}/viral_contigs.txt"
            viral_count=$(wc -l < "${output_dir}/viral_contigs.txt")
            echo "  Found ${viral_count} viral contigs (p < 0.05)"
        fi
        return 0
    else
        echo "  DeepVirFinder failed"
        return 1
    fi
}

####==================================####
####      FUNCTION: RUN CHECKV       ####
####==================================####

run_checkv() {
    echo ""
    echo "--- Running CheckV ---"
    
    local output_dir="${SAMPLE_OUTPUT}/checkv"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/quality_summary.tsv" ]] && [[ -s "${output_dir}/quality_summary.tsv" ]]; then
        echo "  CheckV already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running CheckV..."
    
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        "${CHECKV_CONTAINER}" \
        checkv end_to_end \
            "/input/$(basename "$ASSEMBLY_FILE")" \
            "/output/checkv" \
            -t ${THREADS} \
            2>&1 | tee "${output_dir}/checkv.log"
    
    if [[ $? -eq 0 ]]; then
        echo " CheckV completed successfully"
        return 0
    else
        echo "  CheckV failed"
        return 1
    fi
}

####==================================####
####      FUNCTION: RUN GRAVITY      ####
####==================================####

run_gravity() {
    echo ""
    echo "--- Running GraViTy ---"
    
    local output_dir="${SAMPLE_OUTPUT}/gravity"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/taxonomy_results.txt" ]] && [[ -s "${output_dir}/taxonomy_results.txt" ]]; then
        echo "  GraViTy already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running GraViTy..."
    
    # First, we need to get the viral contigs from VirSorter2 or DeepVirFinder
    local viral_contigs=""
    if [[ -f "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" ]]; then
        viral_contigs="${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa"
    elif [[ -f "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]]; then
        # Extract sequences from assembly file based on DeepVirFinder contig IDs
        python -c "
import sys
with open('${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt') as f:
    contigs = [line.split()[0] for line in f]
with open('${ASSEMBLY_FILE}') as f:
    for line in f:
        if line.startswith('>') and line.strip()[1:] in contigs:
            sys.stdout.write(line + next(f))
" > "${output_dir}/viral_contigs.fa"
        viral_contigs="${output_dir}/viral_contigs.fa"
    else
        echo "  No viral contigs found for GraViTy analysis"
        return 1
    fi
    
    if [[ ! -s "$viral_contigs" ]]; then
        echo "  No viral contigs found for GraViTy analysis"
        return 1
    fi
    
    apptainer exec \
        --bind "$(dirname "$viral_contigs")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        "${GRAVITY_CONTAINER}" \
        gravy -i "/input/$(basename "$viral_contigs")" \
            -o "/output/gravity" \
            -db /gravity_db \
            -t ${THREADS} \
            2>&1 | tee "${output_dir}/gravity.log"
    
    if [[ $? -eq 0 ]]; then
        echo " GraViTy completed successfully"
        return 0
    else
        echo " GraViTy failed"
        return 1
    fi
}

####==================================####
####      FUNCTION: RUN DIAMOND      ####
####==================================####

run_diamond() {
    echo ""
    echo "--- Running DIAMOND ---"
    
    local output_dir="${SAMPLE_OUTPUT}/diamond"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/diamond_results.tsv" ]] && [[ -s "${output_dir}/diamond_results.tsv" ]]; then
        echo "  DIAMOND already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running DIAMOND..."
    
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        --bind "${CONTAINERS}/diamond_db":/db:ro \
        "${DIAMOND_CONTAINER}" \
        diamond blastx \
            -q "/input/$(basename "$ASSEMBLY_FILE")" \
            -d "/db/viral_nr.dmnd" \
            -o "/output/diamond/diamond_results.tsv" \
            --ultra-sensitive \
            --evalue 1e-5 \
            --threads ${THREADS} \
            --outfmt 6 \
            2>&1 | tee "${output_dir}/diamond.log"
    
    if [[ $? -eq 0 ]]; then
        echo " DIAMOND completed successfully"
        return 0
    else
        echo "  DIAMOND failed"
        return 1
    fi
}

####==================================####
####      FUNCTION: CONSOLIDATE       ####
####==================================####

consolidate_results() {
    echo ""
    echo "--- Consolidating Results ---"
    
    local report="${SAMPLE_OUTPUT}/consolidated_report.txt"
    
    cat > "$report" << EOF
===== VIRAL IDENTIFICATION REPORT =====
Mapping Tool: ${MAPPING_TOOL}
Assembler: ${ASSEMBLER}
Sample: ${SAMPLE}
Assembly file: $(basename "$ASSEMBLY_FILE")
Date: $(date)
Job: ${LSB_JOBINDEX}/${TOTAL_JOBS}

--- ASSEMBLY STATISTICS ---
$(if [[ -f "$ASSEMBLY_FILE" ]]; then
    contigs=$(grep -c '^>' "$ASSEMBLY_FILE" 2>/dev/null || echo "0")
    total_len=$(grep -v '^>' "$ASSEMBLY_FILE" | tr -d '\n' | wc -c 2>/dev/null || echo "0")
    echo "Total contigs: $contigs"
    echo "Total length: $total_len bp"
else
    echo "Assembly file not found"
fi)

--- VirSorter2 Results ---
$(if [[ -f "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" ]]; then
    viral_count=$(grep -c '^>' "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" 2>/dev/null || echo "0")
    echo "Viral contigs found: $viral_count"
    echo "Categories:"
    grep '^>' "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" 2>/dev/null | \
        cut -d'|' -f2 | sort | uniq -c | head -20
else
    echo "VirSorter2 results not available"
fi)

--- DeepVirFinder Results ---
$(if [[ -f "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]]; then
    viral_count=$(wc -l < "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" 2>/dev/null || echo "0")
    echo "Viral contigs found (p < 0.05): $viral_count"
    echo "Top hits:"
    head -10 "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" 2>/dev/null || echo "No hits"
else
    echo "DeepVirFinder results not available"
fi)

--- CheckV Completeness ---
$(if [[ -f "${SAMPLE_OUTPUT}/checkv/quality_summary.tsv" ]]; then
    echo "Completeness summary:"
    echo "Contig\tCompleteness\tContamination\tQuality"
    awk -F'\t' 'NR>1 && NR<=10 {print $1 "\t" $2 "\t" $3 "\t" $8}' \
        "${SAMPLE_OUTPUT}/checkv/quality_summary.tsv" 2>/dev/null || echo "No data"
else
    echo "CheckV results not available"
fi)

--- GraViTy Taxonomy ---
$(if [[ -f "${SAMPLE_OUTPUT}/gravity/taxonomy_results.txt" ]]; then
    echo "Taxonomic assignments:"
    head -20 "${SAMPLE_OUTPUT}/gravity/taxonomy_results.txt" 2>/dev/null || echo "No data"
else
    echo "GraViTy results not available"
fi)

--- DIAMOND Annotation ---
$(if [[ -f "${SAMPLE_OUTPUT}/diamond/diamond_results.tsv" ]]; then
    echo "Top viral hits:"
    echo "Query\tHit\tE-value\tBitscore"
    awk -F'\t' 'NR<=10 {print $1 "\t" $2 "\t" $11 "\t" $12}' \
        "${SAMPLE_OUTPUT}/diamond/diamond_results.tsv" 2>/dev/null || echo "No data"
else
    echo "DIAMOND results not available"
fi)

=======================================
EOF
    
    echo "  Report created: $report"
}

####==================================####
####          RUN ALL TOOLS          ####
####==================================####

echo ""
echo "========================================="
echo "  Starting viral identification tools"
echo "  ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
echo "========================================="

# Run each tool
run_virsorter2
run_deepvirfinder
run_checkv
run_gravity
run_diamond

# Consolidate results
consolidate_results

####==================================####
####          COMPLETION              ####
####==================================####

echo ""
echo "========================================="
echo "  Completed job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${SAMPLE_OUTPUT}"
echo "  Report: ${SAMPLE_OUTPUT}/consolidated_report.txt"
echo "========================================="