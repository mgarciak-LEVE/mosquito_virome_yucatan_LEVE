#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to index vOTU fasta files with Bowtie2.
# 06/10/2026
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
OUTPUT_DIR="${OUTPUT_DIR_BASE}/indexed_vOTU"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
BOWTIE_CONTAINER="${CONTAINERS}/bowtie2_2.5.5.sif"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

#==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=36

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  vOTU Indexing Pipeline"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  Bowtie2 container: ${BOWTIE_CONTAINER}"
echo "  Input directories (mapping_tool/assembler/read_type/sample):"
echo "    - STAR: ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output for Bowtie indexing: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting vOTU indexing pipeline" 2>/dev/null || true

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

####==================================####
####          BOWTIE INDEXING         ####
####==================================####

# Directory for every sample
# Output structure: mapping_tool/read_type/assembler/sample
SAMPLE_OUTPUT="${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"

run_bowtie() {
    echo ""
    echo "--- Running Bowtie ---"
    local output_dir="${SAMPLE_OUTPUT}"
    mkdir -p "$output_dir"
    
    # Check if it has been completed already 
    if [[ -f "${output_dir}/vOTUs.1.bt2" ]] && [[ -s "${output_dir}/vOTUs.1.bt2" ]]; then 
        echo "Bowtie indexing already completed for ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi


    # Verify the Bowtie2 container exists
    if [[ ! -f "${BOWTIE_CONTAINER}" ]]; then
        echo "  ERROR: Container not found at ${BOWTIE_CONTAINER}"
        return 1
    fi

    # Run bowtie2 using apptainer container
    apptainer exec \
        --bind "$(dirname "$OTU_FILE")":/input:ro \
        --bind "${output_dir}":/output \
        "${BOWTIE_CONTAINER}" \
        bowtie2-build \
        /input/vOTUs.fa \
        /output/vOTUs \
        --threads "${THREADS}" \
        2>&1 | tee "${output_dir}/bowtie.log"

    local bowtie_exit=${PIPESTATUS[0]}
    
    if [[ ${bowtie_exit} -ne 0 ]]; then
        echo "  Bowtie2 failed (exit ${bowtie_exit})"
        return 1
    fi

    return 0

}


####==================================####
####          RUN ALL TOOLS          ####
####==================================####

echo ""
echo "========================================="
echo "  Starting viral vOTU indexing"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
echo "========================================="

if ! run_bowtie; then
    echo "ERROR: Bowtie failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
    tg_send "ERROR: Bowtie failed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true
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
echo "  Output: ${SAMPLE_OUTPUT}"
echo "============="
