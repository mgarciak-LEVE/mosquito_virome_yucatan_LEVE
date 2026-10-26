#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to index vOTU fasta files with Bowtie2.
# 14/10/2026
# Ver. 2.0.0

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input: merged vOTU tree (no assembler)
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/vOTU"
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/STAR"
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/Bowtie"

# Output: same shape, no assembler
OUTPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/indexed_vOTU"

CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
BOWTIE_CONTAINER="${CONTAINERS}/bowtie2_2.5.5.sif"

SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=36

echo "========================================="
echo "  vOTU Indexing Pipeline (per sample)"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  Bowtie2 container: ${BOWTIE_CONTAINER}"
echo "  Input (merged vOTUs):"
echo "    - STAR:   ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output:     ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting vOTU indexing pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE LIST            ####
####==================================####

declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a READ_TYPE_ARRAY
declare -a OTU_FILE_ARRAY

add_combination() {
    local mapping_tool=$1
    local read_type=$2
    local sample=$3
    local otu_file=$4

    if [[ -f "$otu_file" ]] && [[ -s "$otu_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        OTU_FILE_ARRAY+=("$otu_file")
        echo "  Added: $mapping_tool/$read_type/$sample - $(basename "$otu_file")"
    else
        echo "  WARNING: vOTU fasta missing or empty: $otu_file"
    fi
}

process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2

    [[ -d "$base_dir" ]] || { echo "  WARNING: Directory not found: $base_dir"; return; }

    echo "  Processing mapping tool: $mapping_tool"

    # Two levels: read_type -> sample
    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type
        read_type=$(basename "$read_type_dir")

        for sample_dir in "$read_type_dir"/*/; do
            [[ -d "$sample_dir" ]] || continue
            local sample
            sample=$(basename "$sample_dir")
            local otu_file="${sample_dir}/vOTUs.fa"
            add_combination "$mapping_tool" "$read_type" "$sample" "$otu_file"
        done
    done
}

process_mapping_tool "STAR"   "$STAR_ASSEMBLY_DIR"
process_mapping_tool "Bowtie" "$BOWTIE_ASSEMBLY_DIR"

TOTAL_JOBS=${#SAMPLE_ARRAY[@]}

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No valid vOTU files found"
    tg_send "ERROR: No vOTU files found for indexing" 2>/dev/null || true
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

JOB_INDEX=$((LSB_JOBINDEX - 1))

if [[ $JOB_INDEX -ge ${#SAMPLE_ARRAY[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#SAMPLE_ARRAY[@]})"
    exit 1
fi

SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"
MAPPING_TOOL="${MAPPING_ARRAY[$JOB_INDEX]}"
READ_TYPE="${READ_TYPE_ARRAY[$JOB_INDEX]}"
OTU_FILE="${OTU_FILE_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Read Type:    ${READ_TYPE}"
echo "  Sample:       ${SAMPLE}"
echo "  vOTU file:    ${OTU_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Indexing: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####==================================####
####          BOWTIE INDEXING         ####
####==================================####

SAMPLE_OUTPUT="${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"

run_bowtie() {
    echo ""
    echo "--- Running Bowtie2-build ---"
    local output_dir="${SAMPLE_OUTPUT}"
    mkdir -p "$output_dir"

    if [[ -f "${output_dir}/vOTUs.1.bt2" ]] && [[ -s "${output_dir}/vOTUs.1.bt2" ]]; then
        echo "Bowtie indexing already completed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
        return 0
    fi

    if [[ ! -f "${BOWTIE_CONTAINER}" ]]; then
        echo "  ERROR: Container not found at ${BOWTIE_CONTAINER}"
        return 1
    fi

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
        echo "  Bowtie2-build failed (exit ${bowtie_exit})"
        return 1
    fi

    return 0
}

####==================================####
####          RUN                     ####
####==================================####

echo ""
echo "========================================="
echo "  Starting vOTU indexing"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
echo "========================================="

if ! run_bowtie; then
    echo "ERROR: Bowtie indexing failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
    tg_send "ERROR: Bowtie indexing failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

####==================================####
####          COMPLETION              ####
####==================================####

echo ""
echo "========================================="
echo "  Completed job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Read Type:    ${READ_TYPE}"
echo "  Sample:       ${SAMPLE}"
echo "  Output:       ${SAMPLE_OUTPUT}"
echo "========================================="

tg_send "Completed indexing: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true
