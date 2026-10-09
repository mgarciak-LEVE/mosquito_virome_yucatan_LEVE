#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script estimate vOTU coverage with CoverM.
# 08/10/2026
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

# Input directories (BAMs, flat layout: tool/assembler/read_type/<SAMPLE>_vOTU_sorted.bam)
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/reads_vs_vOTU"
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/STAR"
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/Bowtie"

# vOTU reference tree (per-sample FASTAs: tool/read_type/assembler/<SAMPLE_DIR>/vOTUs.fa)
VOTU_BASE="${PROJECT_SCRATCH}/results/vOTU"

# Output directories (scratch)
OUTPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_DIR_COVER="${OUTPUT_DIR_BASE}/coverage"

# Permanent results directory (NFS) — mirror of scratch coverage tree
PERMANENT_RESULTS="${PERMANENT_BASE}/git_repos/${PROJECT_NAME}/results/genome_coverage"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
COVERM_CONTAINER="${CONTAINERS}/coverm_0.8.0--h750ce8b_0.sif"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=24

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Genome Coverage Estimation Pipeline"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  CoverM container: ${COVERM_CONTAINER}"
echo "  Input directories (mapping_tool/assembler/read_type/<sample>_vOTU_sorted.bam):"
echo "    - STAR:   ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  vOTU reference base:  ${VOTU_BASE}"
echo "  Output (scratch):     ${OUTPUT_DIR_COVER}"
echo "  Output (NFS):         ${PERMANENT_RESULTS}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting viral genomes coverage estimation pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE-ASSEMBLER LIST  ####
####==================================####

declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a READ_TYPE_ARRAY
declare -a ASSEMBLER_ARRAY
declare -a MAP_FILE_ARRAY

echo "Building sample-assembler list..."

add_combination() {
    local mapping_tool=$1
    local read_type=$2
    local assembler=$3
    local sample=$4
    local map_file=$5

    if [[ -f "$map_file" ]] && [[ -s "$map_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        ASSEMBLER_ARRAY+=("$assembler")
        MAP_FILE_ARRAY+=("$map_file")
        echo "  Added: $mapping_tool/$assembler/$read_type/$sample - $(basename "$map_file")"
    else
        echo "  WARNING: BAM file not found or empty for $mapping_tool/$assembler/$read_type/$sample: $map_file"
    fi
}

process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2

    [[ -d "$base_dir" ]] || { echo "  WARNING: Directory not found: $base_dir"; return; }

    echo "  Processing mapping tool: $mapping_tool"

    for assembler_dir in "$base_dir"/*/; do
        [[ -d "$assembler_dir" ]] || continue
        assembler=$(basename "$assembler_dir")

        for read_type_dir in "$assembler_dir"/*/; do
            [[ -d "$read_type_dir" ]] || continue
            read_type=$(basename "$read_type_dir")

            shopt -s nullglob
            for bam in "$read_type_dir"/*_vOTU_sorted.bam; do
                sample=$(basename "$bam" _vOTU_sorted.bam)
                add_combination "$mapping_tool" "$read_type" "$assembler" "$sample" "$bam"
            done
            shopt -u nullglob
        done
    done
}

process_mapping_tool "STAR"   "$STAR_ASSEMBLY_DIR"
process_mapping_tool "Bowtie" "$BOWTIE_ASSEMBLY_DIR"

TOTAL_JOBS=${#SAMPLE_ARRAY[@]}

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No valid BAM files found"
    tg_send "ERROR: No valid BAM files found for coverage estimation" 2>/dev/null || true
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
ASSEMBLER="${ASSEMBLER_ARRAY[$JOB_INDEX]}"
BAM_FILE="${MAP_FILE_ARRAY[$JOB_INDEX]}"

# vOTU reference tree layout: <tool>/<read_type>/<assembler>/<SAMPLE_DIR>/vOTUs.fa
# paired:   SAMPLE_DIR == SAMPLE
# unpaired: SAMPLE_DIR == SAMPLE + "_R{1,2}_unpaired"
case "${READ_TYPE}" in
    paired)      VOTU_SAMPLE_DIR="${SAMPLE}" ;;
    unpaired_R1) VOTU_SAMPLE_DIR="${SAMPLE}_R1_unpaired" ;;
    unpaired_R2) VOTU_SAMPLE_DIR="${SAMPLE}_R2_unpaired" ;;
    *)           VOTU_SAMPLE_DIR="${SAMPLE}" ;;
esac

VOTU_FASTA="${VOTU_BASE}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${VOTU_SAMPLE_DIR}/vOTUs.fa"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Assembler:    ${ASSEMBLER}"
echo "  Read type:    ${READ_TYPE}"
echo "  Sample:       ${SAMPLE}"
echo "  BAM file:     ${BAM_FILE}"
echo "  vOTU FASTA:   ${VOTU_FASTA}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####==================================####
####          COVERM FUNCTION         ####
####==================================####

SAMPLE_OUTPUT_COVERM="${OUTPUT_DIR_COVER}/${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE}"

run_coverm() {
    echo ""
    echo "--- Running CoverM ---"
    local output_dir="${SAMPLE_OUTPUT_COVERM}"
    mkdir -p "$output_dir"

    local out_tsv="${output_dir}/coverm_abundance_${SAMPLE}.tsv"
    if [[ -s "${out_tsv}" ]]; then
        echo "CoverM already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE}"
        return 0
    fi

    if [[ ! -f "${COVERM_CONTAINER}" ]]; then
        echo "  ERROR: CoverM container not found at ${COVERM_CONTAINER}"
        return 1
    fi

    if [[ ! -f "${VOTU_FASTA}" ]]; then
        echo "  ERROR: vOTU reference FASTA not found: ${VOTU_FASTA}"
        tg_send "ERROR: missing vOTU FASTA for ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
        return 1
    fi

    if [[ ! -s "${VOTU_FASTA}" ]]; then
        echo "  ERROR: vOTU reference FASTA is empty: ${VOTU_FASTA}"
        return 1
    fi

    apptainer exec \
        --bind "$(dirname "${BAM_FILE}")":/input:ro \
        --bind "$(dirname "${VOTU_FASTA}")":/ref:ro \
        --bind "${output_dir}":/output \
        "${COVERM_CONTAINER}" \
        coverm genome \
        --bam-files "/input/${SAMPLE}_vOTU_sorted.bam" \
        --genome-fasta-files "/ref/$(basename "${VOTU_FASTA}")" \
        -m mean relative_abundance covered_fraction \
        --min-covered-fraction 0 \
        -t "${THREADS}" \
        --output-file "/output/coverm_abundance_${SAMPLE}.tsv" \
        2>&1 | tee "${output_dir}/coverm.log"

    local coverm_exit=${PIPESTATUS[0]}
    if [[ ${coverm_exit} -ne 0 ]]; then
        echo "  CoverM failed (exit ${coverm_exit})"
        return 1
    fi
    return 0
}

####==================================####
####         RUN + COPY TO NFS        ####
####==================================####

SAMPLE_OUTPUT_NFS="${PERMANENT_RESULTS}/${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE}"

if run_coverm; then
    if [[ -d "${SAMPLE_OUTPUT_COVERM}" ]]; then
        echo "Copying results to permanent storage..."
        mkdir -p "${SAMPLE_OUTPUT_NFS}"
        cp -r "${SAMPLE_OUTPUT_COVERM}/." "${SAMPLE_OUTPUT_NFS}/"
        echo "Results copied to: ${SAMPLE_OUTPUT_NFS}"
        tg_send "Results copied to NFS: ${SAMPLE_OUTPUT_NFS}" 2>/dev/null || true
    else
        echo "WARNING: No output directory found to copy to permanent storage: ${SAMPLE_OUTPUT_COVERM}"
        tg_send "WARNING: No CoverM output to copy to NFS for ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
    fi
    echo ""
    echo "Job ${LSB_JOBINDEX}/${TOTAL_JOBS} finished for ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE}"
    tg_send "Done: ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true
    exit 0
else
    echo ""
    echo "Job ${LSB_JOBINDEX}/${TOTAL_JOBS} FAILED for ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE}"
    tg_send "FAILED: ${MAPPING_TOOL}/${ASSEMBLER}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true
    exit 1
fi
