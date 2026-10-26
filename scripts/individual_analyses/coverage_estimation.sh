#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script estimate vOTU coverage with CoverM.
# 14/10/2026
# Ver. 2.0.0 

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input: BAMs (flat per-sample, no assembler level)
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/reads_vs_vOTU"
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/STAR"
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/Bowtie"

# vOTU reference tree: <mapper>/<read_type>/<sample>/vOTUs.fa
VOTU_BASE="${PROJECT_SCRATCH}/results/vOTU"

# Output
OUTPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_DIR_COVER="${OUTPUT_DIR_BASE}/coverage"

PERMANENT_RESULTS="${PERMANENT_BASE}/git_repos/${PROJECT_NAME}/results/genome_coverage"

CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
COVERM_CONTAINER="${CONTAINERS}/coverm_0.8.0--h750ce8b_0.sif"

SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=24

echo "========================================="
echo "  Viral Genome Coverage Estimation (per sample)"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  CoverM container: ${COVERM_CONTAINER}"
echo "  BAM input:       ${INPUT_DIR_BASE}/<mapper>/<read_type>/<SAMPLE>_vOTU_sorted.bam"
echo "  vOTU reference:  ${VOTU_BASE}/<mapper>/<read_type>/<sample>/vOTUs.fa"
echo "  Output (scratch): ${OUTPUT_DIR_COVER}"
echo "  Output (NFS):     ${PERMANENT_RESULTS}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting coverage estimation pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE LIST            ####
####==================================####

declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a READ_TYPE_ARRAY
declare -a MAP_FILE_ARRAY

add_combination() {
    local mapping_tool=$1
    local read_type=$2
    local sample=$3
    local map_file=$4

    if [[ -f "$map_file" ]] && [[ -s "$map_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        MAP_FILE_ARRAY+=("$map_file")
        echo "  Added: $mapping_tool/$read_type/$sample - $(basename "$map_file")"
    else
        echo "  WARNING: BAM file not found or empty: $map_file"
    fi
}

process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2

    [[ -d "$base_dir" ]] || { echo "  WARNING: Directory not found: $base_dir"; return; }

    echo "  Processing mapping tool: $mapping_tool"

    # Two levels: read_type -> BAM files (flat)
    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type
        read_type=$(basename "$read_type_dir")

        shopt -s nullglob
        for bam in "$read_type_dir"/*_vOTU_sorted.bam; do
            local sample
            sample=$(basename "$bam" _vOTU_sorted.bam)
            add_combination "$mapping_tool" "$read_type" "$sample" "$bam"
        done
        shopt -u nullglob
    done
}

process_mapping_tool "STAR"   "$STAR_ASSEMBLY_DIR"
process_mapping_tool "Bowtie" "$BOWTIE_ASSEMBLY_DIR"

TOTAL_JOBS=${#SAMPLE_ARRAY[@]}

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No valid BAM files found"
    tg_send "ERROR: No BAM files found for coverage" 2>/dev/null || true
    exit 1
fi

echo ""
echo "========================================="
echo "  Total jobs: ${TOTAL_JOBS}"
echo "  Unique samples: $(printf '%s\n' "${SAMPLE_ARRAY[@]}" | sort -u | wc -l)"
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
BAM_FILE="${MAP_FILE_ARRAY[$JOB_INDEX]}"

VOTU_FASTA="${VOTU_BASE}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}/vOTUs.fa"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Read type:    ${READ_TYPE}"
echo "  Sample:       ${SAMPLE}"
echo "  BAM file:     ${BAM_FILE}"
echo "  vOTU FASTA:   ${VOTU_FASTA}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####==================================####
####          COVERM FUNCTION         ####
####==================================####

SAMPLE_OUTPUT_COVERM="${OUTPUT_DIR_COVER}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"

run_coverm() {
    echo ""
    echo "--- Running CoverM ---"
    local output_dir="${SAMPLE_OUTPUT_COVERM}"
    mkdir -p "$output_dir"

    local out_tsv="${output_dir}/coverm_abundance_${SAMPLE}.tsv"
    if [[ -s "${out_tsv}" ]]; then
        echo "CoverM already completed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
        return 0
    fi

    [[ -f "${COVERM_CONTAINER}" ]] || { echo "  ERROR: container not found"; return 1; }
    [[ -s "${VOTU_FASTA}" ]] || { echo "  ERROR: vOTU FASTA missing/empty: ${VOTU_FASTA}"; return 1; }

    # Split multi-FASTA into one file per reference (CoverM's genome mode expects this)
    local split_dir="${output_dir}/vOTU_split"
    local n_refs
    n_refs=$(grep -c '^>' "${VOTU_FASTA}")
    local n_split
    n_split=$(find "${split_dir}" -maxdepth 1 -name '*.fna' 2>/dev/null | wc -l)

    if [[ ! -d "${split_dir}" ]] || [[ "${n_split}" -ne "${n_refs}" ]]; then
        rm -rf "${split_dir}"
        mkdir -p "${split_dir}"
        awk -v dir="${split_dir}" '
            /^>/ {
                if (out) close(out)
                hdr = substr($0, 2)
                sub(/[ \t].*/, "", hdr)
                gsub(/[^A-Za-z0-9_.-]/, "_", hdr)
                out = dir "/" hdr ".fna"
            }
            { print > out }
        ' "${VOTU_FASTA}"
    fi

    apptainer exec \
        --bind "$(dirname "${BAM_FILE}")":/input:ro \
        --bind "${split_dir}":/ref:ro \
        --bind "${output_dir}":/output \
        "${COVERM_CONTAINER}" \
        coverm genome \
            --bam-files "/input/$(basename "${BAM_FILE}")" \
            --genome-fasta-directory /ref \
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
}

####==================================####
####         RUN + COPY TO NFS        ####
####==================================####

SAMPLE_OUTPUT_NFS="${PERMANENT_RESULTS}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"

if run_coverm; then
    if [[ -d "${SAMPLE_OUTPUT_COVERM}" ]]; then
        echo "Copying results to permanent storage..."
        mkdir -p "${SAMPLE_OUTPUT_NFS}"

        cp -v "${SAMPLE_OUTPUT_COVERM}/coverm_abundance_${SAMPLE}.tsv" \
              "${SAMPLE_OUTPUT_NFS}/"
        if [[ -f "${SAMPLE_OUTPUT_COVERM}/coverm.log" ]]; then
            cp -v "${SAMPLE_OUTPUT_COVERM}/coverm.log" \
                  "${SAMPLE_OUTPUT_NFS}/"
        fi

        echo "Results copied to: ${SAMPLE_OUTPUT_NFS}"
        tg_send "Results copied: ${SAMPLE_OUTPUT_NFS}" 2>/dev/null || true
    fi
    exit 0
else
    echo "Job ${LSB_JOBINDEX}/${TOTAL_JOBS} FAILED for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
    tg_send "FAILED: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi
