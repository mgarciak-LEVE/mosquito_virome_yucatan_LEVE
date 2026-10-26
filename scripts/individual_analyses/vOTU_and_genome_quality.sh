#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script assign viral OTUs and assess genome quality. CD-HIT removed (upstream concatenation already dereplicated)
# 14/10/2026
# Ver. 2.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input: the merged contig_concatenation tree (no assembler level)
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/contig_concatenation"
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/STAR"
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/Bowtie"

# Output: also no assembler level
OUTPUT_DIR_OTU="${PROJECT_SCRATCH}/results/vOTU"
OUTPUT_DIR_COMPLETENESS="${PROJECT_SCRATCH}/results/completeness"

# Containers
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
CHECKV_CONTAINER="${CONTAINERS}/checkv.sif"

# Databases
DB_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/databases"
CHECKV_DB="${DB_BASE}/checkv_db/checkv-db-v1.5"

# Scripts
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=12

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  vOTU + Completeness (per sample)"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  Input (merged):"
echo "    - STAR:   ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output for vOTUs:        ${OUTPUT_DIR_OTU}"
echo "  Output for completeness: ${OUTPUT_DIR_COMPLETENESS}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting vOTU + completeness pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE LIST            ####
####==================================####

# One entry per (mapper, read_type, sample). Two levels of walk,
declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a READ_TYPE_ARRAY
declare -a CONTIG_FILE_ARRAY

add_combination() {
    local mapping_tool=$1
    local read_type=$2
    local sample=$3
    local contig_file=$4

    if [[ -f "$contig_file" ]] && [[ -s "$contig_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        CONTIG_FILE_ARRAY+=("$contig_file")
        echo "  Added: $mapping_tool/$read_type/$sample - $(basename "$contig_file")"
    else
        echo "  WARNING: contig file missing or empty: $contig_file"
    fi
}

process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2

    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type; read_type=$(basename "$read_type_dir")

        for sample_dir in "$read_type_dir"/*/; do
            [[ -d "$sample_dir" ]] || continue
            local sample; sample=$(basename "$sample_dir")
            local contig_file="${sample_dir}/all.fa"
            add_combination "$mapping_tool" "$read_type" "$sample" "$contig_file"
        done
    done
}

process_mapping_tool "STAR"   "$STAR_ASSEMBLY_DIR"
process_mapping_tool "Bowtie" "$BOWTIE_ASSEMBLY_DIR"

TOTAL_JOBS=${#SAMPLE_ARRAY[@]}

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No valid input files found"
    tg_send "ERROR: No valid contig files found for vOTU step" 2>/dev/null || true
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
CONTIG_FILE="${CONTIG_FILE_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Read Type:    ${READ_TYPE}"
echo "  Sample:       ${SAMPLE}"
echo "  Contig file:  ${CONTIG_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####=========================================####
####   vOTUs = copy of merged all.fa (no CD-HIT)
####=========================================####

SAMPLE_OUTPUT_OTU="${OUTPUT_DIR_OTU}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"

run_votu_copy() {
    echo ""
    echo "--- Copying merged contigs as vOTUs.fa (no CD-HIT) ---"
    local output_dir="${SAMPLE_OUTPUT_OTU}"
    mkdir -p "$output_dir"

    local out_file="${output_dir}/vOTUs.fa"

    if [[ -s "$out_file" ]]; then
        echo "  vOTUs.fa already present for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
        return 0
    fi

    if [[ ! -s "$CONTIG_FILE" ]]; then
        echo "  ERROR: input all.fa missing or empty: $CONTIG_FILE"
        return 1
    fi

    cp "$CONTIG_FILE" "$out_file"

    if [[ ! -s "$out_file" ]]; then
        echo "  ERROR: copy produced an empty vOTUs.fa"
        return 1
    fi

    local n
    n=$(grep -c '^>' "$out_file" || true)
    echo "  vOTUs.fa written: ${n} contigs"
    return 0
}

####============================================####
####   COMPLETENESS OUTPUT DIR
####============================================####

SAMPLE_OUTPUT_COMPLETENESS="${OUTPUT_DIR_COMPLETENESS}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"

run_checkv() {
    echo ""
    echo "--- Running CheckV ---"
    local output_dir="${SAMPLE_OUTPUT_COMPLETENESS}"
    mkdir -p "$output_dir"

    if [[ -f "${output_dir}/quality_summary.tsv" ]]; then
        echo "  CheckV already completed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
        return 0
    fi

    if [[ ! -d "${CHECKV_DB}" ]]; then
        echo "  WARNING: CheckV database not found at ${CHECKV_DB}"
        return 0
    fi

    if [[ ! -f "${CHECKV_CONTAINER}" ]]; then
        echo "  ERROR: Container not found at ${CHECKV_CONTAINER}"
        return 1
    fi

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
        echo "  CheckV failed (exit ${checkv_exit})"
        return 1
    fi

    return 0
}

####==================================####
####          RUN                  ####
####==================================####

echo ""
echo "========================================="
echo "  Starting vOTU + completeness"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
echo "========================================="

if ! run_votu_copy; then
    echo "ERROR: vOTU copy failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
    tg_send "ERROR: vOTU copy failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

if ! run_checkv; then
    echo "ERROR: CheckV failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
    tg_send "ERROR: CheckV failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
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
echo "  Output: ${SAMPLE_OUTPUT_OTU} & ${SAMPLE_OUTPUT_COMPLETENESS}"
echo "========================================="

tg_send "Completed: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true
