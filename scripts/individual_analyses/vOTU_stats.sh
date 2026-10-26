#!/usr/bin/env bash
# Author: Jorge Alberto Castro Rodríguez (adapted)
# Script: per-OTU stats (seqkit) + CheckV join + NFS copy
# Version: 1.0.0
# Date: 22/10/2026

####==================================####
####           CONFIGURATION          ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
PROJECT_NFS="${PERMANENT_BASE}/git_repos/${PROJECT_NAME}"

# Inputs (outputs of the vOTU + completeness pipeline)
INPUT_DIR_OTU="${PROJECT_SCRATCH}/results/vOTU"
INPUT_DIR_COMPLETENESS="${PROJECT_SCRATCH}/results/completeness"

# Scratch output for stats
OUTPUT_DIR_STATS="${PROJECT_SCRATCH}/docs/vOTU_stats"

# NFS output (mirrored layout)
NFS_OUTPUT_DIR_STATS="${PROJECT_NFS}/docs/vOTU_stats"

# Containers (only used if seqkit is not available natively)
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
SEQKIT_CONTAINER="${CONTAINERS}/seqkit_2.8.2.sif"

# Scripts
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  vOTU per-OTU stats + NFS copy"
echo "  Project: ${PROJECT_NAME}"
echo "  Input (vOTUs):  ${INPUT_DIR_OTU}"
echo "  Input (CheckV): ${INPUT_DIR_COMPLETENESS}"
echo "  Output (scratch): ${OUTPUT_DIR_STATS}"
echo "  Output (NFS):     ${NFS_OUTPUT_DIR_STATS}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting vOTU stats pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE LIST            ####
####==================================####

# Walk the vOTU output tree (vOTUs.fa is the source of truth for "done").
# One entry per (mapper, read_type, sample) with a non-empty vOTUs.fa.

declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a READ_TYPE_ARRAY
declare -a VOTU_FILE_ARRAY
declare -a CHECKV_DIR_ARRAY

add_combination() {
    local mapping_tool=$1
    local read_type=$2
    local sample=$3
    local votu_file=$4
    local checkv_dir=$5

    if [[ -s "$votu_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        VOTU_FILE_ARRAY+=("$votu_file")
        CHECKV_DIR_ARRAY+=("$checkv_dir")
        echo "  Added: $mapping_tool/$read_type/$sample"
    else
        echo "  SKIP (no vOTUs.fa): $mapping_tool/$read_type/$sample"
    fi
}

process_mapping_tool() {
    local mapping_tool=$1
    local otu_base=$2
    local completeness_base=$3

    if [[ ! -d "$otu_base" ]]; then
        echo "  WARNING: OTU directory not found: $otu_base"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    for read_type_dir in "$otu_base"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type; read_type=$(basename "$read_type_dir")

        for sample_dir in "$read_type_dir"/*/; do
            [[ -d "$sample_dir" ]] || continue
            local sample; sample=$(basename "$sample_dir")
            local votu_file="${sample_dir}/vOTUs.fa"
            local checkv_dir="${completeness_base}/${mapping_tool}/${read_type}/${sample}"
            add_combination "$mapping_tool" "$read_type" "$sample" "$votu_file" "$checkv_dir"
        done
    done
}

process_mapping_tool "STAR"   "${INPUT_DIR_OTU}/STAR"   "${INPUT_DIR_COMPLETENESS}/STAR"
process_mapping_tool "Bowtie" "${INPUT_DIR_OTU}/Bowtie" "${INPUT_DIR_COMPLETENESS}/Bowtie"

TOTAL_JOBS=${#SAMPLE_ARRAY[@]}

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No vOTUs.fa files found in ${INPUT_DIR_OTU}"
    tg_send "ERROR: No vOTUs.fa files found for stats step" 2>/dev/null || true
    exit 1
fi

echo ""
echo "========================================="
echo "  Total jobs: ${TOTAL_JOBS}"
echo "========================================="

####==================================####
####          GET ARRAY TASK          ####
####==================================####

if [[ -z "${LSB_JOBINDEX:-}" ]]; then
    echo "ERROR: LSB_JOBINDEX not set."
    exit 1
fi

JOB_INDEX=$((LSB_JOBINDEX - 1))

if [[ $JOB_INDEX -ge ${#SAMPLE_ARRAY[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#SAMPLE_ARRAY[@]})"
    exit 1
fi

SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"
MAPPING_TOOL="${MAPPING_ARRAY[$JOB_INDEX]}"
READ_TYPE="${READ_TYPE_ARRAY[$JOB_INDEX]}"
VOTU_FILE="${VOTU_FILE_ARRAY[$JOB_INDEX]}"
CHECKV_DIR="${CHECKV_DIR_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
echo "  vOTUs.fa:   ${VOTU_FILE}"
echo "  CheckV dir: ${CHECKV_DIR}"
echo "========================================="

tg_send "Stats job: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####==================================####
####       OUTPUT DIRS              ####
####==================================####

SAMPLE_OUTPUT_STATS="${OUTPUT_DIR_STATS}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
NFS_SAMPLE_OUTPUT_STATS="${NFS_OUTPUT_DIR_STATS}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"

mkdir -p "${SAMPLE_OUTPUT_STATS}"
mkdir -p "${NFS_SAMPLE_OUTPUT_STATS}"

####==================================####
####     PER-OTU STATS              ####
####==================================####

run_per_otu_stats() {
    echo ""
    echo "--- Per-OTU stats (seqkit fx2tab) ---"
    local out_file="${SAMPLE_OUTPUT_STATS}/vOTU_stats.tsv"

    if [[ -s "$out_file" ]]; then
        echo "  Already present: $out_file"
        return 0
    fi

    if [[ ! -s "$VOTU_FILE" ]]; then
        echo "  ERROR: input missing or empty: $VOTU_FILE"
        return 1
    fi

    local votu_dir
    votu_dir=$(dirname "$VOTU_FILE")
    
    # -n : name only (drop seq/qual columns)
    # -l : add length
    # -g : add GC content
    # -i : use the sequence ID (first token) instead of the full header
    # -H : header row
    # Result columns: #name  length  GC
    apptainer exec \
        --bind "${votu_dir}:${votu_dir}:ro" \
        "${SEQKIT_CONTAINER}" \
        seqkit fx2tab -n -l -g -i -H "$VOTU_FILE" > "$out_file"

    local rc=$?
    if [[ $rc -ne 0 ]] || [[ ! -s "$out_file" ]]; then
        echo "  ERROR: seqkit fx2tab failed (rc=${rc})"
        return 1
    fi

    local n
    n=$(( $(wc -l < "$out_file") - 1 ))
    echo "  Wrote ${n} per-OTU rows to $out_file"
    return 0
}

####==================================####
####     WHOLE-FILE SUMMARY          ####
####==================================####

run_overall_stats() {
    echo ""
    echo "--- Whole-file summary (seqkit stats) ---"
    local out_file="${SAMPLE_OUTPUT_STATS}/vOTUs_stats_overall.tsv"

    if [[ -s "$out_file" ]]; then
        echo "  Already present: $out_file"
        return 0
    fi

    local votu_dir
    votu_dir=$(dirname "$VOTU_FILE")

    # -a : all stats (incl. N50)
    # -T : tabular
    apptainer exec \
        --bind "${votu_dir}:${votu_dir}:ro" \
        "${SEQKIT_CONTAINER}" \
        seqkit stats -a -T "$VOTU_FILE" > "$out_file"

    local rc=$?
    if [[ $rc -ne 0 ]] || [[ ! -s "$out_file" ]]; then
        echo "  ERROR: seqkit stats failed (rc=${rc})"
        return 1
    fi

    echo "  Wrote overall summary to $out_file"
    return 0
}

####==================================####
####     JOIN WITH CHECKV            ####
####==================================####

run_join_checkv() {
    echo ""
    echo "--- Join with CheckV quality_summary.tsv ---"
    local stats_file="${SAMPLE_OUTPUT_STATS}/vOTU_stats.tsv"
    local checkv_file="${CHECKV_DIR}/quality_summary.tsv"
    local out_file="${SAMPLE_OUTPUT_STATS}/vOTU_master_stats.tsv"

    if [[ -s "$out_file" ]]; then
        echo "  Already present: $out_file"
        return 0
    fi

    if [[ ! -f "$checkv_file" ]]; then
        echo "  WARNING: CheckV output missing ($checkv_file); writing stats only"
        cp "$stats_file" "$out_file"
        return 0
    fi

    # Header-aware join on column 1 (contig_id).
    # - Drops the CheckV contig_id column from the header (avoids duplicate ID col)
    # - Pads with NA for OTUs missing from CheckV
    awk -F'\t' -v OFS='\t' '
        # --- First file: CheckV quality_summary.tsv ---
        FNR == NR {
            if (FNR == 1) {
                cv_header = $0
                sub(/^[^\t]*\t/, "", cv_header)   # drop contig_id header
                cv_ncols = NF - 1
                next
            }
            checkv[$1] = $0
            next
        }
        # --- Second file: seqkit stats ---
        FNR == 1 {
            print $0, cv_header
            next
        }
        {
            if ($1 in checkv) {
                cline = checkv[$1]
                sub(/^[^\t]*\t/, "", cline)       # drop duplicated contig_id
                print $0, cline
            } else {
                pad = ""
                for (i = 1; i <= cv_ncols; i++) pad = pad "\tNA"
                print $0, pad
            }
        }
    ' "$checkv_file" "$stats_file" > "$out_file"

    echo "  Wrote master table: $out_file"
    return 0
}

####==================================####
####     COPY TO NFS                ####
####==================================####

copy_to_nfs() {
    echo ""
    echo "--- Copy to NFS ---"

    local files=()
    [[ -s "${SAMPLE_OUTPUT_STATS}/vOTU_stats.tsv"         ]] && files+=("${SAMPLE_OUTPUT_STATS}/vOTU_stats.tsv")
    [[ -s "${SAMPLE_OUTPUT_STATS}/vOTU_master_stats.tsv"  ]] && files+=("${SAMPLE_OUTPUT_STATS}/vOTU_master_stats.tsv")
    [[ -s "${SAMPLE_OUTPUT_STATS}/vOTUs_stats_overall.tsv" ]] && files+=("${SAMPLE_OUTPUT_STATS}/vOTUs_stats_overall.tsv")

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "  WARNING: no stats files to copy"
        return 0
    fi

    rsync -a "${files[@]}" "${NFS_SAMPLE_OUTPUT_STATS}/"
    echo "  Copied ${#files[@]} file(s) to ${NFS_SAMPLE_OUTPUT_STATS}"
    return 0
}

####==================================####
####          RUN                  ####
####==================================####

echo ""
echo "========================================="
echo "  Starting stats pipeline"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
echo "========================================="

if ! run_per_otu_stats; then
    echo "ERROR: per-OTU stats failed"
    tg_send "ERROR: per-OTU stats failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

if ! run_overall_stats; then
    echo "WARNING: overall stats failed (continuing)"
fi

if ! run_join_checkv; then
    echo "ERROR: join with CheckV failed"
    tg_send "ERROR: join with CheckV failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

if ! copy_to_nfs; then
    echo "ERROR: NFS copy failed"
    tg_send "ERROR: NFS copy failed for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
    exit 1
fi

####==================================####
####          COMPLETION              ####
####==================================####

echo ""
echo "========================================="
echo "  Completed job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
echo "  Scratch: ${SAMPLE_OUTPUT_STATS}"
echo "  NFS:     ${NFS_SAMPLE_OUTPUT_STATS}"
echo "========================================="

tg_send "Stats completed: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true
