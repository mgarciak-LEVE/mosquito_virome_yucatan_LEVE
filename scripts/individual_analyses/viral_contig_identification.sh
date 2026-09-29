#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to identify assembled contigs
# 24/09/2026
# Ver. 1.4.0

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

# Database paths
VIRSORTER2_DB="/lustre/scratch126/tol/teams/lawniczak/users/jr46/databases/virsorter2_db/db"
VITAP_DB="/lustre/scratch126/tol/teams/lawniczak/users/jr46/databases/vitap_db/db/DB_hybrid_MSL37_RefSeq209_IMGVR"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
DEEPVIRFINDER_CONTAINER="${CONTAINERS}/deepvirfinder.sif"
VITAP_CONTAINER="${CONTAINERS}/vitap_1.12--pyhdfd78af_0.sif"

# Conda environment path
CONDA_ENV_VIRSORTER2="/software/treeoflife/conda/users/envs/team222/jr46/vs2"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

# Explicitly initialize conda for non-interactive batch jobs
source /software/treeoflife/conda/miniconda/etc/profile.d/conda.sh
conda activate /software/treeoflife/conda/users/envs/team222/jr46/vs2

# Apptainer load
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
echo "  Conda env: ${CONDA_ENV_VIRSORTER2}"
echo "  Input directories (mapping_tool/assembler/read_type/sample):"
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
declare -a READ_TYPE_ARRAY
declare -a ASSEMBLER_ARRAY
declare -a ASSEMBLY_FILE_ARRAY

echo "Building sample-assembler list..."

# Function to add sample-assembler combination
add_combination() {
    local mapping_tool=$1
    local read_type=$2         
    local assembler=$3
    local sample=$4
    local assembly_file=$5

    if [[ -f "$assembly_file" ]] && [[ -s "$assembly_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        READ_TYPE_ARRAY+=("$read_type")
        ASSEMBLER_ARRAY+=("$assembler")
        ASSEMBLY_FILE_ARRAY+=("$assembly_file")
        echo "  Added: $mapping_tool/$read_type/$assembler/$sample - $(basename "$assembly_file")"
    else
        echo "  WARNING: Assembly file not found or empty for $mapping_tool/$read_type/$assembler/$sample: $assembly_file"
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
                        echo "        WARNING: Unknown assembler: $assembler"
                        continue
                        ;;
                esac

                add_combination "$mapping_tool" "$read_type" "$assembler" "$sample" "$assembly_file"
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
    echo "ERROR: Invlid job index ${LSB_JOBINDEX} (max ${#SAMPLE_ARRAY[@]})"
    exit 1
fi

# Get sample, mapping tool, and assembler for this job
SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"
MAPPING_TOOL="${MAPPING_ARRAY[$JOB_INDEX]}"
READ_TYPE="${READ_TYPE_ARRAY[$JOB_INDEX]}"       
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
SAMPLE_OUTPUT="${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
mkdir -p "${SAMPLE_OUTPUT}"

# Create writable directories for VirSorter2
VS2_CACHE="${SAMPLE_OUTPUT}/virsorter2_cache"
VS2_TMP="${SAMPLE_OUTPUT}/virsorter2_tmp"
CONDA_PKGS="${SAMPLE_OUTPUT}/conda_pkgs"
mkdir -p "$VS2_CACHE" "$VS2_TMP" "$CONDA_PKGS"

# Change to output directory
cd "$SAMPLE_OUTPUT" || exit 1

run_virsorter2() {
    echo ""
    echo "--- Running VirSorter2 (native conda) ---"

    local output_dir="${SAMPLE_OUTPUT}/virsorter2"
    mkdir -p "$output_dir"

    if [[ -s "${output_dir}/final-viral-combined.fa" ]]; then
        echo "  VirSorter2 already completed for ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi

    if [[ ! -f "${VIRSORTER2_DB}/hmm/viral/combined.hmm" ]]; then
        echo "  WARNING: VirSorter2 DB missing at ${VIRSORTER2_DB}"
        return 0
    fi

    if ! command -v virsorter >/dev/null 2>&1; then
        echo "  ERROR: virsorter not on PATH after conda activate"
        conda deactivate
        return 1
    fi

    echo "  virsorter at: $(command -v virsorter)"
    echo "  interpreter:  $(head -1 "$(command -v virsorter)")"

    virsorter run \
        -i "${ASSEMBLY_FILE}" \
        -w "${output_dir}" \
        --db-dir "${VIRSORTER2_DB}" \
        --include-groups "dsDNAphage,ssDNA,RNA" \
        --min-length 1000 \
        --min-score 0.5 \
        --use-conda-off \
        -j "${CORES}" \
        all \
        2>&1 | tee "${output_dir}/virsorter2.log"

    local vs2_exit=${PIPESTATUS[0]}

    if [[ ${vs2_exit} -ne 0 ]]; then
        echo "  VirSorter2 failed (exit ${vs2_exit})"
        return 1
    fi

    if [[ -s "${output_dir}/final-viral-combined.fa" ]]; then
        local n
        n=$(grep -c '^>' "${output_dir}/final-viral-combined.fa")
        echo "  VirSorter2 completed (${n} viral contigs)"
        return 0
    else
        echo "  VirSorter2 finished but no output — check ${output_dir}/virsorter2.log"
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
    local existing_pred
    existing_pred=$(find "${output_dir}" -maxdepth 1 -name "*_dvfpred.txt" 2>/dev/null | head -n1)
    if [[ -n "${existing_pred}" ]] && [[ -s "${existing_pred}" ]]; then
        echo "  DeepVirFinder already completed for ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running DeepVirFinder..."
    
    # Ensure CORES is set
    if [[ -z "${CORES}" ]]; then
        CORES=4
        echo "  WARNING: CORES not set, defaulting to ${CORES}"
    fi
    
    # Create writable directory for DeepVirFinder cache
    local dvf_cache="${SAMPLE_OUTPUT}/dvf_cache"
    mkdir -p "$dvf_cache"
    
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        --bind "${VS2_TMP}":/tmp \
        --bind "${dvf_cache}":/root/.dvf \
        --env "TMPDIR=/tmp" \
        "${DEEPVIRFINDER_CONTAINER}" \
        python /DeepVirFinder/dvf.py \
            -i "/input/$(basename "$ASSEMBLY_FILE")" \
            -o "/output/deepvirfinder" \
            -l 1000 \
            -c "${CORES}" \
            2>&1 | tee "${output_dir}/deepvirfinder.log"

    local deepvir_exit=${PIPESTATUS[0]}

    if [[ ${deepvir_exit} -ne 0 ]]; then
        echo "  DeepVirFinder failed (exit ${deepvir_exit})"
        return 1
    fi

    local dvf_pred
    dvf_pred=$(find "${output_dir}" -maxdepth 1 -name "*_dvfpred.txt" 2>/dev/null | head -n1)

    if [[ -n "${dvf_pred}" ]] && [[ -s "${dvf_pred}" ]]; then
        awk '$4 <= 0.05' "${dvf_pred}" > "${output_dir}/viral_contigs.txt"
        local viral_count
        viral_count=$(wc -l < "${output_dir}/viral_contigs.txt")
        echo "  Found ${viral_count} viral contigs (p < 0.05)"
        return 0
    else
        echo "  DeepVirFinder finished but produced no prediction file"
        echo "  Looked for: ${output_dir}/*_dvfpred.txt"
        return 1
    fi
}

####==================================####
####       FUNCTION: RUN VITAP        ####
####==================================####

run_vitap() {
    echo ""
    echo "--- Running VITAP ---"

    local output_dir="${SAMPLE_OUTPUT}/vitap"
    mkdir -p "$output_dir"

    if [[ -f "${output_dir}/best_determined_lineages.tsv" ]] && \
    [[ -s "${output_dir}/best_determined_lineages.tsv" ]]; then
    echo "  VITAP already completed for ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
    return 0
    fi

    # ---- Input: original assembly ----
    local input_contigs="${ASSEMBLY_FILE}"
    if [[ ! -s "${input_contigs}" ]]; then
        echo "  ERROR: Assembly file missing or empty: ${input_contigs}"
        return 1
    fi

    local contig_count
    contig_count=$(grep -c '^>' "${input_contigs}")
    echo "  Input: $(basename "${input_contigs}") (${contig_count} contigs)"

    # ---- Sanity checks on config ----
    if [[ -z "${VITAP_DB:-}" ]] || [[ ! -d "${VITAP_DB}" ]]; then
        echo "  WARNING: VITAP_DB not set or not a directory: '${VITAP_DB:-<unset>}'"
        echo "  Skipping VITAP"
        return 0
    fi

    if [[ ! -f "${VITAP_CONTAINER}" ]]; then
        echo "  WARNING: VITAP container not found: ${VITAP_CONTAINER}"
        echo "  Skipping VITAP"
        return 0
    fi

    echo "  Running VITAP..."

    # Per-sample tmp (VITAP writes a lot of intermediates)
    local vitap_tmp="${SAMPLE_OUTPUT}/vitap_tmp"
    mkdir -p "${vitap_tmp}"

    apptainer exec \
        --bind "$(dirname "${input_contigs}")":/input:ro \
        --bind "${output_dir}":/output \
        --bind "${VITAP_DB}":/vitap_db:ro \
        --bind "${vitap_tmp}":/tmp \
        --env "TMPDIR=/tmp" \
        "${VITAP_CONTAINER}" \
        VITAP assignment \
            -i "/input/$(basename "${input_contigs}")" \
            -o "/output" \
            -d "/vitap_db" \
        2>&1 | tee "${output_dir}/vitap.log"

    local vitap_exit=${PIPESTATUS[0]}

    if [[ ${vitap_exit} -ne 0 ]]; then
        echo "  VITAP failed (exit ${vitap_exit})"
        return 1
    fi

    local result_tsv="${output_dir}/best_determined_lineages.tsv"
    if [[ -f "${result_tsv}" ]] && [[ -s "${result_tsv}" ]]; then
        local assigned
        assigned=$(tail -n +2 "${result_tsv}" | wc -l)
        echo "  VITAP completed successfully (${assigned} contigs with assignments)"
        return 0
    else
        echo "  VITAP finished but no best_determined_lineages.tsv — check ${output_dir}/vitap.log"
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
Read Type: ${READ_TYPE}
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
else
    echo "VirSorter2 results not available"
fi)

--- DeepVirFinder Results ---
$(if [[ -f "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]]; then
    viral_count=$(wc -l < "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" 2>/dev/null || echo "0")
    echo "Viral contigs found (p < 0.05): $viral_count"
else
    echo "DeepVirFinder results not available"
fi)

--- VITAP Results ---
$(if [[ -f "${SAMPLE_OUTPUT}/vitap/best_determined_lineages.tsv" ]]; then
    viral_count=$(awk -F'\t' 'NR>1 && ($2 ~ /Riboviria|Duplodnaviria|Monodnaviria/) {print $1}' \
        "${SAMPLE_OUTPUT}/vitap/best_determined_lineages.tsv" | wc -l)
    echo "Viral contigs detected (realm filter): $viral_count"
else
    echo "VITAP results not available"
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
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
echo "========================================="

# Run each tool
run_virsorter2
run_deepvirfinder
run_vitap

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

tg_send "Completed: ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true
