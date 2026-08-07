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
INPUT_MEGAHIT_DIR="${INPUT_DIR_BASE}/MEGAhit"
INPUT_RNASPADES_DIR="${INPUT_DIR_BASE}/rnaSPAdes"
INPUT_METASPADES_DIR="${INPUT_DIR_BASE}/metaSPAdes"
INPUT_METAVIRALSPADES_DIR="${INPUT_DIR_BASE}/metaviralSPAdes"

# Output directory
OUTPUT_DIR="${PROJECT_SCRATCH}/results/identification"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
VIRSORTER_CONTAINER="${CONTAINERS}/spades_3.15.5.sif"
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
echo "  Unmapped Reads Assembly Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Assembly base input directory: ${INPUT_DIR_ASSEMBLY}"
echo "    - metaviralSPAdes: ${INPUT_METAVIRALSPADES_DIR}"
echo "    - rnaSPAdes: ${INPUT_RNASPADES_DIR}"
echo "    - MEGAhit: ${INPUT_MEGAHIT_DIR}"
echo "    - metaSPAdes: ${INPUT_METASPADES_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting unmapped reads assembly" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE-ASSEMBLER LIST  ####
####==================================####

# Create arrays to store sample-assembler combinations
declare -a SAMPLE_ARRAY
declare -a ASSEMBLER_ARRAY
declare -a ASSEMBLY_FILE_ARRAY

echo "Building sample-assembler list..."

# Function to add sample-assembler combination
add_combination() {
    local assembler=$1
    local sample=$2
    local assembly_file=$3
    
    if [[ -f "$assembly_file" ]] && [[ -s "$assembly_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        ASSEMBLER_ARRAY+=("$assembler")
        ASSEMBLY_FILE_ARRAY+=("$assembly_file")
        echo "  Added: $sample ($assembler) - $(basename "$assembly_file")"
    else
        echo "  WARNING: Assembly file not found or empty for $sample ($assembler): $assembly_file"
    fi
}

# Process MEGAhit
if [[ -d "$INPUT_MEGAHIT_DIR" ]]; then
    for sample_dir in "$INPUT_MEGAHIT_DIR"/*/; do
        if [[ -d "$sample_dir" ]]; then
            sample=$(basename "$sample_dir")
            assembly_file="${sample_dir}/final.contigs.fa"
            add_combination "MEGAhit" "$sample" "$assembly_file"
        fi
    done
fi

# Process metaSPAdes
if [[ -d "$INPUT_METASPADES_DIR" ]]; then
    for sample_dir in "$INPUT_METASPADES_DIR"/*/; do
        if [[ -d "$sample_dir" ]]; then
            sample=$(basename "$sample_dir")
            assembly_file="${sample_dir}/contigs.fasta"
            add_combination "metaSPAdes" "$sample" "$assembly_file"
        fi
    done
fi

# Process metaviralSPAdes
if [[ -d "$INPUT_METAVIRALSPADES_DIR" ]]; then
    for sample_dir in "$INPUT_METAVIRALSPADES_DIR"/*/; do
        if [[ -d "$sample_dir" ]]; then
            sample=$(basename "$sample_dir")
            assembly_file="${sample_dir}/contigs.fasta"
            add_combination "metaviralSPAdes" "$sample" "$assembly_file"
        fi
    done
fi

# Process rnaSPAdes
if [[ -d "$INPUT_RNASPADES_DIR" ]]; then
    for sample_dir in "$INPUT_RNASPADES_DIR"/*/; do
        if [[ -d "$sample_dir" ]]; then
            sample=$(basename "$sample_dir")
            if [[ -f "${sample_dir}/hard_filtered_transcripts.fasta" ]]; then
                assembly_file="${sample_dir}/hard_filtered_transcripts.fasta"
            else
                assembly_file="${sample_dir}/transcripts.fasta"
            fi
            add_combination "rnaSPAdes" "$sample" "$assembly_file"
        fi
    done
fi

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

# Get sample and assembler for this job
SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"
ASSEMBLER="${ASSEMBLER_ARRAY[$JOB_INDEX]}"
ASSEMBLY_FILE="${ASSEMBLY_FILE_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Sample: ${SAMPLE}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Assembly file: ${ASSEMBLY_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${SAMPLE} (${ASSEMBLER}) [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true