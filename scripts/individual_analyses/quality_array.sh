#!/bin/bash
# Single FastQC task for job array
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 17/07/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

# Directory where scripts are
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/team222/projects"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
INPUT_FASTQC="${PROJECT_SCRATCH}/data/raw/total_RNA/cat_files"
OUTPUT_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_FASTQC="${OUTPUT_BASE}/untrimmed_qc/fastqc"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
FASTQC_CONTAINER="${CONTAINERS}/fastqc_staphb.sif"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          GET FILE LIST          ####
####==================================####

cd "$INPUT_FASTQC" || exit 1

# Create array of files
files=()
while IFS= read -r line; do
    files+=("$line")
done < <(ls -1 *.fastq *.fastq.gz 2>/dev/null | sort)

# Check if files exist
if [[ ${#files[@]} -eq 0 ]]; then
    echo "No FASTQ files found in ${INPUT_FASTQC}"
    exit 1
fi

# Get the file for this array task
FILE_INDEX=$((LSB_JOBINDEX - 1))

if [[ $FILE_INDEX -ge ${#files[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#files[@]})"
    exit 1
fi

file="${files[$FILE_INDEX]}"
file_name=$(basename "$file")

echo "========================================="
echo "  FastQC Array Task ${LSB_JOBINDEX}/${#files[@]}"
echo "  File: ${file_name}"
echo "  Date: $(date)"
echo "========================================="

tg_send "FastQC: ${file_name} (${LSB_JOBINDEX}/${#files[@]})" 2>/dev/null || true

####==================================####
####         RUN FASTQC              ####
####==================================####

# Create output directory
mkdir -p "${OUTPUT_FASTQC}"

# Check if container exists
if [[ ! -f "$FASTQC_CONTAINER" ]]; then
    echo "ERROR: Container ${FASTQC_CONTAINER} not found"
    tg_send "ERROR: FastQC container not found" 2>/dev/null || true
    exit 1
fi

# Run FastQC with container
echo "Processing: ${file_name}"
tg_send "Processing: ${file_name}" 2>/dev/null || true

if apptainer exec \
    --bind "${INPUT_FASTQC}:/input:ro" \
    --bind "${OUTPUT_FASTQC}:/output" \
    --bind /usr/share/fonts:/usr/share/fonts:ro \
    --bind /usr/share/fontconfig:/usr/share/fontconfig:ro \
    --bind /etc/fonts:/etc/fonts:ro \
    "$FASTQC_CONTAINER" \
    fastqc "/input/${file_name}" -o "/output" 2>&1; then
    echo "FastQC completed for ${file_name}"
    tg_send "FastQC: ${file_name} complete" 2>/dev/null || true
else
    echo "FastQC FAILED for ${file_name}"
    tg_send "FastQC: ${file_name} FAILED" 2>/dev/null || true
    exit 1
fi

echo "Done processing: ${file_name}"
