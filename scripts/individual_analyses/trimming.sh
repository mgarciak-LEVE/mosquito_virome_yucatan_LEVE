#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to validate fastq files.
# 21/07/2026
# Version 2.3.5 (farm-ready)

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
# Working directory on Lustre
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
# Input directory - where FASTQ files are
INPUT_DIR="${PROJECT_SCRATCH}/data/raw/total_RNA/cat_files"
# Output directory for validation results (on Lustre)
OUTPUT_DIR="${PROJECT_SCRATCH}/results/trimmed"

# Container configuration
CONTAINERS="${HOME}/git_repos/${PROJECT_NAME}/containers"
TRIMMOMATIC_CONTAINER="${CONTAINERS}/trimmomatic_0.39.sif"

# Scripts directory (on NFS)
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot (source from NFS)
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Sequence Trimming Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Input: ${INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Trimming Container: ${TRIMMOMATIC_CONTAINER}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting trimming for ${PROJECT_NAME}" 2>/dev/null || true

####==================================####
####            CHECK INPUT           ####
####==================================####

# Check if input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo " ERROR: Input directory ${INPUT_DIR} does not exist"
    tg_send " ERROR: MultiQC input directory not found" 2>/dev/null || true
    exit 1
fi

####==================================####
####          GET FILE LIST          ####
####==================================####

cd "$INPUT_DIR" || exit 1

# Create array of R1 files (R2 will be derived from R1)
r1_files=()
while IFS= read -r line; do
    r1_files+=("$line")
done < <(ls -1 *_R1.fastq *_R1.fastq.gz 2>/dev/null | sort)

# Check if files exist
if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "ERROR: No R1 FASTQ files found in ${INPUT_DIR}"
    tg_send "ERROR: No R1 FASTQ files found" 2>/dev/null || true
    exit 1
fi

echo "Found ${#r1_files[@]} R1 files"

# Get the R1 file for this array task
FILE_INDEX=$((LSB_JOBINDEX - 1))

if [[ $FILE_INDEX -ge ${#r1_files[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#r1_files[@]})"
    exit 1
fi

R1_FILE="${r1_files[$FILE_INDEX]}"
R1_NAME=$(basename "$R1_FILE")

# Derive R2 file name
R2_FILE="${R1_FILE/_R1./_R2.}"
R2_NAME=$(basename "$R2_FILE")

# Check if R2 exists
if [[ ! -f "$R2_FILE" ]]; then
    echo "ERROR: R2 file ${R2_FILE} not found for ${R1_NAME}"
    exit 1
fi

echo "========================================="
echo "  Trimmomatic Array Task ${LSB_JOBINDEX}/${#r1_files[@]}"
echo "  R1: ${R1_NAME}"
echo "  R2: ${R2_NAME}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Trimmomatic: ${R1_NAME} and ${R2_NAME} (${LSB_JOBINDEX}/${#r1_files[@]})" 2>/dev/null || true

####=========================####
####         TRIMMING        ####
####=========================####

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check if container exists
if [[ ! -f "$TRIMMOMATIC_CONTAINER" ]]; then
    echo "ERROR: Container ${TRIMMOMATIC_CONTAINER} not found"
    tg_send "ERROR: Trimmomatic container not found" 2>/dev/null || true
    exit 1
fi

# Extract sample name (PMXXXX)
sample_name=$(basename "$R1_FILE" | sed 's/_.*//')
echo "Sample: ${sample_name}"

# Output file names
trimmed_dir="${OUTPUT_DIR}/${sample_name}"
mkdir -p "$trimmed_dir"

# Output files
R1_PAIRED="${trimmed_dir}/${sample_name}_R1_paired.fastq"
R1_UNPAIRED="${trimmed_dir}/${sample_name}_R1_unpaired.fastq"
R2_PAIRED="${trimmed_dir}/${sample_name}_R2_paired.fastq"
R2_UNPAIRED="${trimmed_dir}/${sample_name}_R2_unpaired.fastq"

echo "Running Trimmomatic..."
tg_send "Running Trimmomatic for ${sample_name}" 2>/dev/null || true

if apptainer exec \
    --bind "${INPUT_DIR}:/input:ro" \
    --bind "${trimmed_dir}:/output" \
    "$TRIMMOMATIC_CONTAINER" \
    trimmomatic PE \
    -threads 4 \
    -phred33 \
    "/input/${R1_NAME}" \
    "/input/${R2_NAME}" \
    "/output/${sample_name}_R1_paired.fastq" \
    "/output/${sample_name}_R1_unpaired.fastq" \
    "/output/${sample_name}_R2_paired.fastq" \
    "/output/${sample_name}_R2_unpaired.fastq" \
    ILLUMINACLIP:/usr/local/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:20 \
    TRAILING:20 \
    SLIDINGWINDOW:3:25 \
    MINLEN:50 2>&1; then
    
    echo "Trimmomatic completed for ${sample_name}"
    tg_send "rimmomatic: ${sample_name} complete" 2>/dev/null || true
    
    # Count reads in output files
    if [[ -f "$R1_PAIRED" ]]; then
        reads=$(wc -l < "$R1_PAIRED")
        echo "  Paired reads: $((reads/4))"
    fi
    
else
    echo "Trimmomatic FAILED for ${sample_name}"
    tg_send "Trimmomatic: ${sample_name} FAILED" 2>/dev/null || true
    exit 1
fi

echo "Done processing: ${sample_name}"
