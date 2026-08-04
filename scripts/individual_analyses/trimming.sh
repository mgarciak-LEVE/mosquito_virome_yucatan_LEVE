#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to validate fastq files.
# 05/08/2026
# Version 2.4.0 (farm-ready)

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
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
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

####==================================####
####        DECOMPRESS INPUT         ####
####==================================####

# Create a temporary working directory for this job
TEMP_DIR="${INPUT_DIR}/temp_${LSB_JOBINDEX}_${$}"
mkdir -p "$TEMP_DIR"

# Cleanup function for temp directory
cleanup_temp() {
    if [[ -d "$TEMP_DIR" ]]; then
        echo "Cleaning up temporary directory: ${TEMP_DIR}"
        rm -rf "$TEMP_DIR"
    fi
}

# Set trap to cleanup on script exit (success or failure)
trap cleanup_temp EXIT

echo "Creating temporary working directory: ${TEMP_DIR}"

# Copy and decompress files to temp directory
if [[ "$R1_NAME" == *.gz ]]; then
    echo "Decompressing: ${R1_NAME} to temp..."
    zcat "$R1_FILE" > "${TEMP_DIR}/${R1_NAME%.gz}"
    R1_NAME="${R1_NAME%.gz}"
    R1_FILE="${TEMP_DIR}/${R1_NAME}"
    echo "Decompressed to: ${R1_NAME}"
else
    cp "$R1_FILE" "${TEMP_DIR}/$R1_NAME"
    R1_FILE="${TEMP_DIR}/${R1_NAME}"
fi

if [[ "$R2_NAME" == *.gz ]]; then
    echo "Decompressing: ${R2_NAME} to temp..."
    zcat "$R2_FILE" > "${TEMP_DIR}/${R2_NAME%.gz}"
    R2_NAME="${R2_NAME%.gz}"
    R2_FILE="${TEMP_DIR}/${R2_NAME}"
    echo "Decompressed to: ${R2_NAME}"
else
    cp "$R2_FILE" "${TEMP_DIR}/$R2_NAME"
    R2_FILE="${TEMP_DIR}/${R2_NAME}"
fi

echo "Input files ready for trimming:"
echo "  R1: ${R1_NAME}"
echo "  R2: ${R2_NAME}"

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

# Output directory for this sample
trimmed_dir="${OUTPUT_DIR}/${sample_name}"
mkdir -p "$trimmed_dir"

# Output files
R1_PAIRED="${trimmed_dir}/${sample_name}_R1_paired.fastq"
R1_UNPAIRED="${trimmed_dir}/${sample_name}_R1_unpaired.fastq"
R2_PAIRED="${trimmed_dir}/${sample_name}_R2_paired.fastq"
R2_UNPAIRED="${trimmed_dir}/${sample_name}_R2_unpaired.fastq"

echo "Running Trimmomatic..."
tg_send "Running Trimmomatic for ${sample_name}" 2>/dev/null || true

# FIX 1: Bind mount TEMP_DIR (not INPUT_DIR)
if apptainer exec \
    --bind "${TEMP_DIR}:/input:ro" \
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
    ILLUMINACLIP:/usr/local/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:8:2:keepBothReads \
    LEADING:20 \
    TRAILING:20 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36 2>&1; then
    
    # ILLUMINACLIP 2:30:10:2 simple clip threshold to catch adapters on short inserts [<seedMismatches>:<palindromeClipThreshold>:<simpleClipThreshold>:<minAdapterLength>]
    # TRAILING 3 since quality drops below 3 in the last 97 nucleotides
    # SLIDINGWINDOW 4:15 to remove reads below 15 quality under a 4 nt window
    # MINLEN 36 to capture some small RNAs and remove degradation products

    echo "Trimmomatic completed for ${sample_name}"
    tg_send "Trimmomatic: ${sample_name} complete" 2>/dev/null || true
    
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

# FIX 2: Remove this - no need to recompress since we used temp files
# REMOVE the entire "RECOMPRESS INPUT" section

####==================================####
####        COMPRESS OUTPUT          ####
####==================================####

echo "Compressing output files..."

# Compress paired-end output files
if [[ -f "$R1_PAIRED" ]]; then
    echo "  Compressing: ${sample_name}_R1_paired.fastq"
    gzip "$R1_PAIRED"
fi

if [[ -f "$R2_PAIRED" ]]; then
    echo "  Compressing: ${sample_name}_R2_paired.fastq"
    gzip "$R2_PAIRED"
fi

# Compress unpaired output files if they exist and are not empty
if [[ -f "$R1_UNPAIRED" ]]; then
    if [[ -s "$R1_UNPAIRED" ]]; then
        echo "  Compressing: ${sample_name}_R1_unpaired.fastq"
        gzip "$R1_UNPAIRED"
    else
        echo "  Removing empty file: ${sample_name}_R1_unpaired.fastq"
        rm "$R1_UNPAIRED"
    fi
fi

if [[ -f "$R2_UNPAIRED" ]]; then
    if [[ -s "$R2_UNPAIRED" ]]; then
        echo "  Compressing: ${sample_name}_R2_unpaired.fastq"
        gzip "$R2_UNPAIRED"
    else
        echo "  Removing empty file: ${sample_name}_R2_unpaired.fastq"
        rm "$R2_UNPAIRED"
    fi
fi

# FIX 3: Remove the "RECOMPRESS INPUT" section entirely
# TEMP_DIR will be cleaned up by the trap on EXIT

echo ""
echo "========================================="
echo "  Processing complete for: ${sample_name}"
echo "  Output files in: ${trimmed_dir}"
echo "  End time: $(date)"
echo "========================================="

# List final output files
echo ""
echo "Final output files:"
ls -la "${trimmed_dir}"/* 2>/dev/null || echo "No output files found"

tg_send "Trimmomatic completed: ${sample_name}" 2>/dev/null || true
