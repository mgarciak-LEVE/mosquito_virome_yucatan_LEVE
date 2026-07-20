#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script for concatenation of sequencing lanes per sample, with integrity verification and MD5 checksum generation.
# 08/07/2026
# Version 2.1.5 (farm-ready - handles .gz files)

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
# Permanent storage
PERMANENT_BASE="/nfs/users/nfs_j/jr46"

# Scratch storage
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
# Working directory on Lustre
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directory - where FASTQ files are
INPUT_DIR="${PROJECT_SCRATCH}/data/raw/total_RNA"

# Output directory for validation results (on Lustre)
OUTPUT_DIR="${INPUT_DIR}/cat_files"

# Permanent results directory
PERMANENT_RESULTS="${PERMANENT_BASE}/git_repos/${PROJECT_NAME}/results"

# Scripts directory (on NFS)
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# --- MD5 DIRECTORIES ---
# MD5 on Lustre
MD5_LUSTRE="${PROJECT_SCRATCH}/docs/md5_files"

# MD5 on NFS 
MD5_NFS="${PERMANENT_RESULTS}/../docs/md5_files"

# Telegram bot (source from NFS)
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####           PRINT CONFIG           ####
####==================================####

echo "========================================="
echo "  Lane Concatenation"
echo "  Project: ${PROJECT_NAME}"
echo "  Input directory: ${INPUT_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  MD5 (Lustre): ${MD5_LUSTRE}"
echo "  MD5 (NFS): ${MD5_NFS}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting lane concatenation for ${PROJECT_NAME}" 2>/dev/null || true

####===========================####
####     FILE VERIFICATION     ####
####===========================####

cd "${INPUT_DIR}" || {
    echo "ERROR: Cannot cd to ${INPUT_DIR}"
    tg_send "ERROR: Cannot cd to ${INPUT_DIR}" 2>/dev/null || true
    exit 1
}

echo "Searching for FASTQ files in ${INPUT_DIR}..."
tg_send "Searching for FASTQ files in ${INPUT_DIR}..." 2>/dev/null || true

if ! ls *.fastq.gz *.fastq 2>/dev/null | grep -q .; then
    echo "No FASTQ files found in ${INPUT_DIR}"
    tg_send "No FASTQ files found in ${INPUT_DIR}" 2>/dev/null || true
    exit 1
fi

# List the files (show .fastq.gz if they exist, otherwise .fastq)
ls -la *.fastq.gz 2>/dev/null || ls -la *.fastq 2>/dev/null

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Checking integrity of FASTQ files..."
tg_send "Checking integrity of FASTQ files..." 2>/dev/null || true

for fastq in *.fastq *.fastq.gz; do
    [[ -f "$fastq" ]] || continue
    
    # If file is compressed, check it without decompressing
    if [[ "$fastq" == *.gz ]]; then
        # Check if it's a valid gzip file
        if gunzip -t "$fastq" 2>/dev/null; then
            # Get line count from compressed file
            lines=$(zcat "$fastq" | wc -l)
            echo "OK (compressed): $fastq - $((lines/4)) reads"
        else
            echo "WARNING: $fastq is corrupted (invalid gzip)"
            tg_send "WARNING: $fastq is corrupted (invalid gzip)" 2>/dev/null || true
            continue
        fi
    else
        lines=$(wc -l < "$fastq")
        if (( lines % 4 != 0 )); then
            echo "WARNING: $fastq - $lines is not a multiple of 4 → Potentially corrupted file"
            tg_send "WARNING: $fastq - $lines is not a multiple of 4" 2>/dev/null || true
        else
            echo "OK: $fastq - $((lines/4)) reads"
        fi
    fi
done

echo "FASTQ file check completed."
tg_send "FASTQ file check completed." 2>/dev/null || true

####======================####
####  LANE CONCATENATION  ####
####======================####

# Get unique sample names (PMXXXX_SXX)
for sample in $(ls *.fastq *.fastq.gz 2>/dev/null | grep -E '_L00[12]_' | sed 's/_L00[12]_.*//' | sort -u); do

    echo "Processing sample ${sample}"
    tg_send "Processing sample ${sample}" 2>/dev/null || true

    # Initialize arrays for R1 and R2 files
    r1_files=()
    r2_files=()

    # Search for L001 and L002 files (both .fastq and .fastq.gz)
    for lane in L001 L002; do
        for ext in fastq fastq.gz; do
            r1_file="${sample}_${lane}_R1_001.${ext}"
            r2_file="${sample}_${lane}_R2_001.${ext}"

            if [[ -f "$r1_file" ]]; then
                r1_files+=("$r1_file")
                echo "  $r1_file found"
            fi

            if [[ -f "$r2_file" ]]; then
                r2_files+=("$r2_file")
                echo "  $r2_file found"
            fi
        done
    done

    # Concatenate R1 files
    if [[ ${#r1_files[@]} -gt 0 ]]; then
        # Determine if files are compressed
        if [[ "${r1_files[0]}" == *.gz ]]; then
            # Decompress, concatenate, then recompress
            echo "  Decompressing ${#r1_files[@]} R1 files..."
            temp_r1="${OUTPUT_DIR}/${sample}_R1_temp.fastq"
            for f in "${r1_files[@]}"; do
                zcat "$f" >> "$temp_r1"
            done
            mv "$temp_r1" "${OUTPUT_DIR}/${sample}_R1.fastq"
            echo "  Created ${sample}_R1.fastq from ${#r1_files[@]} lane(s)"
        else
            cat "${r1_files[@]}" > "${OUTPUT_DIR}/${sample}_R1.fastq"
            echo "  Created ${sample}_R1.fastq from ${#r1_files[@]} lane(s)"
        fi
        tg_send "  Created ${sample}_R1.fastq from ${#r1_files[@]} lane(s)" 2>/dev/null || true
    else
        echo "  WARNING: No R1 files found for sample ${sample}"
    fi

    # Concatenate R2 files
    if [[ ${#r2_files[@]} -gt 0 ]]; then
        if [[ "${r2_files[0]}" == *.gz ]]; then
            echo "  Decompressing ${#r2_files[@]} R2 files..."
            temp_r2="${OUTPUT_DIR}/${sample}_R2_temp.fastq"
            for f in "${r2_files[@]}"; do
                zcat "$f" >> "$temp_r2"
            done
            mv "$temp_r2" "${OUTPUT_DIR}/${sample}_R2.fastq"
            echo "  Created ${sample}_R2.fastq from ${#r2_files[@]} lane(s)"
        else
            cat "${r2_files[@]}" > "${OUTPUT_DIR}/${sample}_R2.fastq"
            echo "  Created ${sample}_R2.fastq from ${#r2_files[@]} lane(s)"
        fi
        tg_send "  Created ${sample}_R2.fastq from ${#r2_files[@]} lane(s)" 2>/dev/null || true
    else
        echo "  WARNING: No R2 files found for sample ${sample}"
    fi

done

echo "Lane concatenation completed."
tg_send "Lane concatenation completed." 2>/dev/null || true

####==================####
####    COMPRESSION   ####
####==================####

# Move to output directory
cd "${OUTPUT_DIR}" || exit 1

# Compress concatenated files
echo "Compressing concatenated files..."
tg_send "Compressing concatenated files..." 2>/dev/null || true

for file in *.fastq; do
    if [[ -f "$file" ]]; then
        echo "  Compressing: $(basename "$file")"
        gzip "$file"
    fi
done

echo "Compression completed."
tg_send "Compression completed." 2>/dev/null || true

####==================================####
####      REMOVE ORIGINAL FILES      ####
####==================================####

# Remove original lane files (OPTIONAL)
KEEP_ORIGINAL="no"  # Change to "yes" if you want to keep originals

if [[ "$KEEP_ORIGINAL" == "no" ]]; then
    echo "Removing original lane files..."
    cd "${INPUT_DIR}" || exit 1
    rm -f -- *_L001_*.fastq *_L001_*.fastq.gz *_L002_*.fastq *_L002_*.fastq.gz
    echo "Original lane files removed."
    tg_send "Original files removed from ${INPUT_DIR}" 2>/dev/null || true
else
    echo "Original files preserved."
    tg_send "Original files preserved in ${INPUT_DIR}" 2>/dev/null || true
fi

####=============================####
####      MD5 GENERATION         ####
####=============================####

echo "Generating MD5 checksums..."
tg_send "Generating MD5 checksums..." 2>/dev/null || true

# Create MD5 directories
mkdir -p "${MD5_LUSTRE}"
mkdir -p "${MD5_NFS}"

# Generate MD5 checksums
cd "${OUTPUT_DIR}" || exit 1

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
MD5_FILENAME="concatenated_md5sums_${TIMESTAMP}.txt"

# Generate MD5 on Lustre
echo "  Writing MD5 to Lustre: ${MD5_LUSTRE}/${MD5_FILENAME}"
md5sum *.fastq.gz > "${MD5_LUSTRE}/${MD5_FILENAME}"

# Copy to NFS (permanent storage)
echo "  Copying MD5 to NFS: ${MD5_NFS}/${MD5_FILENAME}"
cp "${MD5_LUSTRE}/${MD5_FILENAME}" "${MD5_NFS}/"

# Verify MD5 files exist
if [[ -f "${MD5_LUSTRE}/${MD5_FILENAME}" ]] && [[ -f "${MD5_NFS}/${MD5_FILENAME}" ]]; then
    echo "MD5 checksums successfully written to both locations"
    tg_send "MD5 checksums saved to Lustre and NFS" 2>/dev/null || true
else
    echo "ERROR: Failed to write MD5 checksums"
    tg_send "ERROR: MD5 checksum generation failed" 2>/dev/null || true
fi

####=============================####
####            SUMMARY          ####
####=============================####

echo "========================================="
echo "  Lane Concatenation Summary"
echo "  Date: $(date)"
echo "  Input directory: ${INPUT_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  MD5 (Lustre): ${MD5_LUSTRE}/${MD5_FILENAME}"
echo "  MD5 (NFS): ${MD5_NFS}/${MD5_FILENAME}"
echo "========================================="

tg_send "Lane concatenation complete: ${OUTPUT_DIR}" 2>/dev/null || true
