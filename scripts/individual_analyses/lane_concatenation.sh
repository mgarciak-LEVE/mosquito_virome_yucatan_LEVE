#!/bin/bash

# Script for concatenation of sequencing lanes per sample, with integrity verification and MD5 checksum generation.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 2.1.1 
# 15/04/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw/total_RNA"  
OUTPUT_DIR="${PROJECT_ROOT}/data/raw/total_RNA/cat_files"

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"  # Assuming same directory

# Create output directory
mkdir -p "${OUTPUT_DIR}"

cd "${RAW_DATA_DIR}" || exit 1

echo "Directory: $(pwd)"
tg_send "Directory: $(pwd)"

echo "Output directory: ${OUTPUT_DIR}"
tg_send "Output directory: ${OUTPUT_DIR}"

####==================================####
####       File Verification          ####
####==================================####

echo "Searching for FASTQ files in ${RAW_DATA_DIR}..."
tg_send "Searching for FASTQ files in ${RAW_DATA_DIR}..."

ls -la *.fastq 2>/dev/null || { echo "No FASTQ files found"; exit 1; }

# Check basic integrity of FASTQ files

echo "Checking integrity of FASTQ files..."
tg_send "Checking integrity of FASTQ files..."

for fastq in *.fastq; do
    [[ -f "$fastq" ]] || continue
    lines=$(wc -l < "$fastq")
    if (( lines % 4 != 0 )); then
        echo "$fastq - $lines is not a multiple of 4 → Potentially corrupted file"
        tg_send "$fastq - $lines is not a multiple of 4 → Potentially corrupted file"
    else
        echo "OK: $fastq - $((lines/4)) reads"
        tg_send "OK: $fastq - $((lines/4)) reads"
    fi
done

echo "FASTQ file check completed."
tg_send "FASTQ file check completed."


####======================####
####  Lane Concatenation  ####
####======================####

# Get unique sample names (PMXXXX_SXX)
for sample in $(ls *.fastq | grep -E '_L00[12]_' | sed 's/_L00[12]_.*//' | sort -u); do

    echo "Processing sample ${sample}"
    tg_send "Processing sample ${sample}"

    # Initialize arrays for R1 and R2 files
    r1_files=()
    r2_files=()

    # Search for L001 and L002 files
    for lane in L001 L002; do
        r1_file="${sample}_${lane}_R1_001.fastq"
        r2_file="${sample}_${lane}_R2_001.fastq"

        if [[ -f "$r1_file" ]]; then
            r1_files+=("$r1_file")
            echo "$r1_file found"
            tg_send "$r1_file found"
        fi

        if [[ -f "$r2_file" ]]; then
            r2_files+=("$r2_file")
            echo "$r2_file found"
            tg_send "$r2_file found"
        fi
    done

    # Concatenate R1 files
    if [[ ${#r1_files[@]} -gt 0 ]]; then
        cat "${r1_files[@]}" > "${OUTPUT_DIR}/${sample}_R1.fastq"
        echo "Created ${sample}_R1.fastq from ${#r1_files[@]} lane(s)"
        tg_send "Created ${sample}_R1.fastq from ${#r1_files[@]} lane(s)"
    else
        echo "No R1 files found for sample ${sample}"
        tg_send "No R1 files found for sample ${sample}"
    fi

    # Concatenate R2 files
    if [[ ${#r2_files[@]} -gt 0 ]]; then
        cat "${r2_files[@]}" > "${OUTPUT_DIR}/${sample}_R2.fastq"
        echo "Created ${sample}_R2.fastq from ${#r2_files[@]} lane(s)"
        tg_send "Created ${sample}_R2.fastq from ${#r2_files[@]} lane(s)"
    else
        echo "No R2 files found for sample ${sample}"
        tg_send "No R2 files found for sample ${sample}"
    fi

done

####=============================####
####    Cleanup and Compression  ####
####=============================####

# Optional: Remove original lane files
echo "Do you want to remove the original files? (y/n)"
read -r response
if [[ "$response" == "y" || "$response" == "Y" ]]; then
    rm -f -- *_L001_*.fastq *_L002_*.fastq
    echo "Original lane files removed."
    tg_send "Original files removed from ${RAW_DATA_DIR}"

else
    echo "Original files preserved."
    tg_send "Original files preserved in ${RAW_DATA_DIR}"
fi

# Compress concatenated files
echo "Compressing files."
for file in "${OUTPUT_DIR}"/*.fastq; do
    if [ -f "$file" ]; then
        echo "Compressing: $(basename "$file")"
        tg_send "Compressing: $(basename "$file")"
        gzip "$file"
    fi
done

####=============================####
####      MD5 Generation         ####
####=============================####

# Create directory for MD5 if it doesn't exist
MD5_DIR="${PROJECT_ROOT}/../docs/archivos_md5"
mkdir -p "${MD5_DIR}"

cd "${OUTPUT_DIR}" || exit 1
md5sum *.fastq.gz > "${MD5_DIR}/md5sums_concatenados.txt"
echo ""
echo "MD5 checksums saved to ${MD5_DIR}/md5sums_concatenados.txt"
tg_send "MD5 checksums saved to ${MD5_DIR}/md5sums_concatenados.txt"


echo "Concatenated files are located in: ${OUTPUT_DIR}"
tg_send "Concatenated files are located in: ${OUTPUT_DIR}"