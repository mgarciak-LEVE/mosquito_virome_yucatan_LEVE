#!/bin/bash

# Script for lane concatenation
# Author: Jorge Alberto Castro Rodríguez
# Ver. 2.1.0
# 15/04/2026


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw/total_RNA"

cd ${RAW_DATA_DIR}
output_dir="${RAW_DATA_DIR}/cat_files"
mkdir -p "${output_dir}"

# Gunzip files and keep originals
for file in *.gz; do
    echo "Extracting: $file"
    gunzip -k "$file"
done

# Checking if file corruption ocurred.

echo "Checking for corrupt FASTQ files..."
for fastq in *.fastq; do
    lines=$(wc -l < "$fastq")
    if (( lines % 4 != 0 )); then
        echo "❌ CORRUPT: $fastq - $lines lines (not multiple of 4)"
    fi
done

echo "Corrupt file checking DONE"

# Get unique sample names (everything before _L00[12]_)
for sample in $(ls *.fastq | grep -E '_L00[12]_' | cut -d'_' -f1-2 | sort -u); do
    
    # Initialize arrays to collect files
    r1_files=()
    r2_files=()
    
    # Check for L001 and L002 files
    for lane in L001 L002; do
        r1_file="${sample}_${lane}_R1_001.fastq"
        r2_file="${sample}_${lane}_R2_001.fastq"
        
        if [[ -f "$r1_file" ]]; then
            r1_files+=("$r1_file")
        fi
        
        if [[ -f "$r2_file" ]]; then
            r2_files+=("$r2_file")
        fi
    done
    
    # Concatenate R1 files if any exist
    if [[ ${#r1_files[@]} -gt 0 ]]; then
        cat "${r1_files[@]}" > "${output_dir}/$(basename ${sample})_R1.fastq"
        echo "Created: ${output_dir}/$(basename ${sample})_R1.fastq from ${#r1_files[@]} lane(s)"
    else
        echo "Warning: No R1 files found for sample ${sample}"
    fi
    
    # Concatenate R2 files if any exist
    if [[ ${#r2_files[@]} -gt 0 ]]; then
        cat "${r2_files[@]}" > "${output_dir}/$(basename ${sample})_R2.fastq"
        echo "Created: ${output_dir}/$(basename ${sample})_R2.fastq from ${#r2_files[@]} lane(s)"
    else
        echo "Warning: No R2 files found for sample ${sample}"
    fi
done

# Optional: Remove original lane files

rm -f -- *_L001_*.fastq *_L002_*.fastq

# Compressing of .fastq files
for file in "${output_dir}"/*.fastq; do
    if [ -f "$file" ]; then
        echo "Compressing: $file" 	
        gzip "$file"
    fi
done

echo "Process completed..."