#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to validate fastq files.
# 13/07/2026
# Version 2.3.5 (farm-ready)

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
# Permanent storage
PERMANENT_BASE="/nfs/team222/projects"

# Scratch storage
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
# Working directory on Lustre
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directory - where FASTQ files are
INPUT_DIR="${PROJECT_SCRATCH}/data/raw/total_RNA/cat_files"

# Output directory for validation results (on Lustre)
OUTPUT_DIR="${PROJECT_SCRATCH}/results/untrimmed_qc/stats"

# Permanent results directory
PERMANENT_RESULTS="${PERMANENT_BASE}/${PROJECT_NAME}/docs/fastq_validation"

# Scripts directory (on NFS)
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot (source from NFS)
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

# Print configuration
echo "========================================="
echo "  FASTQ Validation (OPTIMIZED)"
echo "  Project: ${PROJECT_NAME}"
echo "  Input directory: ${INPUT_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting FASTQ validation (optimized) for ${PROJECT_NAME}" 2>/dev/null || true

####==================================####
####     FASTQ INTEGRITY FUNCTION     ####
####==================================####

fastq_val() {
    local file="$1"
    local file_name=$(basename "$file")
    
    # Stats file
    stats_file="${OUTPUT_DIR}/stats_validation.csv"
    
    # Create directory if it doesn't exist
    mkdir -p "$(dirname "$stats_file")"
    
    # Create CSV with headers only once
    if [[ ! -f "$stats_file" ]]; then
        echo "file,reads,lines,errors,gc_percentage,total_n,size_mb" > "$stats_file"
    fi
    
    # File size
    local bytes=$(stat -c%s "$file" 2>/dev/null || echo "0")
    local size_mb=$(echo "scale=2; $bytes / 1048576" | bc 2>/dev/null || echo "0")
    

    awk_script=$(cat <<'EOF'
BEGIN {
    lines=0; reads=0; errors=0; gc_total=0; bases_total=0; total_n=0; seq_len=0
}
{
    lines++
    pos = lines % 4
    
    if (pos == 1) {
        if ($0 !~ /^@/) {
            print fname " - Line " lines ": Must start with '@'" > "/dev/stderr"
            errors++
        }
    }
    else if (pos == 2) {
        if ($0 !~ /^[ATGCNRYSWKMBDHV]+$/) {
            print fname " - Line " lines ": Invalid characters in sequence" > "/dev/stderr"
            errors++
        }
        seq_len = length($0)
        bases_total += seq_len
        
        gc_count = 0
        n_count = 0
        for (i = 1; i <= seq_len; i++) {
            char = substr($0, i, 1)
            if (char == "G" || char == "C") gc_count++
            if (char == "N") n_count++
        }
        gc_total += gc_count
        total_n += n_count
    }
    else if (pos == 3) {
        if ($0 !~ /^\+/) {
            print fname " - Line " lines ": Must have '+' (separator)" > "/dev/stderr"
        }
    }
    else if (pos == 0) {
        if (length($0) != seq_len) {
            print fname " - Line " lines ": Quality length " length($0) " doesn't match sequence " seq_len > "/dev/stderr"
            errors++
        }
        reads++
    }
}
END {
    if (lines % 4 != 0) {
        printf "%s: Corrupt - %d lines is not a multiple of 4\n", fname, lines > "/dev/stderr"
        errors++
    }
    
    gc_percentage = (bases_total > 0) ? (gc_total * 100 / bases_total) : 0
    
    printf "%s: %d reads, %d lines, errors=%d, GC=%.2f%%, N=%d, size=%.2fMB\n",
           fname, reads, lines, errors, gc_percentage, total_n, size
    
    printf "%s,%d,%d,%d,%.2f,%d,%.2f\n",
           fname, reads, lines, errors, gc_percentage, total_n, size >> stats
}
EOF
)
    
    # Execute awk with the script
    awk -v fname="$file_name" -v size="$size_mb" -v stats="$stats_file" "$awk_script" "$file" 2>&1
    
    # Return error count
    local error_count=$(tail -n1 "$stats_file" 2>/dev/null | cut -d',' -f4 || echo "0")
    
    # Send Telegram notification with summary
    local last_line=$(tail -n1 "$stats_file" 2>/dev/null)
    if [[ -n "$last_line" ]]; then
        tg_send "${file_name}: $(echo "$last_line" | cut -d',' -f2-7 | tr ',' ' ')" 2>/dev/null || true
    fi
    
    return ${error_count:-0}
}

####=====================================####
####   Existence and integrity of file   ####
####=====================================####

# Check if input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then 
    echo "ERROR: Directory ${INPUT_DIR} does not exist"
    tg_send "ERROR: Directory ${INPUT_DIR} does not exist" 2>/dev/null || true
    exit 1
else 
    echo "Input directory: ${INPUT_DIR}"
    tg_send "Input directory: ${INPUT_DIR}" 2>/dev/null || true
fi

# Change to input directory
cd "$INPUT_DIR" || exit 1

# Check if there are files to process
if ! ls *.fastq.gz *.fastq 2>/dev/null | grep -q .; then
    echo "No FASTQ files found in ${INPUT_DIR}"
    tg_send "No FASTQ files found in ${INPUT_DIR}" 2>/dev/null || true
    exit 1
fi

echo "Files to process:"
ls -la *.fastq.gz 2>/dev/null || ls -la *.fastq 2>/dev/null

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Initialize counters
total_files=0
valid_files=0
error_files=0

# Stats file for validation results on Lustre
STATS_FILE="${OUTPUT_DIR}/stats_validation.csv"

# Process all FASTQ files
for file in *.fastq.gz *.fastq; do
    [[ ! -f "$file" ]] && continue
    
    ((total_files++))
    file_name=$(basename "$file")
    compressed=0
    
    echo "-----------------------------------------"
    echo "Processing: ${file_name}"
    tg_send "Processing: ${file_name}" 2>/dev/null || true
    
    if [[ -s $file ]]; then
        echo "  File ${file_name} is not empty"
        
        # If file is compressed...
        if [[ "$file_name" == *.gz ]]; then
            echo "  Decompressing ${file_name}..."
            tg_send "  Decompressing ${file_name}..." 2>/dev/null || true
            
            gunzip "$file"
            file="${file%.gz}"
            file_name=$(basename "$file")
            compressed=1
            
            echo "  Decompression done: ${file_name}"
            tg_send "  Decompression done: ${file_name}" 2>/dev/null || true
        fi
        
        # Validate filename pattern
        if [[ $file_name =~ ^PM[0-9]{4}_S[0-9]{1,2}_R[12]\.fastq$ ]]; then
            echo "  File ${file_name} name format is valid"
            
            echo "  Validating content: ${file_name}..."
            if fastq_val "$file"; then
                ((valid_files++))
                echo "   Validation PASSED for ${file_name}"
                tg_send " ${file_name} validation PASSED" 2>/dev/null || true
            else
                ((error_files++))
                echo "   Validation FAILED for ${file_name} with errors"
                tg_send " ${file_name} validation FAILED" 2>/dev/null || true
            fi
        else
            echo "   WARNING: ${file_name} filename pattern is not valid"
            ((error_files++))
        fi
        
        # Re-compress if it was compressed
        if [[ $compressed -eq 1 ]]; then
            echo "  Compressing ${file_name}..."
            tg_send "  Compressing ${file_name}..." 2>/dev/null || true
            gzip "$file"
            echo "  Compression done: ${file_name}.gz"
            tg_send "  Compression done: ${file_name}.gz" 2>/dev/null || true
        fi
        
    else
        echo "   File ${file_name} is empty"
        ((error_files++))
    fi
done

# Summary
echo "========================================="
echo "  FASTQ Validation Summary"
echo "  Date: $(date)"
echo "  Total files processed: ${total_files}"
echo "  Valid files: ${valid_files}"
echo "  Files with errors: ${error_files}"
echo "  Results saved to: ${STATS_FILE}"
echo "========================================="

tg_send "FASTQ Validation Summary: ${total_files} files, ${valid_files} valid, ${error_files} errors" 2>/dev/null || true

# Copy results to permanent storage (NFS)
if [[ -f "${STATS_FILE}" ]]; then
    echo "Copying results to permanent storage..."
    mkdir -p "${PERMANENT_RESULTS}"
    cp "${STATS_FILE}" "${PERMANENT_RESULTS}/"
    echo "Results copied to: ${PERMANENT_RESULTS}/stats_validation.csv"
    tg_send "Results copied to NFS: ${PERMANENT_RESULTS}" 2>/dev/null || true
else
    echo "WARNING: No stats file found to copy to permanent storage"
fi
