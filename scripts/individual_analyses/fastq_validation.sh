#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to validate fastq files.
# 02/07/2026
# Version 2.2.1 (farm-ready)

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_pipeline}"

# --- STORAGE LOCATIONS ---
# Permanent storage
PERMANENT_BASE="/nfs/team222/projects"

# Scratch storage
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46"

# --- PROJECT DIRECTORIES ---
# Working directory on Lustre
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directory - where FASTQ files are
# Uses symlink by default, or can be overridden with $2
INPUT_DIR="${2:-${HOME}/git_repos/mosquito_virome_pipeline/raw_data}"

# Output directory for validation results (on Lustre)
OUTDIR="${PROJECT_SCRATCH}/results/untrimmed_qc"

# Permanent results directory (backed up - copy final results here)
PERMANENT_RESULTS="${PERMANENT_BASE}/${PROJECT_NAME}/results"

# Scripts directory (on NFS)
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot (source from NFS)
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

# Print configuration
echo "========================================="
echo "  FASTQ Validation"
echo "  Project: ${PROJECT_NAME}"
echo "  Input directory: ${INPUT_DIR}"
echo "  Output directory: ${OUTDIR}"
echo "  Date: $(date)"
echo "========================================="

####==================================####
####     FASTQ INTEGRITY FUNCTION     ####
####==================================####

fastq_val() {
    local file="$1"
    local file_name=$(basename "$file")
    local lines=0
    local reads=0
    local errors=0
    local gc_total=0
    local bases_total=0
    local total_n=0

    stats_file="${OUTDIR}/stats/stats_validation.csv"

    # Stats 
    mkdir -p "$(dirname "$stats_file")"

    # CSV with headers only one time
    if [[ ! -f "$stats_file" ]]; then
        echo "file,reads,lines,errors,gc_percentage,total_n,size_mb" > "$stats_file"
    fi

    # File size
    local bytes=$(stat -c%s "$file" 2>/dev/null || echo "0")
    local size_mb=$(echo "scale=2; $bytes / 1048576" | bc 2>/dev/null || echo "0")

    while IFS= read -r line; do
        ((lines++))
        local pos=$((lines % 4))
        
        case $pos in
            1)  # ID line
                if [[ ! "$line" =~ ^@ ]]; then
                    echo "${file_name} - Line $lines: Must start with '@'"
                    ((errors++))
                fi
                ;;
            2)  # Sequence line
                if [[ ! "$line" =~ ^[ATGCNRYSWKMBDHV]+$ ]]; then
                    echo "${file_name} - Line $lines: No valid characters on sequence"
                    ((errors++))
                fi
                seq_len=${#line}
                ((bases_total += seq_len))
                
                # GC content
                gc_count=$(echo "$line" | tr -cd 'GC' | wc -c)
                ((gc_total += gc_count))
                
                # Amount of N's
                n_count=$(echo "$line" | tr -cd 'N' | wc -c)
                ((total_n += n_count))
                ;;
            3)  # Separator
                if [[ ! "$line" =~ ^\+ ]]; then
                    echo "${file_name} - Line $lines: Must have '+' (separator)"
                fi
                ;;
            0)  # Quality line
                if [[ ${#line} -ne $seq_len ]]; then
                    echo "${file_name} -Line $lines: Quality (${#line}) doesn't match sequence length ($seq_len)"
                    ((errors++))
                fi
                ((reads++))
                ;;
        esac
    done < "$file"
    
    # Final validation
    if [[ $((lines % 4)) -ne 0 ]]; then
        echo "${file_name}: Corrupt - $lines lines (is not a multiple of 4)"
        ((errors++))
    fi

    # GC content 
    local gc_percentage=0
    [[ $bases_total -gt 0 ]] && gc_percentage=$(echo "scale=2; $gc_total * 100 / $bases_total" | bc 2>/dev/null || echo "0")
    
    echo "${file_name}: $reads reads, $lines lines, errors=$errors, GC=${gc_percentage}%, N=${total_n}, size=${size_mb}MB"
    tg_send "${file_name}: $reads reads, $lines lines, errors=$errors, GC=${gc_percentage}%, N=${total_n}, size=${size_mb}MB" 2>/dev/null || true

    # Fixed: gc_percenatge -> gc_percentage
    echo "${file_name},${reads},${lines},${errors},${gc_percentage},${total_n},${size_mb}" >> "$stats_file"
    return $errors
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
    tg_send "Starting FASTQ validation for ${INPUT_DIR}" 2>/dev/null || true
fi

# Create output directories
mkdir -p "${OUTDIR}/fastqc"
mkdir -p "${OUTDIR}/multiqc"
mkdir -p "${OUTDIR}/stats"

# Initialize counters
total_files=0
valid_files=0
error_files=0

# Verification for all fastq files
for file in "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fastq.gz; do
    [[ ! -f "$file" ]] && continue  # Ignores if it is not a file
    
    ((total_files++))
    file_name=$(basename "$file")
    compressed=0  
    
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
            fastq_val "$file"
            if [[ $? -eq 0 ]]; then
                ((valid_files++))
                echo "  Validation PASSED for ${file_name}"
                tg_send "${file_name} validation PASSED" 2>/dev/null || true
            else
                ((error_files++))
                echo "Validation FAILED for ${file_name} with errors"
                tg_send "${file_name} validation FAILED" 2>/dev/null || true
            fi
        else
            echo "WARNING: ${file_name} filename pattern is not valid"
            ((error_files++))
        fi

        if [[ $compressed -eq 1 ]]; then
            echo "  Compressing ${file_name}..."
            tg_send "  Compressing ${file_name}..." 2>/dev/null || true
            gzip "$file"
            echo "  Compression done: ${file_name}.gz"
            tg_send "  Compression done: ${file_name}.gz" 2>/dev/null || true
        fi

    else 
        echo "File ${file_name} is empty"
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
echo "  Results saved to: ${OUTDIR}/stats/stats_validation.csv"
echo "========================================="

tg_send "FASTQ Validation Summary: ${total_files} files, ${valid_files} valid, ${error_files} errors" 2>/dev/null || true

# Copy results to permanent storage (NFS)
if [[ -f "${OUTDIR}/stats/stats_validation.csv" ]]; then
    echo "Copying results to permanent storage..."
    mkdir -p "${PERMANENT_RESULTS}/untrimmed_qc/stats/"
    cp "${OUTDIR}/stats/stats_validation.csv" "${PERMANENT_RESULTS}/untrimmed_qc/stats/"
    echo "Results copied to: ${PERMANENT_RESULTS}/untrimmed_qc/stats/"
fi

echo ""
echo "   - Input data: ${INPUT_DIR}"
echo "   - Results (Lustre): ${OUTDIR}"
echo "   - Permanent storage (NFS): ${PERMANENT_RESULTS}"
echo "   - LSF logs: ~/lsf_logs/$(date +%d_%m_%Y)/"
echo ""
