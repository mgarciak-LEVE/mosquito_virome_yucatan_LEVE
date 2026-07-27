#!/bin/bash
#BSUB -o fastq_validate_%J.out
#BSUB -e fastq_validate_%J.err
#BSUB -q normal
#BSUB -n 1                     # single core
#BSUB -M 2000                  # 2 GB memory
#BSUB -R "select[mem>2000] rusage[mem=2000] span[hosts=1]"
#BSUB -G team222          

# Author: Jorge Alberto Castro Rodríguez (adapted for Sanger farm)
# Script to validate fastq files.
# 25/06/2026
# Version 2.2.0 (farm‑ready)

####==================================####
####           CONFIGURATION          ####
####==================================####

# --- Permanent storage locations ---
# Set these to your actual team directories
PERM_INPUT_BASE="/nfs/team222/input"    # where your original fastq dirs live
PERM_OUTPUT_BASE="/nfs/team222/results" # where final results should go

# --- Scratch location ---
SCRATCH_BASE="/lustre/scratch126/${USER}"   # use your scratch filesystem
JOB_SCRATCH="${SCRATCH_BASE}/fastq_val_${LSB_JOBID}"

# --- Input directory from command line (relative to permanent storage) ---
INPUT_DIR_NAME="$1"
if [[ -z "$INPUT_DIR_NAME" ]]; then
    echo "ERROR: Please provide a directory name"
    exit 1
fi
PERM_INPUT_DIR="${PERM_INPUT_BASE}/${INPUT_DIR_NAME}"

# --- Script directory and Telegram bot  ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/bot_telegram.sh" 

# --- Output directory on Lustre (temporary) ---
OUTDIR="${JOB_SCRATCH}/results"

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

    stats_file="${OUTDIR}/untrimmed_qc/stats/stats_validation.csv"

    mkdir -p "$(dirname "$stats_file")"

    if [[ ! -f "$stats_file" ]]; then
        echo "file,reads,lines,errors,gc_percentage,total_n,size_mb" > "$stats_file"
    fi

    local bytes=$(stat -c%s "$file")
    local size_mb=$(echo "scale=2; $bytes / 1048576" | bc)

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
                
                gc_count=$(echo "$line" | tr -cd 'GC' | wc -c)
                ((gc_total += gc_count))
                
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
    
    if [[ $((lines % 4)) -ne 0 ]]; then
        echo "${file_name}: Corrupt - $lines lines (is not a multiple of 4)"
        ((errors++))
    fi

    local gc_percentage=0
    [[ $bases_total -gt 0 ]] && gc_percentage=$(echo "scale=2; $gc_total * 100 / $bases_total" | bc)
    
    echo "${file_name}: $reads reads, $lines lines, errors=$errors, GC=${gc_percentage}%, N=${total_n}, size=${size_mb}MB"
    tg_send "${file_name}: $reads reads, $lines lines, errors=$errors, GC=${gc_percentage}%, N=${total_n}, size=${size_mb}MB"

    echo "${file_name},${reads},${lines},${errors},${gc_percentage},${total_n},${size_mb}" >> "$stats_file"
    return $errors
}

####=====================================####
####   Existence and integrity of file   ####
####=====================================####

# --- Step 1: Copy input data from permanent storage to Lustre ---
echo "Copying input from ${PERM_INPUT_DIR} to Lustre scratch ${JOB_SCRATCH}..."
mkdir -p "${JOB_SCRATCH}"
rsync -av "${PERM_INPUT_DIR}/" "${JOB_SCRATCH}/input/"
if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to copy input data."
    tg_send "ERROR: Failed to copy input data from ${PERM_INPUT_DIR}"
    exit 1
fi
tg_send "Input copied to Lustre scratch."

# --- Step 2: Work on Lustre ---
cd "${JOB_SCRATCH}/input" || { echo "Cannot cd to scratch"; exit 1; }

if [[ ! -d . ]]; then 
    echo "Directory ${JOB_SCRATCH}/input does not exist"
    tg_send "Directory ${JOB_SCRATCH}/input does not exist"
    exit 1
else 
    echo "Working directory: ${JOB_SCRATCH}/input"
fi

# Verification for all fastq files
for file in *.fastq *.fastq.gz; do
    [[ ! -f "$file" ]] && continue
    
    file_name=$(basename "$file")
    compressed=0
    
    echo "File ${file_name} exists in scratch"
    tg_send "File ${file_name} exists in scratch"

    if [[ -s $file ]]; then 
        echo "File ${file_name} is not empty"
        tg_send "File ${file_name} is not empty"

        if [[ "$file_name" == *.gz ]]; then
            echo "Decompressing ${file_name}..."
            tg_send "Decompressing ${file_name}..."
            gunzip "$file"
            file="${file%.gz}"
            file_name=$(basename "$file")
            compressed=1
            echo "Decompression done: ${file_name}"
            tg_send "Decompression done: ${file_name}"
        fi

        if [[ $file_name =~ ^PM[0-9]{4}_S[0-9]{1,2}_R[12]\.fastq$ ]]; then
            echo "File ${file_name} is valid"
            tg_send "${file_name} is valid"
            echo ""
            echo "Validating content: ${file_name}..."
            fastq_val "$file"
            echo ""
        else
            echo "Error: ${file_name} is not valid"
            tg_send "${file_name} is not valid"
        fi

        if [[ $compressed -eq 1 ]]; then
            echo "Compressing ${file_name}..."
            tg_send "Compressing ${file_name}..."
            gzip "$file"
            echo "Compression done: ${file_name}.gz"
            tg_send "Compression done ${file_name}.gz"
        fi

    else 
        echo "File ${file_name} is empty"
        tg_send "${file_name} is empty"
    fi
done

# --- Copy important results back to permanent storage ---
echo "Copying results to permanent storage: ${PERM_OUTPUT_BASE}/${INPUT_DIR_NAME}/"
mkdir -p "${PERM_OUTPUT_BASE}/${INPUT_DIR_NAME}/untrimmed_qc/stats/"
cp -r "${OUTDIR}/untrimmed_qc/stats/" "${PERM_OUTPUT_BASE}/${INPUT_DIR_NAME}/untrimmed_qc/stats/"
# Also copy any other outputs you want to keep

# --- Clean up Lustre scratch ---
echo "Cleaning up scratch directory ${JOB_SCRATCH}"
cd /
rm -rf "${JOB_SCRATCH}"

echo "Job finished successfully."
tg_send "Validation job for ${INPUT_DIR_NAME} completed."