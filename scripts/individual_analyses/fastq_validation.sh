#!/bin/bash

# Authr: Jorge Alberto Castro Rodríguez
# Script to validate fastq files.
# 25/06/2026
# Version 2.2.0

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTDIR="${PROJECT_ROOT}/results"

# Telegram bot.
source "${SCRIPT_DIR}/bot_telegram.sh"

# Variables
directory=$1           # Directory where fastq files are

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


    stats_file=${OUTDIR}/untrimmed_qc/stats/stats_validation.csv

    # Stats 
    mkdir -p "$(dirname "$stats_file")"

    # CSV with headers only one time
    if [[ ! -f "$stats_file" ]]; then
        echo "file,reads,lines,errors,gc_percentage,total_n,size_mb" > "$stats_file"
    fi

    # File size
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
    [[ $bases_total -gt 0 ]] && gc_percentage=$(echo "scale=2; $gc_total * 100 / $bases_total" | bc)
    
    echo "${file_name}: $reads reads, $lines lines, errors=$errors, GC=${gc_percentage}%, N=${total_n}, size=${size_mb}MB"
    tg_send "${file_name}: $reads reads, $lines lines, errors=$errors, GC=${gc_percentage}%, N=${total_n}, size=${size_mb}MB"

    # Write on each CSV row
    echo "${file_name},${reads},${lines},${errors},${gc_percenatge},${total_n},${size_mb}" >> "$stats_file"
    return $errors
}

####=====================================####
####   Existence and integrity of file   ####
####=====================================####

if [[ ! -d $directory ]]; then 
    echo "Directory ${directory} does not exist"
    tg_send "Directory ${directory} does not exist"
    exit 1
else 
    echo "Directory ${directory} already exists"
    tg_send "Directorio ${directory} already exists"
fi

# Verification for all fastq files
for file in "$directory"/*.fastq "$directory"/*.fastq.gz; do
    [[ ! -f "$file" ]] && continue  # Ignores if it is not a file
    
    file_name=$(basename "$file")
    
    echo "File ${file_name} exists in ${directory}"
    tg_send "File ${file_name} exists in ${directory}"

    if [[ -s $file ]]; then 
        echo "File ${file_name} is not empty"
        tg_send "File ${file_name} is not empty"

        # If file is compressed...
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

       # File compression
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