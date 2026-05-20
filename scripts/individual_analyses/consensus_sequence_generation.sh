#!/bin/bash

# Script for consensus genome generation from samples that had coverage for viral contigs
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 19/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw/czid_raw/mapping"
CONSENSUS_DIR="${RAW_DATA_DIR}/consensus_genomes"

# Virus names
VIRUS_NAMES=(
    "OSFV"
    "TAV"
    "GMMLV"
)

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

source "$(conda info --base)/etc/profile.d/conda.sh"

echo "=========================================="
echo "Consensus Genome Generation"
echo "=========================================="
echo "Consensus directory: ${CONSENSUS_DIR}"
tg_send "Starting consesus genome..."

####==================================####
####     PROCESS EACH VIRUS TYPE      ####
####==================================####

for VIRUS_NAME in "${VIRUS_NAMES[@]}"; do

    VIRUS_INPUT="${RAW_DATA_DIR}/${VIRUS_NAME}"
    VIRUS_OUTPUT="${CONSENSUS_DIR}/${VIRUS_NAME}"
    

    echo "=========================================="
    echo "Virus: ${VIRUS_NAME}"
    echo "Directory: ${VIRUS_INPUT}"
    echo "=========================================="

    # Check if virus directory exists
    if [ ! -d "$VIRUS_INPUT" ]; then
        echo "WARNING: Directory not found for ${VIRUS_NAME}"
        tg_send "WARNING: No  directory for ${VIRUS_NAME} consesnus genome. Skipping."
        continue
    fi

####==================================####
####   PROCESS EACH SAMPLE DIRECTORY  ####
####==================================####

    for sample_dir in "${VIRUS_INPUT}"/*/; do

        sample_name=$(basename "$sample_dir")
        bam_file="${sample_dir}/${sample_name}_sorted.bam.bai"

        # Skip if BAM file doesn't exist
        if [ ! -f "$bam_file" ]; then
            echo "WARNING: No BAM file for ${sample_name}. Skipping."
            continue
        fi
        
        echo "mpileup generation..."
        
        conda activate samtools_env

        # Generate per-base depth file
        samtools mpileup \
            -r \
            -o ${VIRUS_OUTPUT}/${VIRUS_NAME}

    done

done