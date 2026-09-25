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
CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes"
CONTIG_BASE="${PROJECT_ROOT}/data/raw/czid_raw/viral_contigs"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
SAMTOOLS_SIF="${CONTAINER_DIR}/samtools_1.20--h50ea8bc_1.sif"
IVAR_SIF="${CONTAINER_DIR}/ivar_1.4.4--h077b44d_0.sif"

# Consensus parameters
MIN_DEPTH=3        # Minimum depth to call a base
ALLELE_FREQ=0.5    # Minimum allele frequency (50%)

# Virus names
VIRUS_NAMES=(
    "OSFV"
    "TAV"
    "GMMLV"
)

VIRUS_CONTIG_SUBDIRS=(
    "mosquito_viromes_leve_OSFV"
    "mosquito_viromes_leve_TAV"
    "mosquito_viromes_leve_GMMLV"
)

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

echo "=========================================="
echo "Consensus Genome Generation"
echo "=========================================="
echo "Consensus directory: ${CONSENSUS_DIR}"
tg_send "Starting consensus genome..."

####==================================####
####     PROCESS EACH VIRUS TYPE      ####
####==================================####

for i in "${!VIRUS_NAMES[@]}"; do

    VIRUS_NAME="${VIRUS_NAMES[$i]}"
    CONTIG_SUBDIR="${VIRUS_CONTIG_SUBDIRS[$i]}"
    
    VIRUS_INPUT="${RAW_DATA_DIR}/${VIRUS_NAME}"
    VIRUS_OUTPUT="${CONSENSUS_DIR}/${VIRUS_NAME}"
    VIRUS_CONTIG_DIR="${CONTIG_BASE}/${CONTIG_SUBDIR}"

    echo "=========================================="
    echo "Virus: ${VIRUS_NAME}"
    echo "Input:       ${VIRUS_INPUT}"
    echo "Output:      ${VIRUS_OUTPUT}"
    echo "Contig dir:  ${VIRUS_CONTIG_DIR}"
    echo "=========================================="

    # Check if virus directory exists
    if [ ! -d "$VIRUS_INPUT" ]; then
        echo "WARNING: Directory not found for ${VIRUS_NAME}"
        tg_send "WARNING: No  directory for ${VIRUS_NAME}. Skipping."
        continue
    fi

    # Create output directory
    mkdir -p "${VIRUS_OUTPUT}"

####==================================####
####   PROCESS EACH SAMPLE DIRECTORY  ####
####==================================####

    for sample_dir in "${VIRUS_INPUT}"/*/; do

        sample_name=$(basename "$sample_dir")
        bam_file="${sample_dir}/${sample_name}_sorted.bam"

        # Skip if BAM file doesn't exist
        if [ ! -f "$bam_file" ]; then
            echo "WARNING: No BAM file for ${sample_name}. Skipping."
            continue
        fi

        # Find the original contig file
        contig_file=$(ls "${VIRUS_CONTIG_DIR}/${sample_name}"* 2>/dev/null | head -1)

        if [ -z "$contig_file" ] || [ ! -f "$contig_file" ]; then
            echo "  WARNING: No contig FASTA found for ${sample_name}. Skipping."
            tg_send "WARNING: No contig for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        if [ ! -s "$contig_file" ]; then
            echo "  WARNING: Contig file is empty for ${sample_name}. Skipping."
            tg_send "WARNING: Empty contig for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        echo "----------------------------------------"
        echo "Sample: ${sample_name}"
        echo "BAM:    $(basename "$bam_file")"
        echo "Contig: $(basename "$contig_file")"
        echo "----------------------------------------"
        
        ####============================####
        ####      MPILEUP + CONSENSUS   ####
        ####============================####

        echo "Generating mpileup and calling consensus..."
        tg_send "Generating consensus for ${sample_name}..."

        apptainer exec "${SAMTOOLS_SIF}" samtools mpileup \
            -aa \
            -A \
            -d 0 \
            -f "$contig_file" \
            "$bam_file" \
        | apptainer exec "${IVAR_SIF}" ivar consensus \
            -p "${VIRUS_OUTPUT}/${sample_name}_consensus" \
            -m "${MIN_DEPTH}" \
            -t "${ALLELE_FREQ}" \
            -n N

        if [ $? -ne 0 ]; then
            echo "  ERROR: Consensus failed for ${sample_name}"
            tg_send "ERROR: Consensus failed for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        echo "  Consensus: ${sample_name}_consensus.fa"
        tg_send "Consensus: ${VIRUS_NAME}/${sample_name}"

    done

    echo "Completed ${VIRUS_NAME}"
    tg_send "Completed consensus for ${VIRUS_NAME}"

done

echo "=========================================="
echo "Consensus Generation Complete"
echo "=========================================="
echo ""
echo "Output: ${CONSENSUS_DIR}/<virus>/<sample>_consensus.fa"

tg_send "Consensus generation complete."