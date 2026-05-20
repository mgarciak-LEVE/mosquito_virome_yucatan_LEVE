#!/bin/bash

# Script for non-host read mapping to viral contigs for validation.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 18/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw/czid_raw"  

READS_DIR="${RAW_DATA_DIR}/non_host_reads"
OUTPUT_DIR="${PROJECT_ROOT}/data/raw/czid_raw/mapping"

# Virus-specific directories
VIRUS1_CONTIG_DIR="${RAW_DATA_DIR}/non_host_contigs/mosquito_viromes_leve_OSFV"
VIRUS2_CONTIG_DIR="${RAW_DATA_DIR}/non_host_contigs/mosquito_viromes_leve_TAV"
VIRUS3_CONTIG_DIR="${RAW_DATA_DIR}/non_host_contigs/mosquito_viromes_leve_GMMLV"

# Array of virus directories for iteration
VIRUS_DIRS=(
    "${VIRUS1_CONTIG_DIR}"
    "${VIRUS2_CONTIG_DIR}"
    "${VIRUS3_CONTIG_DIR}"
)

# Virus names (for output organization and logging)
VIRUS_NAMES=(
    "OSFV"
    "TAV"
    "GMMLV"
)

THREADS=14

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

source "$(conda info --base)/etc/profile.d/conda.sh" 

# Create output directory
mkdir -p "${OUTPUT_DIR}"

cd "${RAW_DATA_DIR}" || exit 1

echo "Directory: $(pwd)"
tg_send "Directory: $(pwd)"

echo "Output directory: ${OUTPUT_DIR}"
tg_send "Output directory: ${OUTPUT_DIR}"

####==================================####
####     PROCESS EACH VIRUS TYPE      ####
####==================================####

for i in "${!VIRUS_DIRS[@]}"; do

    VIRUS_DIR="${VIRUS_DIRS[$i]}"
    VIRUS_NAME="${VIRUS_NAMES[$i]}"
    VIRUS_OUTPUT="${OUTPUT_DIR}/${VIRUS_NAME}"

    echo "Virus: ${VIRUS_NAME}"
    echo "Contig directory: ${VIRUS_DIR}"
    echo "Output directory: ${VIRUS_OUTPUT}"

    # Check if virus directory exists
    if [ ! -d "$VIRUS_DIR" ]; then
        echo "WARNING: Directory not found for ${VIRUS_NAME}"
        tg_send "WARNING: Directory not found for ${VIRUS_NAME}. Skipping."
        continue
    fi

    # Create output directory for this virus
    mkdir -p "${VIRUS_OUTPUT}"

    # Count contig files
    contig_count=$(ls "${VIRUS_DIR}"/*.fasta 2>/dev/null | wc -l)
    
    if [ "$contig_count" -eq 0 ]; then
        echo "No contig files (.fasta) found for ${VIRUS_NAME}"
        tg_send "WARNING: No contig files for ${VIRUS_NAME}"
        continue
    fi

    echo "Found ${contig_count} contig file(s) for ${VIRUS_NAME}"

####==================================####
####   PROCESS EACH CONTIG (SAMPLE)   ####
####==================================####

    for contig_file in "${VIRUS_DIR}"/*.fasta; do

        # Skip if no files match the pattern
        if [ ! -f "$contig_file" ]; then
            continue
        fi

        full_filename=$(basename "$contig_file" .fasta)

        # Extract sample name from contig filename
        sample_name=$(echo "$full_filename" | sed 's/_contigs.*$//')

        echo "Sample: ${sample_name}"
        echo "Contig: ${full_filename}.fasta"

        # Create sample-specific output directory
        sample_out="${VIRUS_OUTPUT}/${sample_name}"
        mkdir -p "${sample_out}"

        ####============================####
        ####       CONTIG INDEXING       ####
        ####============================####

        echo "Indexing contig..."
        
        conda activate minimap2_env

        minimap2 -d "${sample_out}/${sample_name}_contigs.mmi" "$contig_file"

        if [ $? -eq 0 ]; then
            echo "Index created: ${sample_name}_contigs.mmi"
            tg_send "Index created: ${sample_name}_contigs.mmi"

            # Optional: get contig stats
            contig_size=$(wc -c < "$contig_file" | tr -d ' ')
            contig_seqs=$(grep -c "^>" "$contig_file")
            echo "  Contig size: ${contig_size} bytes, ${contig_seqs} sequence(s)"
            
        else
            echo " ERROR: Indexing failed for ${sample_name}"
            tg_send "ERROR: Indexing failed for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        ####============================####
        ####    FIND NON-HOST READS     ####
        ####============================####

        echo "Finding non-host reads..."

        reads_file=""
        reads_file=$(ls "${READS_DIR}/${sample_name}"*.fasta.gz 2>/dev/null | head -1)
        if [ -z "$reads_file" ]; then
            reads_file=$(ls "${READS_DIR}/${sample_name}"*.fasta 2>/dev/null | head -1)
        fi

        if [ -z "$reads_file" ]; then
            echo "  ERROR: No read file found for ${sample_name}"
            tg_send "ERROR: No reads for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        echo "  Read file: $(basename "$reads_file")"

        ####============================####
        ####        READ MAPPING         ####
        ####============================####

        echo "Mapping reads to contig..."

        minimap2 -ax sr \
            -t "${THREADS}" \
            "${sample_out}/${sample_name}_contigs.mmi" \
            "$reads_file" \
            > "${sample_out}/${sample_name}_aligned.sam"

        if [ $? -ne 0 ]; then
            echo "  ERROR: Mapping failed for ${sample_name}"
            tg_send "ERROR: Mapping failed for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        conda deactivate

        conda activate samtools_env

        # Count mapped reads
        mapped_reads=$(samtools view -c "${sample_out}/${sample_name}_aligned.sam" 2>/dev/null)
        echo "  Reads mapped: ${mapped_reads}"
        tg_send "Mapped: ${VIRUS_NAME}/${sample_name} (${mapped_reads} reads)"

        conda deactivate

        ####============================####
        ####      SAM TO BAM + INDEX     ####
        ####============================####

        echo "Converting to sorted BAM..."

        conda activate samtools_env

        samtools sort \
            -@ "${THREADS}" \
            "${sample_out}/${sample_name}_aligned.sam" \
            -o "${sample_out}/${sample_name}_sorted.bam"

        if [ $? -ne 0 ]; then
            echo "  ERROR: SAM to BAM conversion failed"
            tg_send "ERROR: BAM conversion failed for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        samtools index "${sample_out}/${sample_name}_sorted.bam"

        if [ $? -ne 0 ]; then
            echo "  ERROR: BAM indexing failed"
            tg_send "ERROR: BAM indexing failed for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        conda deactivate
        echo "BAM created and indexed"

        ####============================####
        ####          CLEANUP            ####
        ####============================####

        rm -f "${sample_out}/${sample_name}_aligned.sam"
        echo "SAM file removed"


    done

    echo "Completed ${VIRUS_NAME}: ${contig_count} contig file(s) processed"
    tg_send "Completed ${VIRUS_NAME}: ${contig_count} file(s)"

done

echo "Mapping pipeline complete."
tg_send "Mapping pipeline complete."