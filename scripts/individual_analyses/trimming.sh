#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to validate fastq files.
# 15/07/2026
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

# TRIMMING FUNCTION
run_trimming() {
    #Runs trimming for paired-end reads. It distinguishes between the R1 and R2 reads
    echo "Starting trimming process..."
    mkdir -p "${TRIMMED_DIR}" # -p flag means to create parent directories if they dont exist 

    for R1_file in "$WORKDIR"/*_R1.fastq; do   # Does a loop. Finds every file in the working directory that matches the pattern "*_R1.fastq" and assigns each file path to the variable R1
        R2_file="${R1_file/_R1.fastq/_R2.fastq}" # Parameter expansion: ${variable/pattern/replacement}
    
    newR1=$(sample_names "$R1_file" "full")
    newR2=$(sample_names "$R2_file" "full")


    conda activate trimmomatic 
    echo "Input files: $(basename "$R1_file") and $(basename "$R2_file")"
    echo "Processing: $newR1 and $newR2"
    adapters="/Users/Parsimony/miniconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa" 
    
    # For the forward and reverse files ($R1 and $R2)
    trimmomatic PE -threads 12 -phred33 \
        "$R1_file" "$R2_file" \
        "${TRIMMED_DIR}/${newR1}_paired.fastq" \
        "${TRIMMED_DIR}/${newR1}_unpaired.fastq" \
        "${TRIMMED_DIR}/${newR2}_paired.fastq" \
        "${TRIMMED_DIR}/${newR2}_unpaired.fastq" \
        ILLUMINACLIP:${adapters}:2:30:10:2:keepBothReads \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:3:25 MINLEN:50 # ILUMMINACLIP is for adapter clipping
    done   
    conda deactivate
}

####=================================####
####        FUNCTION EXECUTION       ####
####=================================####

# Create output directories
mkdir -p "${OUTPUT_FASTQC}"
mkdir -p "${OUTPUT_MULTIQC}"

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Script is being executed directly
    echo "Sequence quality inspection module"
    
    OUTPUT_FASTQC="${OUTPUT_BASE}/untrimmed_qc/fastqc"
    OUTPUT_MULTIQC="${OUTPUT_BASE}/untrimmed_qc/multiqc"

    # Run the pipeline steps
    run_fastqc "$INPUT_FASTQC" "$OUTPUT_FASTQC"
    run_multiqc "$INPUT_MULTIQC" "$OUTPUT_MULTIQC"

fi

echo ""
echo "========================================="
echo "  Copying results to permanent storage"
echo "========================================="

# Copy FastQC results
if [[ -d "$OUTPUT_FASTQC" ]]; then
    echo "Copying FastQC results..."
    mkdir -p "${PERMANENT_RESULTS}/fastqc"
    cp -r "${OUTPUT_FASTQC}"/* "${PERMANENT_RESULTS}/fastqc/" 2>/dev/null || true
    echo "FastQC results copied to: ${PERMANENT_RESULTS}/fastqc/"
fi

# Copy MultiQC results
if [[ -d "$OUTPUT_MULTIQC" ]]; then
    echo "Copying MultiQC results..."
    mkdir -p "${PERMANENT_RESULTS}/multiqc"
    cp -r "${OUTPUT_MULTIQC}"/* "${PERMANENT_RESULTS}/multiqc/" 2>/dev/null || true
    echo "MultiQC results copied to: ${PERMANENT_RESULTS}/multiqc/"
fi

tg_send "Sequence Quality module completed for ${PROJECT_NAME}" 2>/dev/null || true

echo ""
echo "========================================="
echo "  Sequence Quality Module Complete"
echo "  FastQC: ${OUTPUT_FASTQC}"
echo "  MultiQC: ${OUTPUT_MULTIQC}"
echo "  Permanent storage: ${PERMANENT_RESULTS}"
echo "========================================="
