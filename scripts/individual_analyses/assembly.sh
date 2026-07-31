#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to assemble 
# 30/07/2026
# Ver. 1.1.0 (farm-ready)

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directories
STAR_INPUT_DIR="${PROJECT_SCRATCH}/results/aligned/STAR_alignment"
BOWTIE_INPUT_DIR="${PROJECT_SCRATCH}/results/aligned/Bowtie2_alignment"

# Output directory
OUTPUT_DIR="${PROJECT_SCRATCH}/assembly"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Unmapped Reads Assembly Module"
echo "  Project: ${PROJECT_NAME}"
echo "  STAR unmapped reads input: ${STAR_INPUT_DIR}"
echo "  Bowtie unmapped reads input: ${BOWTIE_INPUT_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting unampped reads assembly" 2>/dev/null || true

# Create output directory
mkdir -p "${OUTPUT_DIR}/statistics"

    echo "Starting assembly process..."
    mkdir -p "${ASSEMBLY_DIR}/fastq" # Output directory for assembly

    for R1 in "${ALIGNED_DIR}"/*_Unmapped.out.mate1; do
        R2="${R1/_Unmapped.out.mate1/_Unmapped.out.mate2}"
        sample=$(basename $R1 | cut -d'_' -f1)

        # Sample specific output directories
        mkdir -p "${ASSEMBLY_DIR}/rnaSPAdes/${sample}"
        mkdir -p "${ASSEMBLY_DIR}/metaSPAdes/${sample}"

        echo "Cleaning up previous MEGAhit results for $sample..."
        rm -rf "$OUTDIR/assembly/MEGAhit/$sample" 2>/dev/null
        

        R1_fastq="${ASSEMBLY_DIR}/fastq/${sample}_R1.fastq"
        R2_fastq="${ASSEMBLY_DIR}/fastq/${sample}_R2.fastq"

        echo "Extracting unmapped reads from BAM format to FASTQ"
        
        conda activate samtools_env
        samtools fastq "$R1" > "$R1_fastq"
        samtools fastq "$R2" > "$R2_fastq"
        conda deactivate

        # rnaSPAdes --------
        # rnaSPAdes is good for transcriptome analyses.
        echo "Starting assembly through rnaSPAdes..."
        # Input data as every corresponding unmapped read file.

        echo "Starting rnaSPAdes assembly for sample $sample..."
        conda activate SPADES_env
        rnaspades.py \
        -1 "$R1_fastq" \
        -2 "$R2_fastq" \
        -o "$OUTDIR"/assembly/rnaSPAdes/$sample \
        -t $threads
        conda deactivate
        echo "rnaSPAdes assembly finished for sample $sample"

        # metaSPAdes. --------
        # metaSPAdes is good for small genomes.
        echo "Starting assembly through metaSPAdes..." 

        # Input files were the same unmapped reads with the same format. 
        echo "Starting metaSPAdes assembly for sample $sample..."
        conda activate SPADES_env
        metaspades.py \
        -1 "$R1_fastq" \
        -2 "$R2_fastq" \
        -o "$OUTDIR"/assembly/metaSPAdes/$sample \
        -t $threads
        conda deactivate
        echo "metaSPAdes assembly finished for sample $sample"

        # MEGAhit. --------
        echo "Starting assembly through MEGAhit..."

        echo "Starting MEGAhit assembly for sample $sample..."
        conda activate MEGAhit_env 
        megahit \
        -1 "$R1_fastq" \
        -2 "$R2_fastq" \
        -o "$OUTDIR"/assembly/MEGAhit/$sample \
        -t $threads \
        -m 30
        conda deactivate
        echo "MEGAhit assembly finished for sample $sample"

    done

    echo "All assembly processes completed"
}
