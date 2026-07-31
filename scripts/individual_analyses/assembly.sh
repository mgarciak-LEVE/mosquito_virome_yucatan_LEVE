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

####===================================================####
####          GET SAMPLE LIST FOR STAR ASSEMBLY        ####
####===================================================####

cd "$STAR_INPUT_DIR" || exit 1

# Get list of samples
samples=()
while IFS= read -r line; do
    samples+=("$line")
done < <(find . -maxdepth 1 -type d -name "PM[0-9]{4}_paired_Unmapped.out.mate[12]" | sed 's/\.\///' | sort)

# Check if samples exist
if [[ ${#samples[@]} -eq 0 ]]; then
    echo "No samples found in ${STAR_INPUT_DIR}"
    tg_send "No samples found" 2>/dev/null || true
    exit 1
fi

echo "Found ${#samples[@]} samples"

# Get the sample for this array task
SAMPLE_INDEX=$((LSB_JOBINDEX - 1))

if [[ $SAMPLE_INDEX -ge ${#samples[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#samples[@]})"
    exit 1
fi

sample="${samples[$SAMPLE_INDEX]}"

echo "========================================="
echo "  Sequence Assembly Array ${LSB_JOBINDEX}/${#samples[@]}"
echo "  Sample: ${sample}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Assembly: ${sample} (${LSB_JOBINDEX}/${#samples[@]})" 2>/dev/null || true


####================================####
####     DECOMPRESS INPUT FILES      ####
####================================####

echo "Checking for compressed input files..."

# Decompress R1 paired if compressed
if [[ -f "${STAR_INPUT_DIR}/${sample}_Unmapped.out.mate1.gz" ]]; then
    echo "Decompressing: ${sample}_Unmapped.out.mate1.gz"
    gzip -d "${INPUT_DIR}/${sample}_Unmapped.out.mate1.gz"
fi

# Decompress R2 paired if compressed
if [[ -f "${STAR_INPUT_DIR}/${sample}_Unmapped.out.mate2.gz" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq" ]]; then
    echo "Decompressing: ${sample}_R2_paired.fastq.gz"
    gzip -d "${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq.gz"
fi

# Decompress R1 unpaired if compressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R1_unpaired.fastq.gz" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R1_unpaired.fastq" ]]; then
    echo "Decompressing: ${sample}_R1_unpaired.fastq.gz"
    gzip -d "${INPUT_DIR}/${sample}/${sample}_R1_unpaired.fastq.gz"
fi

# Decompress R2 unpaired if compressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R2_unpaired.fastq.gz" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R2_unpaired.fastq" ]]; then
    echo "Decompressing: ${sample}_R2_unpaired.fastq.gz"
    gzip -d "${INPUT_DIR}/${sample}/${sample}_R2_unpaired.fastq.gz"
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

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
