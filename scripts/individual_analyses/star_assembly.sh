#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to assemble unmapped reads from STAR
# 07/08/2026
# Ver. 1.0.0 (STAR-only)

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

# Input directory
STAR_INPUT_DIR="${PROJECT_SCRATCH}/results/aligned/STAR_alignment"

# Output directory
OUTPUT_DIR="${PROJECT_SCRATCH}/results/assembly/star_assembly"
FASTQ_DIR="${PROJECT_SCRATCH}/results/assembly/fastq"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
SPADES_CONTAINER="${CONTAINERS}/spades_3.15.5.sif"
MEGAHIT_CONTAINER="${CONTAINERS}/megahit_1.2.9.sif"

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
echo "  STAR Unmapped Reads Assembly Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Input: ${STAR_INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting STAR unmapped reads assembly" 2>/dev/null || true

####===================================####
####     GET SAMPLES FROM STAR         ####
####===================================####

cd "$STAR_INPUT_DIR" || exit 1

# STAR paired samples
star_paired_samples=()
while IFS= read -r line; do
    sample_name=$(basename "$line" "_paired_Aligned.sortedByCoord.out.bam")
    star_paired_samples+=("$sample_name")
done < <(ls -1 *_paired_Aligned.sortedByCoord.out.bam 2>/dev/null | sort)

# STAR R1 unpaired samples
star_r1_unpaired_samples=()
while IFS= read -r line; do
    sample_name=$(basename "$line" "_R1_unpaired_Aligned.sortedByCoord.out.bam")
    star_r1_unpaired_samples+=("$sample_name")
done < <(ls -1 *_R1_unpaired_Aligned.sortedByCoord.out.bam 2>/dev/null | sort)

# STAR R2 unpaired samples
star_r2_unpaired_samples=()
while IFS= read -r line; do
    sample_name=$(basename "$line" "_R2_unpaired_Aligned.sortedByCoord.out.bam")
    star_r2_unpaired_samples+=("$sample_name")
done < <(ls -1 *_R2_unpaired_Aligned.sortedByCoord.out.bam 2>/dev/null | sort)

# Count STAR samples
star_total=$(( ${#star_paired_samples[@]} + ${#star_r1_unpaired_samples[@]} + ${#star_r2_unpaired_samples[@]} ))

if [[ $star_total -eq 0 ]]; then
    echo "No STAR unmapped samples found"
    tg_send "No STAR unmapped samples found" 2>/dev/null || true
    exit 1
fi

echo ""
echo "Summary:"
echo "  - STAR samples: ${star_total}"
echo "    - Paired: ${#star_paired_samples[@]}"
echo "    - R1 unpaired: ${#star_r1_unpaired_samples[@]}"
echo "    - R2 unpaired: ${#star_r2_unpaired_samples[@]}"

####===================================================####
####          GET ARRAY TASK SAMPLE                    ####
####===================================================####

# Combine all samples into a single array with type information
all_samples=()
sample_types=()

# Add paired samples
for sample in "${star_paired_samples[@]}"; do
    all_samples+=("$sample")
    sample_types+=("paired")
done

# Add R1 unpaired samples
for sample in "${star_r1_unpaired_samples[@]}"; do
    all_samples+=("${sample}_R1_unpaired")
    sample_types+=("unpaired_R1")
done

# Add R2 unpaired samples
for sample in "${star_r2_unpaired_samples[@]}"; do
    all_samples+=("${sample}_R2_unpaired")
    sample_types+=("unpaired_R2")
done

# Total samples for array job
total_samples=${#all_samples[@]}

if [[ $total_samples -eq 0 ]]; then
    echo "No STAR samples found"
    tg_send "No STAR samples found" 2>/dev/null || true
    exit 1
fi

# Get the sample for this array task
SAMPLE_INDEX=$((LSB_JOBINDEX - 1))

if [[ $SAMPLE_INDEX -ge $total_samples ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${total_samples})"
    exit 1
fi

sample="${all_samples[$SAMPLE_INDEX]}"
sample_type="${sample_types[$SAMPLE_INDEX]}"

echo "========================================="
echo "  Processing sample: ${sample}"
echo "  Type: ${sample_type}"
echo "  Source: STAR"
echo "  Array task: ${LSB_JOBINDEX}/${total_samples}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${sample} (${sample_type}) from STAR (${LSB_JOBINDEX}/${total_samples})" 2>/dev/null || true

####=====================================####
####        DETERMINE INPUT FILES        ####
####=====================================####

# Create FASTQ directory
mkdir -p "${FASTQ_DIR}/${sample}"

# Determine input files based on sample type
case "$sample_type" in
    "paired")
        R1_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${sample}_paired_Unmapped.out.mate1"
        R2_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${sample}_paired_Unmapped.out.mate2"
        R1_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}_R1.fastq"
        R2_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}_R2.fastq"
        IS_PAIRED=true
        ;;
    "unpaired_R1")
        # Remove _R1_unpaired suffix to get original sample name
        orig_sample="${sample%_R1_unpaired}"
        R1_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${orig_sample}_R1_unpaired_Unmapped.out.mate1"
        R1_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}.fastq"
        IS_PAIRED=false
        ;;
    "unpaired_R2")
        # Remove _R2_unpaired suffix to get original sample name
        orig_sample="${sample%_R2_unpaired}"
        R2_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${orig_sample}_R2_unpaired_Unmapped.out.mate1"
        R2_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}.fastq"
        IS_PAIRED=false
        ;;
esac

# Check if input files exist
if [[ "$IS_PAIRED" == true ]]; then
    if [[ ! -s "$R1_UNMAPPED_ORIG" ]] || [[ ! -s "$R2_UNMAPPED_ORIG" ]]; then
        echo "WARNING: Unmapped files for ${sample} are empty or missing"
        echo "  Skipping assembly for ${sample}"
        exit 0
    fi
    echo "  Copying STAR paired unmapped reads to FASTQ directory..."
    cp -p "$R1_UNMAPPED_ORIG" "$R1_UNMAPPED"
    cp -p "$R2_UNMAPPED_ORIG" "$R2_UNMAPPED"
    gzip -f "$R1_UNMAPPED_ORIG"
    gzip -f "$R2_UNMAPPED_ORIG"
    
    R1_READS=$(($(wc -l < "$R1_UNMAPPED") / 4))
    R2_READS=$(($(wc -l < "$R2_UNMAPPED") / 4))
    echo "  Using STAR unmapped reads (copied to FASTQ directory):"
    echo "  R1: $(basename "$R1_UNMAPPED") ($R1_READS reads)"
    echo "  R2: $(basename "$R2_UNMAPPED") ($R2_READS reads)"
else
    if [[ ! -s "$R1_UNMAPPED_ORIG" ]]; then
        echo "WARNING: Unmapped file for ${sample} is empty or missing"
        echo "  Skipping assembly for ${sample}"
        exit 0
    fi
    echo "  Copying STAR single-end unmapped reads to FASTQ directory..."
    cp -p "$R1_UNMAPPED_ORIG" "$R1_UNMAPPED"
    gzip -f "$R1_UNMAPPED_ORIG"
    
    READS=$(($(wc -l < "$R1_UNMAPPED") / 4))
    echo "  Using STAR unmapped reads (copied to FASTQ directory):"
    echo "  $(basename "$R1_UNMAPPED") ($READS reads)"
fi

INPUT_DIR_FOR_BIND="${FASTQ_DIR}/${sample}"

####================================####
####          RUN ASSEMBLY          ####
####================================####

# Parameters
THREADS=32
MEMORY=128

# Create base output directory with type subdirectory
mkdir -p "${OUTPUT_DIR}/${sample_type}"

# Clean output directories before running
rm -rf "${OUTPUT_DIR}/${sample_type}/rnaSPAdes/${sample}"/*
rm -rf "${OUTPUT_DIR}/${sample_type}/metaSPAdes/${sample}"/*
rm -rf "${OUTPUT_DIR}/${sample_type}/metaviralSPAdes/${sample}"/*
rm -rf "${OUTPUT_DIR}/${sample_type}/MEGAhit/${sample}"/*

# Create output directories
mkdir -p "${OUTPUT_DIR}/${sample_type}/rnaSPAdes/${sample}"
mkdir -p "${OUTPUT_DIR}/${sample_type}/metaSPAdes/${sample}"
mkdir -p "${OUTPUT_DIR}/${sample_type}/metaviralSPAdes/${sample}"
mkdir -p "${OUTPUT_DIR}/${sample_type}/MEGAhit/${sample}"

echo "Starting assembly process for ${sample} (${sample_type}) from STAR..."

# Build input arguments based on type
if [[ "$IS_PAIRED" == true ]]; then
    R1_BASENAME=$(basename "$R1_UNMAPPED")
    R2_BASENAME=$(basename "$R2_UNMAPPED")
    ASSEMBLY_ARGS="-1 /input/${R1_BASENAME} -2 /input/${R2_BASENAME}"
    echo "  Paired-end assembly: ${R1_BASENAME} + ${R2_BASENAME}"
else
    R1_BASENAME=$(basename "$R1_UNMAPPED")
    ASSEMBLY_ARGS="-s /input/${R1_BASENAME}"
    echo "  Single-end assembly: ${R1_BASENAME}"
fi

####================================####
####        rnaSPAdes Assembly      ####
####================================####

echo "Starting rnaSPAdes assembly for sample ${sample} (${sample_type})..."

apptainer exec \
    --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
    --bind "${OUTPUT_DIR}:/output" \
    "$SPADES_CONTAINER" \
    rnaspades.py \
    ${ASSEMBLY_ARGS} \
    -o "/output/${sample_type}/rnaSPAdes/${sample}" \
    -t "$THREADS"

if [[ $? -eq 0 ]]; then
    echo "rnaSPAdes assembly completed for ${sample}"
else
    echo "rnaSPAdes assembly FAILED for ${sample}"
fi

####================================####
####        metaSPAdes Assembly      ####
####================================####

echo "Starting metaSPAdes assembly for sample ${sample} (${sample_type})..."

apptainer exec \
    --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
    --bind "${OUTPUT_DIR}:/output" \
    "$SPADES_CONTAINER" \
    metaspades.py \
    ${ASSEMBLY_ARGS} \
    -o "/output/${sample_type}/metaSPAdes/${sample}" \
    -t "$THREADS"

if [[ $? -eq 0 ]]; then
    echo "metaSPAdes assembly completed for ${sample}"
else
    echo "metaSPAdes assembly FAILED for ${sample}"
fi

####================================####
####     metaviralSPAdes Assembly   ####
####================================####

echo "Starting metaviralSPAdes assembly for sample ${sample} (${sample_type})..."

apptainer exec \
    --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
    --bind "${OUTPUT_DIR}:/output" \
    "$SPADES_CONTAINER" \
    metaviralspades.py \
    ${ASSEMBLY_ARGS} \
    -o "/output/${sample_type}/metaviralSPAdes/${sample}" \
    -t "$THREADS"

if [[ $? -eq 0 ]]; then
    echo "metaviralSPAdes assembly completed for ${sample}"
else
    echo "metaviralSPAdes assembly FAILED for ${sample}"
fi

####================================####
####        MEGAhit Assembly        ####
####================================####

echo "Starting MEGAhit assembly for sample ${sample} (${sample_type})..."

rm -rf "${OUTPUT_DIR}/${sample_type}/MEGAhit/${sample}" 2>/dev/null

apptainer exec \
    --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
    --bind "${OUTPUT_DIR}:/output" \
    "$MEGAHIT_CONTAINER" \
    megahit \
    ${ASSEMBLY_ARGS} \
    -o "/output/${sample_type}/MEGAhit/${sample}" \
    -t "$THREADS" \
    -m "${MEMORY}"

if [[ $? -eq 0 ]]; then
    echo "MEGAhit assembly completed for ${sample}"
else
    echo "MEGAhit assembly FAILED for ${sample}"
fi

####================================####
####        SUMMARY                  ####
####================================####

echo ""
echo "========================================="
echo "  STAR Assembly completed for ${sample} (${sample_type})"
echo "  Output directories:"
echo "    - rnaSPAdes: ${OUTPUT_DIR}/${sample_type}/rnaSPAdes/${sample}"
echo "    - metaSPAdes: ${OUTPUT_DIR}/${sample_type}/metaSPAdes/${sample}"
echo "    - metaviralSPAdes: ${OUTPUT_DIR}/${sample_type}/metaviralSPAdes/${sample}"
echo "    - MEGAhit: ${OUTPUT_DIR}/${sample_type}/MEGAhit/${sample}"
echo "========================================="

tg_send "Completed STAR: ${sample} (${sample_type})" 2>/dev/null || true
