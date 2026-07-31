#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to map trimmed sequences to mosquito genome superreference
# 28/07/2026
# Ver. 1.0.3 (farm-ready with compression support)

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

# STAR index directory
STAR_SUPER_REF="${PROJECT_SCRATCH}/data/references/aedes_super_index/STAR_index"
BOWTIE_SUPER_REF="${PROJECT_SCRATCH}/data/references/aedes_super_index/bowtie2_index/superreference"

# Input: trimmed FASTQ files (in sample subdirectories)
INPUT_DIR="${PROJECT_SCRATCH}/results/trimmed"

# Output: alignment results
OUTPUT_DIR="${PROJECT_SCRATCH}/results/aligned"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
STAR_CONTAINER="${CONTAINERS}/star_2.7.10a.sif"
BOWTIE2_CONTAINER="${CONTAINERS}/bowtie2_2.5.5.sif"

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
echo "  Sequence Mapping to Superreference Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Superreference: ${STAR_SUPER_REF}"
echo "  Input: ${INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Containers: ${STAR_CONTAINER} & ${BOWTIE2_CONTAINER}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting mapping for ${PROJECT_NAME}" 2>/dev/null || true

####==================================####
####            CHECK INPUT           ####
####==================================####

# Check if input directories exist
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory ${INPUT_DIR} does not exist"
    tg_send "ERROR: Input directory not found" 2>/dev/null || true
    exit 1
fi

if [[ ! -d "$STAR_SUPER_REF" ]]; then
    echo "ERROR: Star superreference directory ${STAR_SUPER_REF} does not exist"
    tg_send "ERROR: Star superreference directory not found" 2>/dev/null || true
    exit 1
fi

if [[ ! -f "${STAR_SUPER_REF}/Genome" ]] || [[ ! -f "${STAR_SUPER_REF}/SA" ]]; then
    echo "ERROR: STAR index files not found in ${STAR_SUPER_REF}"
    echo "Expected files: Genome, SA, SAindex, chrLength.txt ..."
    tg_send "ERROR: STAR index files not found" 2>/dev/null || true
    exit 1
fi

if [[ ! -f "${BOWTIE_SUPER_REF}.1.bt2" ]] && [[ ! -f "${BOWTIE_SUPER_REF}.1.bt2l" ]]; then
    echo "ERROR: Bowtie2 superreference index files not found at ${BOWTIE_SUPER_REF}"
    echo "Expected files: ${BOWTIE_SUPER_REF}.1.bt2[1] or .1.bt2"
    tg_send "ERROR: Bowtie2 superreference index files not found" 2>/dev/null || true
    exit 1
fi

####==================================####
####          GET SAMPLE LIST         ####
####==================================####

cd "$INPUT_DIR" || exit 1

# Get list of sample directories (e.g., PM2486, PM2494, etc.)
samples=()
while IFS= read -r line; do
    samples+=("$line")
done < <(find . -maxdepth 1 -type d -name "PM*" | sed 's/\.\///' | sort)

# Check if samples exist
if [[ ${#samples[@]} -eq 0 ]]; then
    echo "No sample directories found in ${INPUT_DIR}"
    tg_send "No sample directories found" 2>/dev/null || true
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
echo "  Sequence STAR Mapping Array Task ${LSB_JOBINDEX}/${#samples[@]}"
echo "  Sample: ${sample}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Mapping: ${sample} (${LSB_JOBINDEX}/${#samples[@]})" 2>/dev/null || true

####================================####
####     DECOMPRESS INPUT FILES      ####
####================================####

echo "Checking for compressed input files..."

# Decompress R1 paired if compressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R1_paired.fastq.gz" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R1_paired.fastq" ]]; then
    echo "Decompressing: ${sample}_R1_paired.fastq.gz"
    gzip -d "${INPUT_DIR}/${sample}/${sample}_R1_paired.fastq.gz"
fi

# Decompress R2 paired if compressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq.gz" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq" ]]; then
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

####================================####
####     RUN MAPPING WITH STAR      ####
####================================####

# Create STAR output directory
mkdir -p "${OUTPUT_DIR}/STAR_alignment"

# Check if container exists
if [[ ! -f "$STAR_CONTAINER" ]]; then
    echo "ERROR: Container ${STAR_CONTAINER} not found"
    tg_send "ERROR: STAR container not found" 2>/dev/null || true
    exit 1
fi

# Get R1 paired file
R1_PAIRED="${INPUT_DIR}/${sample}/${sample}_R1_paired.fastq"
R2_PAIRED="${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq"
R1_UNPAIRED="${INPUT_DIR}/${sample}/${sample}_R1_unpaired.fastq"
R2_UNPAIRED="${INPUT_DIR}/${sample}/${sample}_R2_unpaired.fastq"

# Parameters
THREADS=8
RAM=64000000000  # 32GB

echo "Aligning sample with STAR: ${sample}"
tg_send "Aligning ${sample} with STAR" 2>/dev/null || true

# Check if paired files exist
if [[ -f "$R1_PAIRED" ]] && [[ -f "$R2_PAIRED" ]]; then
    echo "Found paired-end files:"
    echo "  R1: ${R1_PAIRED}"
    echo "  R2: ${R2_PAIRED}"
    
    # Paired-end alignment
    if apptainer exec \
        --bind "${INPUT_DIR}:/input:ro" \
        --bind "${OUTPUT_DIR}/STAR_alignment:/output" \
        --bind "${STAR_SUPER_REF}:/genome:ro" \
        "$STAR_CONTAINER" \
        STAR \
        --runMode alignReads \
        --genomeDir "/genome" \
        --readFilesIn "/input/${sample}/${sample}_R1_paired.fastq" "/input/${sample}/${sample}_R2_paired.fastq" \
        --outFileNamePrefix "/output/${sample}_paired_" \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outFilterMismatchNoverLmax 0.1 \
        --runThreadN "$THREADS" \
        --limitBAMsortRAM "$RAM" 2>&1; then
        echo "STAR paired-end alignment completed for ${sample}"
        tg_send "STAR paired-end: ${sample} complete" 2>/dev/null || true
    else
        echo "STAR paired-end alignment FAILED for ${sample}"
        tg_send "STAR paired-end: ${sample} FAILED" 2>/dev/null || true
    fi
fi

# Single-end alignment for unpaired reads
if [[ -f "$R1_UNPAIRED" ]]; then
    echo "Found R1 unpaired file: ${R1_UNPAIRED}"
    
    if apptainer exec \
        --bind "${INPUT_DIR}:/input:ro" \
        --bind "${OUTPUT_DIR}/STAR_alignment:/output" \
        --bind "${STAR_SUPER_REF}:/genome:ro" \
        "$STAR_CONTAINER" \
        STAR \
        --runMode alignReads \
        --genomeDir "/genome" \
        --readFilesIn "/input/${sample}/${sample}_R1_unpaired.fastq" \
        --outFileNamePrefix "/output/${sample}_R1_unpaired_" \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outFilterMismatchNoverLmax 0.1 \
        --runThreadN "$THREADS" \
        --limitBAMsortRAM "$RAM" 2>&1; then
        echo "STAR R1 unpaired alignment completed for ${sample}"
    else
        echo "STAR R1 unpaired alignment FAILED for ${sample}"
    fi
fi

if [[ -f "$R2_UNPAIRED" ]]; then
    echo "Found R2 unpaired file: ${R2_UNPAIRED}"
    
    if apptainer exec \
        --bind "${INPUT_DIR}:/input:ro" \
        --bind "${OUTPUT_DIR}/STAR_alignment:/output" \
        --bind "${STAR_SUPER_REF}:/genome:ro" \
        "$STAR_CONTAINER" \
        STAR \
        --runMode alignReads \
        --genomeDir "/genome" \
        --readFilesIn "/input/${sample}/${sample}_R2_unpaired.fastq" \
        --outFileNamePrefix "/output/${sample}_R2_unpaired_" \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outFilterMismatchNoverLmax 0.1 \
        --runThreadN "$THREADS" \
        --limitBAMsortRAM "$RAM" 2>&1; then
        echo "STAR R2 unpaired alignment completed for ${sample}"
    else
        echo "STAR R2 unpaired alignment FAILED for ${sample}"
    fi
fi

echo "Done STAR processing: ${sample}"

# Compress STAR output files
echo "Compressing STAR output files..."

# Compress STAR SAM files to BAM if samtools available, else gzip
if command -v samtools &> /dev/null; then
    for sam_file in "${OUTPUT_DIR}/STAR_alignment/${sample}"*.sam; do
        if [[ -f "$sam_file" ]]; then
            bam_file="${sam_file%.sam}.bam"
            if [[ ! -f "$bam_file" ]] || [[ "$sam_file" -nt "$bam_file" ]]; then
                echo "  Converting to BAM: $(basename $sam_file)"
                samtools view -bS "$sam_file" > "$bam_file" && rm "$sam_file"
            fi
        fi
    done
else
    # If samtools not available, gzip SAM files
    for sam_file in "${OUTPUT_DIR}/STAR_alignment/${sample}"*.sam; do
        if [[ -f "$sam_file" ]] && [[ ! -f "${sam_file}.gz" ]]; then
            echo "  Compressing: $(basename $sam_file)"
            gzip "$sam_file"
        fi
    done
fi

# Compress STAR FASTQ outputs
for fastq_file in "${OUTPUT_DIR}/STAR_alignment/${sample}"*.fastq; do
    if [[ -f "$fastq_file" ]] && [[ ! -f "${fastq_file}.gz" ]]; then
        echo "  Compressing: $(basename $fastq_file)"
        gzip "$fastq_file"
    fi
done

# List output files
echo ""
echo "STAR output files for ${sample}:"
ls -la "${OUTPUT_DIR}"/STAR_alignment/*"${sample}"* 2>/dev/null || echo "No output files found"

####===================================####
####     RUN MAPPING WITH BOWTIE2      ####
####===================================####

# Create Bowtie2 output directory
mkdir -p "${OUTPUT_DIR}/Bowtie2_alignment"

# Check if container exists
if [[ ! -f "$BOWTIE2_CONTAINER" ]]; then
    echo "ERROR: Container ${BOWTIE2_CONTAINER} not found"
    tg_send "ERROR: Bowtie2 container not found" 2>/dev/null || true
    exit 1
fi

echo "Aligning sample with Bowtie2: ${sample}"
tg_send "Aligning ${sample} with Bowtie2" 2>/dev/null || true

# Check if paired files exist
if [[ -f "$R1_PAIRED" ]] && [[ -f "$R2_PAIRED" ]]; then
    echo "Found paired-end files:"
    echo "  R1: ${R1_PAIRED}"
    echo "  R2: ${R2_PAIRED}"
    
    if apptainer exec \
        --bind "${INPUT_DIR}:/input:ro" \
        --bind "${OUTPUT_DIR}/Bowtie2_alignment:/output" \
        --bind "$(dirname "${BOWTIE_SUPER_REF}"):/bowtie_index:ro" \
        "$BOWTIE2_CONTAINER" \
        bowtie2 \
        -x "/bowtie_index/$(basename "${BOWTIE_SUPER_REF}")" \
        -1 "/input/${sample}/${sample}_R1_paired.fastq" \
        -2 "/input/${sample}/${sample}_R2_paired.fastq" \
        -S "/output/${sample}_paired.sam" \
        --un-conc "/output/${sample}_both_unmapped.fastq" \
        --un "/output/${sample}_unmapped_mixed.fastq" \
        --al "/output/${sample}_aligned_concordant.fastq" \
        --threads "$THREADS" \
        --sensitive \
        --rg-id "${sample}" \
        --rg "SM:${sample}" \
        --rg "PL:ILLUMINA" \
        2>&1; then
        echo "Bowtie2 paired-end alignment completed for ${sample}"
        tg_send "Bowtie2 paired-end: ${sample} complete" 2>/dev/null || true
    else
        echo "Bowtie2 paired-end alignment FAILED for ${sample}"
        tg_send "Bowtie2 paired-end: ${sample} FAILED" 2>/dev/null || true
    fi
fi

# Single-end alignment for unpaired reads
if [[ -f "$R1_UNPAIRED" ]]; then
    echo "Found R1 unpaired file: ${R1_UNPAIRED}"
    
    if apptainer exec \
        --bind "${INPUT_DIR}:/input:ro" \
        --bind "${OUTPUT_DIR}/Bowtie2_alignment:/output" \
        --bind "$(dirname "${BOWTIE_SUPER_REF}"):/bowtie_index:ro" \
        "$BOWTIE2_CONTAINER" \
        bowtie2 \
        -x "/bowtie_index/$(basename "${BOWTIE_SUPER_REF}")" \
        -U "/input/${sample}/${sample}_R1_unpaired.fastq" \
        -S "/output/${sample}_R1_unpaired.sam" \
        --un "/output/${sample}_R1_unpaired_unmapped.fastq" \
        --threads "$THREADS" \
        --sensitive \
        --rg-id "${sample}_R1_unpaired" \
        --rg "SM:${sample}" \
        --rg "PL:ILLUMINA" \
        2>&1; then
        echo "Bowtie2 R1 unpaired alignment completed for ${sample}"
    else
        echo "Bowtie2 R1 unpaired alignment FAILED for ${sample}"
    fi
fi

if [[ -f "$R2_UNPAIRED" ]]; then
    echo "Found R2 unpaired file: ${R2_UNPAIRED}"

    if apptainer exec \
        --bind "${INPUT_DIR}:/input:ro" \
        --bind "${OUTPUT_DIR}/Bowtie2_alignment:/output" \
        --bind "$(dirname "${BOWTIE_SUPER_REF}"):/bowtie_index:ro" \
        "$BOWTIE2_CONTAINER" \
        bowtie2 \
        -x "/bowtie_index/$(basename "${BOWTIE_SUPER_REF}")" \
        -U "/input/${sample}/${sample}_R2_unpaired.fastq" \
        -S "/output/${sample}_R2_unpaired.sam" \
        --un "/output/${sample}_R2_unpaired_unmapped.fastq" \
        --threads "$THREADS" \
        --sensitive \
        --rg-id "${sample}_R2_unpaired" \
        --rg "SM:${sample}" \
        --rg "PL:ILLUMINA" \
        2>&1; then
        echo "Bowtie2 R2 unpaired alignment completed for ${sample}"
    else
        echo "Bowtie2 R2 unpaired alignment FAILED for ${sample}"
    fi
fi

echo "Done Bowtie2 processing: ${sample}"

####================================####
####     RECOMPRESS INPUT FILES      ####
####================================####

echo "Recompressing input FASTQ files..."

# Recompress R1 paired if it was decompressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R1_paired.fastq" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R1_paired.fastq.gz" ]]; then
    echo "Recompressing: ${sample}_R1_paired.fastq"
    gzip "${INPUT_DIR}/${sample}/${sample}_R1_paired.fastq"
fi

# Recompress R2 paired if it was decompressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq.gz" ]]; then
    echo "Recompressing: ${sample}_R2_paired.fastq"
    gzip "${INPUT_DIR}/${sample}/${sample}_R2_paired.fastq"
fi

# Recompress R1 unpaired if it was decompressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R1_unpaired.fastq" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R1_unpaired.fastq.gz" ]]; then
    echo "Recompressing: ${sample}_R1_unpaired.fastq"
    gzip "${INPUT_DIR}/${sample}/${sample}_R1_unpaired.fastq"
fi

# Recompress R2 unpaired if it was decompressed
if [[ -f "${INPUT_DIR}/${sample}/${sample}_R2_unpaired.fastq" ]] && [[ ! -f "${INPUT_DIR}/${sample}/${sample}_R2_unpaired.fastq.gz" ]]; then
    echo "Recompressing: ${sample}_R2_unpaired.fastq"
    gzip "${INPUT_DIR}/${sample}/${sample}_R2_unpaired.fastq"
fi

echo ""
echo "========================================="
echo "  Sample ${sample} processing complete"
echo "  End time: $(date)"
echo "========================================="

tg_send "Completed: ${sample}" 2>/dev/null || true
