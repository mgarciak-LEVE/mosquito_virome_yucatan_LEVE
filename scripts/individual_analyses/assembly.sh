  GNU nano 6.2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      assembly.sh                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to assemble unmapped reads from STAR and Bowtie2
# 03/08/2026
# Ver. 1.4.0 (farm-ready)

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
OUTPUT_DIR="${PROJECT_SCRATCH}/results/assembly"
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
echo "  Unmapped Reads Assembly Module"
echo "  Project: ${PROJECT_NAME}"
echo "  STAR input: ${STAR_INPUT_DIR}"
echo "  Bowtie input: ${BOWTIE_INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting unmapped reads assembly" 2>/dev/null || true

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

####==================================####
####     GET SAMPLES FROM BOWTIE      ####
####==================================####

cd "$BOWTIE_INPUT_DIR" || exit 1

# Bowtie paired samples (interleaved file with both reads)
bowtie_paired_samples=()
while IFS= read -r line; do
    sample_name=$(basename "$line" "_unmapped_mixed.fastq")
    bowtie_paired_samples+=("$sample_name")
done < <(ls -1 *_unmapped_mixed.fastq 2>/dev/null | sort)

# Bowtie R1 unpaired samples
bowtie_r1_unpaired_samples=()
while IFS= read -r line; do
    sample_name=$(basename "$line" "_R1_unpaired_unmapped.fastq")
    bowtie_r1_unpaired_samples+=("$sample_name")
done < <(ls -1 *_R1_unpaired_unmapped.fastq 2>/dev/null | sort)

# Bowtie R2 unpaired samples
bowtie_r2_unpaired_samples=()
while IFS= read -r line; do
    sample_name=$(basename "$line" "_R2_unpaired_unmapped.fastq")
    bowtie_r2_unpaired_samples+=("$sample_name")
done < <(ls -1 *_R2_unpaired_unmapped.fastq 2>/dev/null | sort)

# Count Bowtie samples
bowtie_total=$(( ${#bowtie_paired_samples[@]} + ${#bowtie_r1_unpaired_samples[@]} + ${#bowtie_r2_unpaired_samples[@]} ))

####===================================================####
####          COMBINE ALL SAMPLES INTO ONE ARRAY       ####
####===================================================####

all_samples=()
all_types=()
all_sources=()

# Add STAR paired samples
for sample in "${star_paired_samples[@]}"; do
    all_samples+=("$sample")
    all_types+=("paired")
    all_sources+=("STAR")
done

# Add STAR R1 unpaired samples
for sample in "${star_r1_unpaired_samples[@]}"; do
    all_samples+=("$sample")
    all_types+=("R1_unpaired")
    all_sources+=("STAR")
done

# Add STAR R2 unpaired samples
for sample in "${star_r2_unpaired_samples[@]}"; do
    all_samples+=("$sample")
    all_types+=("R2_unpaired")
    all_sources+=("STAR")
done

# Add Bowtie paired samples
for sample in "${bowtie_paired_samples[@]}"; do
    all_samples+=("$sample")
    all_types+=("paired_interleaved")
    all_sources+=("Bowtie")
done

# Add Bowtie R1 unpaired samples
for sample in "${bowtie_r1_unpaired_samples[@]}"; do
    all_samples+=("$sample")
    all_types+=("R1_unpaired")
    all_sources+=("Bowtie")
done

# Add Bowtie R2 unpaired samples
for sample in "${bowtie_r2_unpaired_samples[@]}"; do
    all_samples+=("$sample")
    all_types+=("R2_unpaired")
    all_sources+=("Bowtie")
done

total_samples=${#all_samples[@]}

if [[ $total_samples -eq 0 ]]; then
    echo "No unmapped samples found"
    tg_send "No unmapped samples found" 2>/dev/null || true
    exit 1
fi

echo ""
echo "Summary:"
echo "  - STAR samples: ${star_total}"
echo "    - Paired: ${#star_paired_samples[@]}"
echo "    - R1 unpaired: ${#star_r1_unpaired_samples[@]}"
echo "    - R2 unpaired: ${#star_r2_unpaired_samples[@]}"
echo "  - Bowtie samples: ${bowtie_total}"
echo "    - Paired (interleaved): ${#bowtie_paired_samples[@]}"
echo "    - R1 unpaired: ${#bowtie_r1_unpaired_samples[@]}"
echo "    - R2 unpaired: ${#bowtie_r2_unpaired_samples[@]}"
echo "  - Total: ${total_samples}"

####===================================================####
####          GET ARRAY TASK SAMPLE                    ####
####===================================================####

SAMPLE_INDEX=$((LSB_JOBINDEX - 1))

if [[ $SAMPLE_INDEX -ge ${#all_samples[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#all_samples[@]})"
    exit 1
fi

sample="${all_samples[$SAMPLE_INDEX]}"
sample_type="${all_types[$SAMPLE_INDEX]}"
sample_source="${all_sources[$SAMPLE_INDEX]}"

echo "========================================="
echo "  Processing sample: ${sample}"
echo "  Type: ${sample_type}"
echo "  Source: ${sample_source}"
echo "  Array task: ${LSB_JOBINDEX}/${#all_samples[@]}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${sample} (${sample_type}) from ${sample_source} (${LSB_JOBINDEX}/${#all_samples[@]})" 2>/dev/null || true

####=====================================####
####        DETERMINE INPUT FILES        ####
####=====================================####

# Create FASTQ directory
mkdir -p "${FASTQ_DIR}/${sample}"

case "${sample_source}_${sample_type}" in
    "STAR_paired")
        # STAR unmapped reads
        R1_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${sample}_paired_Unmapped.out.mate1"
        R2_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${sample}_paired_Unmapped.out.mate2"

        # Copy with .fastq extension
        R1_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}_R1.fastq"
        R2_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}_R2.fastq"

        if [[ ! -s "$R1_UNMAPPED_ORIG" ]] || [[ ! -s "$R2_UNMAPPED_ORIG" ]]; then
            echo " WARNING: Unmapped files for ${sample} are empty or missing"
            echo "  Skipping assembly for ${sample}"
            exit 0
        fi

        echo "  Copying unmapped reads to FASTQ directory..."

        # Copy files to FASTQ directory with .fastq extension
        cp -p "$R1_UNMAPPED_ORIG" "$R1_UNMAPPED"
        cp -p "$R2_UNMAPPED_ORIG" "$R2_UNMAPPED"

        # Compress original files to save space
        echo "  Compressing original unmapped files..."
        gzip -f "$R1_UNMAPPED_ORIG"
        gzip -f "$R2_UNMAPPED_ORIG"

        INPUT_DIR_FOR_BIND="${FASTQ_DIR}/${sample}"

        R1_READS=$(($(wc -l < "$R1_UNMAPPED") / 4))
        R2_READS=$(($(wc -l < "$R2_UNMAPPED") / 4))
        echo "  Using STAR unmapped reads (copied to FASTQ directory):"
        echo "  R1: $(basename "$R1_UNMAPPED") ($R1_READS reads)"
        echo "  R2: $(basename "$R2_UNMAPPED") ($R2_READS reads)"
        echo "  Original files compressed: $(basename "$R1_UNMAPPED_ORIG").gz and $(basename "$R2_UNMAPPED_ORIG").gz"
        ;;

    "STAR_R1_unpaired")
        R1_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${sample}_R1_unpaired_Unmapped.out.mate1"
        R1_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}_R1.fastq"

        if [[ ! -s "$R1_UNMAPPED_ORIG" ]]; then
            echo " WARNING: Unmapped file for ${sample} is empty or missing"
            echo "  Skipping assembly for ${sample}"
            exit 0
        fi

        echo "  Copying R1 unmapped reads to FASTQ directory..."

        # Copy file to FASTQ directory with .fastq extension
        cp -p "$R1_UNMAPPED_ORIG" "$R1_UNMAPPED"

        # Compress original file to save space
        echo "  Compressing original R1 unmapped file..."
        gzip -f "$R1_UNMAPPED_ORIG"

        INPUT_DIR_FOR_BIND="${FASTQ_DIR}/${sample}"

        R1_READS=$(($(wc -l < "$R1_UNMAPPED") / 4))
        echo "  Using STAR R1 unmapped reads (copied to FASTQ directory):"
        echo "  R1: $(basename "$R1_UNMAPPED") ($R1_READS reads)"
        echo "  Original file compressed: $(basename "$R1_UNMAPPED_ORIG").gz"
        ;;

    "STAR_R2_unpaired")
        R2_UNMAPPED_ORIG="${STAR_INPUT_DIR}/${sample}_R2_unpaired_Unmapped.out.mate1"
        R2_UNMAPPED="${FASTQ_DIR}/${sample}/${sample}_R2.fastq"

        if [[ ! -s "$R2_UNMAPPED_ORIG" ]]; then
            echo " WARNING: Unmapped file for ${sample} is empty or missing"
            echo "  Skipping assembly for ${sample}"
            exit 0
        fi

        echo "  Copying R2 unmapped reads to FASTQ directory..."

        # Copy file to FASTQ directory with .fastq extension
        cp -p "$R2_UNMAPPED_ORIG" "$R2_UNMAPPED"

        # Compress original file to save space
        echo "  Compressing original R2 unmapped file..."
        gzip -f "$R2_UNMAPPED_ORIG"

        INPUT_DIR_FOR_BIND="${FASTQ_DIR}/${sample}"

        R2_READS=$(($(wc -l < "$R2_UNMAPPED") / 4))
        echo "  Using STAR R2 unmapped reads (copied to FASTQ directory):"
        echo "  R2: $(basename "$R2_UNMAPPED") ($R2_READS reads)"
        echo "  Original file compressed: $(basename "$R2_UNMAPPED_ORIG").gz"
        ;;

    "Bowtie_paired_interleaved")
        R1_UNMAPPED="${BOWTIE_INPUT_DIR}/${sample}_both_unmapped.fastq"
        if [[ ! -f "$R1_UNMAPPED" ]]; then
            R1_UNMAPPED="${BOWTIE_INPUT_DIR}/${sample}_unmapped_mixed.fastq"
        fi
        R2_UNMAPPED=""
        INPUT_DIR_FOR_BIND="${BOWTIE_INPUT_DIR}"
        echo "Bowtie paired-end unmapped reads: $(basename "$R1_UNMAPPED")"
        ;;

    "Bowtie_R1_unpaired")
        R1_UNMAPPED="${BOWTIE_INPUT_DIR}/${sample}_R1_unpaired_unmapped.fastq"
        INPUT_DIR_FOR_BIND="${BOWTIE_INPUT_DIR}"
        echo "Bowtie R1 unpaired: $(basename "$R1_UNMAPPED")"
        ;;

    "Bowtie_R2_unpaired")
        R2_UNMAPPED="${BOWTIE_INPUT_DIR}/${sample}_R2_unpaired_unmapped.fastq"
        INPUT_DIR_FOR_BIND="${BOWTIE_INPUT_DIR}"
        echo "Bowtie R2 unpaired: $(basename "$R2_UNMAPPED")"
        ;;

    *)
        echo "ERROR: Unknown combination: ${sample_source}_${sample_type}"
        exit 1
        ;;
esac

# Check if input files exist
if [[ "$sample_type" == "paired" ]] || [[ "$sample_type" == "paired_interleaved" ]]; then
    if [[ ! -f "$R1_UNMAPPED" ]]; then
        echo "ERROR: Input file not found: ${R1_UNMAPPED}"
        exit 1
    fi
    if [[ "$sample_type" == "paired" ]] && [[ ! -f "$R2_UNMAPPED" ]]; then
        echo "ERROR: Input file not found: ${R2_UNMAPPED}"
        exit 1
    fi
else
    if [[ ! -f "$R1_UNMAPPED" ]] && [[ ! -f "$R2_UNMAPPED" ]]; then
        echo "ERROR: Input file not found"
        exit 1
    fi
fi

####================================####
####          RUN ASSEMBLY          ####
####================================####

# Parameters
THREADS=32
MEMORY=128  # GB for MEGAhit

# Clean output directories before running
rm -rf "${OUTPUT_DIR}/rnaSPAdes/${sample}"/*
rm -rf "${OUTPUT_DIR}/metaSPAdes/${sample}"/*
rm -rf "${OUTPUT_DIR}/MEGAhit/${sample}"/*

# Create output directories
mkdir -p "${OUTPUT_DIR}/rnaSPAdes/${sample}"
mkdir -p "${OUTPUT_DIR}/metaSPAdes/${sample}"
mkdir -p "${OUTPUT_DIR}/MEGAhit/${sample}"

echo "Starting assembly process for ${sample} (${sample_type}) from ${sample_source}..."
tg_send "Starting assembly for ${sample}" 2>/dev/null || true




####================================####
####        rnaSPAdes Assembly      ####
####================================####

echo "Starting rnaSPAdes assembly for sample ${sample}..."

if [[ "$sample_type" == "paired" ]]; then
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        rnaspades.py \
        -1 "/input/$(basename "$R1_UNMAPPED")" \
        -2 "/input/$(basename "$R2_UNMAPPED")" \
        -o "/output/rnaSPAdes/${sample}" \
        -t "$THREADS"
elif [[ "$sample_type" == "paired_interleaved" ]]; then
    # Bowtie paired-end (interleaved file)
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        rnaspades.py \
        --interleaved "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/rnaSPAdes/${sample}" \
        -t "$THREADS"
else
    # Single-end
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        rnaspades.py \
        -s "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/rnaSPAdes/${sample}" \
        -t "$THREADS"
fi

if [[ $? -eq 0 ]]; then
    echo "rnaSPAdes assembly completed for ${sample}"
else
    echo "rnaSPAdes assembly FAILED for ${sample}"
fi

####================================####
####        metaSPAdes Assembly      ####
####================================####

echo "Starting metaSPAdes assembly for sample ${sample}..."

if [[ "$sample_type" == "paired" ]]; then
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        metaspades.py \
        -1 "/input/$(basename "$R1_UNMAPPED")" \
        -2 "/input/$(basename "$R2_UNMAPPED")" \
        -o "/output/metaSPAdes/${sample}" \
        -t "$THREADS"
elif [[ "$sample_type" == "paired_interleaved" ]]; then
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        metaspades.py \
        --interleaved "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/metaSPAdes/${sample}" \
        -t "$THREADS"
else
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        metaspades.py \
        -s "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/metaSPAdes/${sample}" \
        -t "$THREADS"
fi

if [[ $? -eq 0 ]]; then
    echo "metaSPAdes assembly completed for ${sample}"
else
    echo "metaSPAdes assembly FAILED for ${sample}"
fi

####================================####
####     metaviralSPAdes Assembly   ####
####================================####

echo "Starting metaviralSPAdes assembly for sample ${sample}..."

if [[ "$sample_type" == "paired" ]]; then
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        metaviralspades.py \
        -1 "/input/$(basename "$R1_UNMAPPED")" \
        -2 "/input/$(basename "$R2_UNMAPPED")" \
        -o "/output/metaviralSPAdes/${sample}" \
        -t "$THREADS"
elif [[ "$sample_type" == "paired_interleaved" ]]; then
    # Bowtie paired-end (interleaved file)
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        metaviralspades.py \
        --interleaved "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/metaviralSPAdes/${sample}" \
        -t "$THREADS"
else
    # Single-end
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$SPADES_CONTAINER" \
        metaviralspades.py \
        -s "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/metaviralSPAdes/${sample}" \
        -t "$THREADS"
fi

if [[ $? -eq 0 ]]; then
    echo "metaviralSPAdes assembly completed for ${sample}"
else
    echo "metaviralSPAdes assembly FAILED for ${sample}"
fi


####================================####
####        MEGAhit Assembly        ####
####================================####

echo "Starting MEGAhit assembly for sample ${sample}..."

rm -rf "${OUTPUT_DIR}/MEGAhit/${sample}" 2>/dev/null

if [[ "$sample_type" == "paired" ]]; then
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$MEGAHIT_CONTAINER" \
        megahit \
        -1 "/input/$(basename "$R1_UNMAPPED")" \
        -2 "/input/$(basename "$R2_UNMAPPED")" \
        -o "/output/MEGAhit/${sample}" \
        -t "$THREADS" \
        -m "${MEMORY}"
elif [[ "$sample_type" == "paired_interleaved" ]]; then
    # MEGAhit doesn't support interleaved directly, convert to separate files first
    # For now, use single-end mode with the interleaved file
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$MEGAHIT_CONTAINER" \
        megahit \
        -r "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/MEGAhit/${sample}" \
        -t "$THREADS" \
        -m "${MEMORY}"
else
    apptainer exec \
        --bind "${INPUT_DIR_FOR_BIND}:/input:ro" \
        --bind "${OUTPUT_DIR}:/output" \
        "$MEGAHIT_CONTAINER" \
        megahit \
        -r "/input/$(basename "$R1_UNMAPPED")" \
        -o "/output/MEGAhit/${sample}" \
        -t "$THREADS" \
        -m "${MEMORY}"
fi

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
echo "  Assembly completed for ${sample}"
echo "  Output directories:"
echo "    - rnaSPAdes: ${OUTPUT_DIR}/rnaSPAdes/${sample}"
echo "    - metaSPAdes: ${OUTPUT_DIR}/metaSPAdes/${sample}"
echo "    - MEGAhit: ${OUTPUT_DIR}/MEGAhit/${sample}"
echo "    - metaviralSPAdes: ${OUTPUT_DIR}/metaviralSPAdes/${sample}"
echo "========================================="

tg_send "Completed: ${sample}" 2>/dev/null || true
