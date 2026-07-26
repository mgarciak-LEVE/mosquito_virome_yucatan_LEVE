#!/bin/bash
# Sequence Quality Module
# Author: Jorge Alberto Castro Rodríguez
# Ver. 2.1.0 (farm-ready)
# 14/07/2026 

####==================================####
####          CONFIGURATION           ####
####==================================####

# Directory where scripts are
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
# Permanent storage
PERMANENT_BASE="/nfs/users/nfs_j/jr46"

# Scratch storage
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
# Working directory on Lustre
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directory - where FASTQ files are 
INPUT_FASTQC="${PROJECT_SCRATCH}/data/raw/total_RNA/cat_files"

# OUTPUT DIRECTORIES
OUTPUT_BASE="${PROJECT_SCRATCH}/results" # Output base directory
OUTPUT_FASTQC="${OUTPUT_BASE}/untrimmed_qc/fastqc" # FastQC output directory
OUTPUT_MULTIQC="${OUTPUT_BASE}/untrimmed_qc/multiqc" # MultiQC output directory

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
FASTQC_CONTAINER="${CONTAINERS}/fastqc_0.11.2.sif"
MULTIQC_CONTAINER="${CONTAINERS}/multiqc_1.35.sif"

# Permanent results directory
PERMANENT_RESULTS="${PERMANENT_BASE}/${PROJECT_NAME}/results/untrimmed_qc"

# Scripts directory (on NFS)
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot (source from NFS)
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Sequence Quality Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Input directory: ${INPUT_FASTQC}"
echo "  FastQC output: ${OUTPUT_FASTQC}"
echo "  MultiQC output: ${OUTPUT_MULTIQC}"
echo "  FastQC container: ${FASTQC_CONTAINER}"
echo "  MultiQC container: ${MULTIQC_CONTAINER}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting Sequence Quality module for ${PROJECT_NAME}" 2>/dev/null || true

####==================================####
####         FASTQC FUNCTION          ####
####==================================####

run_fastqc() {
    local input_dir=$1
    local output_dir=$2

    echo "========================================="
    echo "  Running FastQC"
    echo "  Input: ${input_dir}"
    echo "  Output: ${output_dir}"
    echo "========================================="
    
    tg_send "Running FastQC on ${input_dir}" 2>/dev/null || true

    mkdir -p "${output_dir}"

    # Input directory validation
    if [ ! -d "$input_dir" ]; then
        echo "Error: directory ${input_dir} does not exist"
        return 1
    fi

    # Check if container exists
    if [[ ! -f "$FASTQC_CONTAINER" ]]; then
        echo "ERROR: Container ${FASTQC_CONTAINER} not found"
        tg_send "ERROR: Container ${FASTQC_CONTAINER} not found" 2>/dev/null || true
        return 1
    fi

    # File counting
    local total_files=$(find "$input_dir" \( -name "*.fastq" -o -name "*.fastq.gz" \) -type f | wc -l)
    local current=0

    if [[ $total_files -eq 0 ]]; then
        echo "No FASTQ files found in ${input_dir}"
        tg_send "No FASTQ files found" 2>/dev/null || true
        return 1
    fi
    
    echo "Found ${total_files} FASTQ files to process"
    tg_send "Found ${total_files} FASTQ files" 2>/dev/null || true
    
    while IFS= read -r file; do
        ((current++))
        local file_basename=$(basename "$file")
        echo "[${current}/${total_files}] Processing: ${file_basename}"
        tg_send "[${current}/${total_files}] Processing: ${file_basename}"

        # Container with FastQC
        apptainer exec \
            --bind "${input_dir}:/input:ro" \
            --bind "${output_dir}:/output" \
            --bind /usr/share/fonts:/usr/share/fonts:ro \
            --bind /usr/share/fontconfig:/usr/share/fontconfig:ro \
            --bind /etc/fonts:/etc/fonts:ro \
            "$FASTQC_CONTAINER" \
            fastqc "/input/${file_basename}" -o "/output"
    done < <(find "$input_dir" \( -name "*.fastq" -o -name "*.fastq.gz" \) -type f)
    
    echo "FastQC completed successfully"
    tg_send "FastQC completed: ${total_files} files processed"
}

####==================================####
####         MULTIQC FUNCTION        ####
####==================================####

run_multiqc() {
    local input_dir=$1
    local output_dir=$2

    echo "========================================="
    echo "  Running MultiQC"
    echo "  Input: ${input_dir}"
    echo "  Output: ${output_dir}"
    echo "========================================="

    tg_send "Running MultiQC" 2>/dev/null || true
    
    mkdir -p "${output_dir}"

    # Check if input directory exists and has FastQC reports
    if [[ ! -d "$input_dir" ]]; then
        echo "ERROR: Input directory ${input_dir} does not exist"
        tg_send "ERROR: MultiQC input directory not found" 2>/dev/null || true
        return 1
    fi
    
    # Check if there are any FastQC reports
    local fastqc_reports=$(find "$input_dir" -maxdepth 1 -name "*fastqc.zip" | wc -l)
    
    if [[ $fastqc_reports -eq 0 ]]; then
        echo "No FastQC reports found in ${input_dir}"
        tg_send "No FastQC reports found" 2>/dev/null || true
        return 1
    fi
    
    echo "Found ${fastqc_reports} FastQC reports"

    # Run MultiQC in container
    apptainer exec \
        --bind "${input_dir}:/input:ro" \
        --bind "${output_dir}:/output" \
        "$MULTIQC_CONTAINER" \
        multiqc "/input" -o "/output"
    
    echo "MultiQC completed successfully"
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
    run_multiqc "$OUTPUT_FASTQC" "$OUTPUT_MULTIQC"

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

