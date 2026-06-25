#!/bin/bash
# Sequence Quality Module
# Author: Jorge Alberto Castro Rodríguez
# Ver. 2.0.0
# 06/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Load Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

# Working directories
WORKDIR="${1:-${PROJECT_ROOT}/data/raw/total_RNA/cat_files}"
OUTPUT_BASE="${2:-${PROJECT_ROOT}/results}"

# Container configuration
CONTAINER="${PROJECT_ROOT}/containers/pipeline_calidad.sif"

####==================================####
####         FASTQC FUNCTION          ####
####==================================####

run_fastqc() {
    local input_dir=$1
    local output_dir=$2

    echo "Running FastQC"
    tg_send "Running FastQC" 
    echo "Input: ${input_dir}"
    tg_send "Input: ${input_dir}"
    echo "Output: ${output_dir}"
    tg_send "Output: ${output_dir}"

    mkdir -p "${output_dir}"

    # Input directory validation
    if [ ! -d "$input_dir" ]; then
        echo "Error: directory ${input_dir} does not exist"
        return 1
    fi

    # File counting
    local total_files=$(find "$input_dir" \( -name "*.fastq" -o -name "*.fastq.gz" \) -type f | wc -l)
    local current=0
    
    while IFS= read -r file; do
        ((current++))
        local file_basename=$(basename "$file")
        echo "[${current}/${total_files}] Processing: ${file_basename}"
        tg_send "[${current}/${total_files}] Processing: ${file_basename}"

        # Container with FastQC
        apptainer exec \
            --bind "${input_dir}:/input:ro" \
            --bind "${output_dir}:/output" \
            "$CONTAINER" \
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

    echo "  Running MultiQC"
    tg_send "Running MultiQC" 
    echo "Input: ${input_dir}"
    tg_send "Input: ${input_dir}"
    echo "Output: ${output_dir}"
    tg_send "Output: ${output_dir}"

    
    mkdir -p "${output_dir}"

    # Run MultiQC in container
    apptainer exec \
        --bind "${input_dir}:/input:ro" \
        --bind "${output_dir}:/output" \
        "$CONTAINER" \
        multiqc "/input" -o "/output"
    
    echo "MultiQC completed successfully"
}

####=================================####
####        FUNCTION EXECUTION       ####
####=================================####

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Script is being executed directly
    echo "Sequence quality inspection module"
    
    FASTQC_OUTPUT="${OUTPUT_BASE}/untrimmed_qc/fastqc"
    MULTIQC_OUTPUT="${OUTPUT_BASE}/untrimmed_qc/multiqc"

    # Run the pipeline steps
    run_fastqc "$WORKDIR" "$FASTQC_OUTPUT"
    run_multiqc "$FASTQC_OUTPUT" "$MULTIQC_OUTPUT"

fi