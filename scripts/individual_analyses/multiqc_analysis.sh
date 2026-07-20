#!/bin/bash
# MultiQC script
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 20/07/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

# Directory where scripts are
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/team222/projects"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
# Permanent results directory
PERMANENT_RESULTS="${PERMANENT_BASE}/${PROJECT_NAME}/results/untrimmed_qc"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
INPUT_MULTIQC="${PROJECT_SCRATCH}/results/untrimmed_qc/fastqc"
OUTPUT_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_MULTIQC="${OUTPUT_BASE}/untrimmed_qc/multiqc"

# Container configuration
CONTAINERS="${HOME}/git_repos/${PROJECT_NAME}/containers"
MULTIQC_CONTAINER="${CONTAINERS}/multiqc_1.35.sif"

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
echo "  MultiQC Analysis"
echo "  Project: ${PROJECT_NAME}"
echo "  Input (FastQC reports): ${INPUT_MULTIQC}"
echo "  Output: ${OUTPUT_MULTIQC}"
echo "  Container: ${MULTIQC_CONTAINER}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting MultiQC for ${PROJECT_NAME}" 2>/dev/null || true

####==================================####
####            CHECK INPUT           ####
####==================================####

# Check if input directory exists
if [[ ! -d "$INPUT_MULTIQC" ]]; then
    echo " ERROR: Input directory ${INPUT_MULTIQC} does not exist"
    tg_send " ERROR: MultiQC input directory not found" 2>/dev/null || true
    exit 1
fi

# Check if there are FastQC reports
fastqc_reports=$(find "$INPUT_MULTIQC" -maxdepth 1 -name "*fastqc.zip" 2>/dev/null | wc -l)

if [[ $fastqc_reports -eq 0 ]]; then
    echo " ERROR: No FastQC reports found in ${INPUT_MULTIQC}"
    tg_send " ERROR: No FastQC reports found" 2>/dev/null || true
    exit 1
fi

echo "Found ${fastqc_reports} FastQC reports"
tg_send "Found ${fastqc_reports} FastQC reports" 2>/dev/null || true

# Check if container exists
if [[ ! -f "$MULTIQC_CONTAINER" ]]; then
    echo " ERROR: Container ${MULTIQC_CONTAINER} not found"
    tg_send " ERROR: MultiQC container not found" 2>/dev/null || true
    exit 1
fi


####==================================####
####          RUN MULTIQC             ####
####==================================####

# Create output directory
mkdir -p "${OUTPUT_MULTIQC}"

echo "Running MultiQC..."
tg_send "Running MultiQC" 2>/dev/null || true

if apptainer exec \
    --bind "${INPUT_MULTIQC}:/input:ro" \
    --bind "${OUTPUT_MULTIQC}:/output" \
    "$MULTIQC_CONTAINER" \
    multiqc "/input" -o "/output" 2>&1; then
    echo "MultiQC completed"
    tg_send "MultiQC completed" 2>/dev/null || true
else
    echo "MultiQC FAILED"
    tg_send "MultiQC FAILED" 2>/dev/null || true
    exit 1
fi

####==================================####
####       COPY RESULTS TO NFS        ####
####==================================####

echo ""
echo "========================================="
echo "  Copying results to permanent storage"
echo "========================================="

mkdir -p "${PERMANENT_RESULTS}"
mkdir -p "${PERMANENT_RESULTS}/multiqc"

if [[ -d "$OUTPUT_MULTIQC" ]]; then
    echo "Copying MultiQC results..."
    mkdir -p "${PERMANENT_RESULTS}/multiqc"
    cp -r "${OUTPUT_MULTIQC}"/* "${PERMANENT_RESULTS}/multiqc/" 2>/dev/null || true
    echo "MultiQC results copied to: ${PERMANENT_RESULTS}/multiqc/"
    tg_send "Results copied to NFS: ${PERMANENT_RESULTS}/multiqc/" 2>/dev/null || true
fi

echo ""
echo "========================================="
echo "  MultiQC Complete"
echo "  Input: ${INPUT_MULTIQC}"
echo "  Output: ${OUTPUT_MULTIQC}"
echo "  Permanent storage: ${PERMANENT_RESULTS}/multiqc/"
echo "========================================="

tg_send "MultiQC complete for ${PROJECT_NAME}" 2>/dev/null || true
