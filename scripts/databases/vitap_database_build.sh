#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to build VITAP database
# 13/09/2026
# Ver. 1.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"

# --- DATABASE DIRECTORIES ---
DB_BASE="${PERMANENT_BASE}/databases"

# Input directory -  where VMR csv is
INPUT_DIR="${DB_BASE}/vitap_db"
# VMR ICTV file
VMR_FILE="${INPUT_DIR}/VMR_MSL41_vitap.csv"

# Output directory
OUTPUT_DIR="${DB_BASE}/vitap_db/db"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
VITAP_CONTAINER="${CONTAINERS}/vitap_1.12--pyhdfd78af_0.sif"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  VITAP database construction"
echo "  VMR:    ${VMR_FILE}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date:   $(date)"
echo "========================================="

tg_send "Starting VITAP database build" 2>/dev/null || true

####===============================####
####       BUILD VITAP DATABASE    ####
####===============================####

mkdir -p "${INPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

echo "Building VITAP database"

echo "Y" | apptainer exec \
    --bind "${INPUT_DIR}":/input:ro \
    --bind "${OUTPUT_DIR}":/output \
    --env "TMPDIR=/tmp" \
    "${VITAP_CONTAINER}" \
    VITAP upd \
        --vmr "/input/$(basename "${VMR_FILE}")" \
        -o "/output/ICTV_VMR_reformat.csv" \
        -d VMR-MSL

EXIT_CODE=$?

if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "VITAP database built successfully in ${OUTPUT_DIR}"
    [[ -n "${TG_SCRIPT:-}" ]] && tg_send "VITAP DB build completed" 2>/dev/null || true
else
    echo "VITAP database build FAILED (exit ${EXIT_CODE})"
    [[ -n "${TG_SCRIPT:-}" ]] && tg_send "VITAP DB build FAILED" 2>/dev/null || true
    exit ${EXIT_CODE}
fi
