#!/bin/bash
# Author: Jorge Alberto Castro Rodríguez
# GENOMAD  database construction
# Ver. 1.1.0

####==================================####
####           CONFIGURATION          ####
####==================================####
set -euo pipefail

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- Storage ---
PERMANENT_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46"
DB_BASE="${PERMANENT_BASE}/databases"
GENOMAD_DB="${DB_BASE}/genomad_db"
SCRATCH_BASE="${PERMANENT_BASE}/projects"

# --- Container ---
CONTAINERS="${PERMANENT_BASE}/containers"
GENOMAD_CONTAINER="${CONTAINERS}/genomad_1.12.0--pyhdfd78af_0.sif"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  GENOMAD database construction"
echo "  Output: ${GENOMAD_DB}"
echo "  Date:   $(date)"
echo "========================================="

tg_send "Starting GENOMAD database build" 2>/dev/null || true

####===============================####
####       BUILD VITAP DATABASE    ####
####===============================####

mkdir -p "${GENOMAD_DB}"

echo "Building GENOMAD database"

apptainer exec \
    --bind "${GENOMAD_DB}":/output \
    "${GENOMAD_CONTAINER}" \
    genomad download-database \
        "/output"

EXIT_CODE=$?

if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "GENOMAD database built successfully in ${GENOMAD_DB}"
    [[ -n "${TG_SCRIPT:-}" ]] && tg_send "GENOMAD DB build completed" 2>/dev/null || true
else
    echo "GENOMAD database build FAILED (exit ${EXIT_CODE})"
    [[ -n "${TG_SCRIPT:-}" ]] && tg_send "GENOMAD DB build FAILED" 2>/dev/null || true
    exit ${EXIT_CODE}
fi
