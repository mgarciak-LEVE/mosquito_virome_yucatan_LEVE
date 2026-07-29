#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Build Bowtie2 index from superreference
# 28/07/2026
# Ver. 1.0.1 (farm-ready)

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

SUPERREF="${PROJECT_SCRATCH}/data/references/aedes_super_index/superreference.fna"
BOWTIE2_INDEX="${PROJECT_SCRATCH}/data/references/aedes_super_index/bowtie2_index/superreference"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
BOWTIE2_CONTAINER="${CONTAINERS}/bowtie2_2.5.5.sif"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Bowtie2 Superreference Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Superreference: ${SUPERREF}"
echo "  Output: ${BOWTIE2_INDEX}"
echo "  Container: ${BOWTIE2_CONTAINER}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting bowtie superreference creation for ${PROJECT_NAME}" 2>/dev/null || true

mkdir -p "$(dirname "$BOWTIE2_INDEX")"

# Build index
apptainer exec \
    --bind "$(dirname "$SUPERREF"):/input:ro" \
    --bind "$(dirname "$BOWTIE2_INDEX"):/output" \
    "$BOWTIE2_CONTAINER" \
    bowtie2-build \
    "/input/$(basename "$SUPERREF")" \
    "/output/$(basename "$BOWTIE2_INDEX")"

echo "Bowtie2 index built at: ${BOWTIE2_INDEX}"