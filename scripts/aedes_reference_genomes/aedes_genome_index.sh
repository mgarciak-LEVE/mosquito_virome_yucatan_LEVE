#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to make a super reference based on nine Aedes species: 
# Ae. albopictus, Ae. aegypti, Ae. mascarensis, Ae. japonicus, Ae. notoscriptus, Ae. (Ochlerotatus) sierrensis, Ae. taeniorhynchus, Ae. koreicus, and Ae. polynesiensis.
# For details on the genomes used see .csv file in the folder
# 23/07/2026
# Version 2.0.0 (farm-ready)

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
# Working directory on Lustre
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
# Input directory - where mosquito superreference is
SUPER_REF="${PROJECT_SCRATCH}/data/references/aedes_super_index/superreference.fna"
# Output directory for results
INDEX_DIR="${PROJECT_SCRATCH}/data/references/mosquito_genomes/aedes_super_index"

mkdir -p "$(dirname "${SUPER_REF}")" "${INDEX_DIR}"

# Container configuration
CONTAINERS="${HOME}/git_repos/${PROJECT_NAME}/containers"
STAR_CONTAINER="${CONTAINERS}/star_2.7.10a.sif"

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
echo "  Mosquito Indexing Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Input: ${SUPER_REF}"
echo "  Output: ${INDEX_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting genome indexing for ${PROJECT_NAME}" 2>/dev/null || true

####==================================####
####          CHECK INPUT FILE        ####
####==================================####

if [[ ! -f "$SUPER_REF" ]]; then
    echo "ERROR: Superreference file not found: ${SUPER_REF}"
    tg_send "ERROR: Superreference file not found" 2>/dev/null || true
    exit 1
fi

if [[ ! -s "$SUPER_REF" ]]; then
    echo "ERROR: Superreference file is empty: ${SUPER_REF}"
    tg_send "ERROR: Superreference file is empty" 2>/dev/null || true
    exit 1
fi

echo "Superreference file found and not empty"

####==================================####
####          GENOME INDEXING         ####
####==================================####

GENOME_CHR_BIN_BITS=14
SA_INDEX_BASES=13
THREADS=8
STAR_RAM=64000000000  # 64GB in bytes

# Remove any previous incomplete index
rm -rf "$INDEX_DIR"/*

echo "Building STAR index..."
echo "Parameters:"
echo "  • Threads: $THREADS"
echo "  • RAM limit: $STAR_RAM bytes (64GB)"
echo "  • GenomeChrBinNbits: $GENOME_CHR_BIN_BITS"
echo "  • SAindexNbases: $SA_INDEX_BASES"

tg_send "Building STAR index" 2>/dev/null || true

if apptainer exec \
    --bind "${SUPER_REF}:/input:ro" \
    --bind "${INDEX_DIR}:/output" \
    "$STAR_CONTAINER" \
    STAR \
    --runMode genomeGenerate \
    --genomeDir "/output" \
    --genomeFastaFiles "/input" \
    --runThreadN "$THREADS" \
    --genomeChrBinNbits "$GENOME_CHR_BIN_BITS" \
    --limitGenomeGenerateRAM "$STAR_RAM" \
    --genomeSAindexNbases "$SA_INDEX_BASES" 2>&1; then
    echo "STAR index completed"
    tg_send "STAR index completed" 2>/dev/null || true
else
    echo "STAR index FAILED"
    tg_send "STAR index FAILED" 2>/dev/null || true
    exit 1
fi

echo ""
echo "========================================="
echo "  STAR Index Complete"
echo "  Index directory: ${INDEX_DIR}"
echo "  Index files:"
ls -la "$INDEX_DIR/"
echo "========================================="

tg_send "STAR indexing complete for ${PROJECT_NAME}" 2>/dev/null || true
