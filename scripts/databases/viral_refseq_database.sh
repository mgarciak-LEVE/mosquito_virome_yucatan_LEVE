#!/bin/bash
# Author: Jorge Alberto Castro Rodríguez
# Viral refseq database construction
# Ver. 1.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####
set -euo pipefail

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- Storage ---
PERMANENT_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46"
DB_BASE="${PERMANENT_BASE}/databases"
VIRALREFSEQ_DB="${DB_BASE}/blast_db/viral_refseq"
SCRATCH_BASE="${PERMANENT_BASE}/projects"
DB_NAME="refseq_viruses_rep_genomes"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/databases"
# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load blast/2.15.0--pl5321h6f7f691_1 || echo "BLAST module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Refseq BLAST database construction"
echo "  Output: ${VIRALREFSEQ_DB}"
echo "  Date:   $(date)"
echo "========================================="

tg_send "Starting Viral Refseq BLAST database build" 2>/dev/null || true

####======================================####
####       BUILD VIRAL REFSEQ DATABASE    ####
####======================================####

mkdir -p "${VIRALREFSEQ_DB}"
cd "${VIRALREFSEQ_DB}"

echo "Downloading ${DB_NAME} into $(pwd)"

if update_blastdb.pl --decompress "${DB_NAME}"; then
    echo "Viral refseq database built successfully in ${VIRALREFSEQ_DB}"
    tg_send "Viral refseq DB build completed" 2>/dev/null || true
else
    EXIT_CODE=$?
    echo "BLAST Viral refseq database build FAILED (exit ${EXIT_CODE})" >&2
    tg_send "BLAST Viral refseq DB build FAILED" 2>/dev/null || true
    exit "${EXIT_CODE}"
fi

# Verification for BLAST database 
if ! blastdbcmd -db "${DB_NAME}" -info; then
    echo "WARNING: database verification failed" >&2
fi
