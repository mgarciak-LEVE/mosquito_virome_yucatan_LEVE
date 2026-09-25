#!/bin/bash
# Build UniRef90 DIAMOND database with taxonomy for VITAP v1.12
# Author: Jorge Alberto Castro Rodriguez
# Ver. 1.0.0
# 20/09/2026

set -euo pipefail

####==================================####
####           CONFIGURATION          ####
####==================================####

# Output directory (must have >=300GB free)
WORK_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/databases/vitap_db/db"
OUTPUT_DB="${WORK_DIR}/uniref90.dmnd"

# VITAP container (provides diamond 2.1.16)
VITAP_CONTAINER="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers/vitap_1.12--pyhdfd78af_0.sif"

# Threads for diamond makedb
THREADS=12

# Disk space check (GB)
MIN_FREE_GB=300

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          DISK CHECK             ####
####==================================####

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

AVAIL_KB=$(df -P . | awk 'NR==2 {print $4}')
AVAIL_GB=$((AVAIL_KB / 1024 / 1024))
echo "Available disk space: ${AVAIL_GB} GB"

if (( AVAIL_GB < MIN_FREE_GB )); then
    echo "ERROR: Need at least ${MIN_FREE_GB} GB free, only ${AVAIL_GB} GB available."
    echo "Aborting."
    exit 1
fi

####==================================####
####       STEP 1: DOWNLOAD          ####
####==================================####

echo "========================================="
echo "  Step 1: Downloading source files"
echo "========================================="

# UniRef90 FASTA
if [[ ! -f "uniref90.fasta.gz" ]]; then
    echo "Downloading UniRef90 FASTA (~43 GB)..."
    wget -c -O uniref90.fasta.gz \
        "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
else
    echo "uniref90.fasta.gz already exists, skipping download."
fi

# NCBI protein accession2taxid
if [[ ! -f "prot.accession2taxid.gz" ]]; then
    echo "Downloading prot.accession2taxid.gz..."
    wget -c -O prot.accession2taxid.gz \
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
else
    echo "prot.accession2taxid.gz already exists, skipping download."
fi

# NCBI taxonomy dump (nodes.dmp, names.dmp)
if [[ ! -f "nodes.dmp" ]] || [[ ! -f "names.dmp" ]]; then
    echo "Downloading NCBI taxonomy dump..."
    wget -c -O taxdmp.zip \
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
    unzip -o taxdmp.zip nodes.dmp names.dmp
else
    echo "nodes.dmp and names.dmp already exist, skipping."
fi

####==================================####
####       STEP 2: DECOMPRESS        ####
####==================================####

echo "========================================="
echo "  Step 2: Decompressing UniRef90 FASTA"
echo "========================================="

if [[ ! -f "uniref90.fasta" ]]; then
    echo "Decompressing uniref90.fasta.gz (this takes a while)..."
    gunzip -k uniref90.fasta.gz
else
    echo "uniref90.fasta already exists, skipping."
fi

####==================================####
####       STEP 3: BUILD DB          ####
####==================================####

echo "========================================="
echo "  Step 3: Building DIAMOND database"
echo "========================================="

apptainer exec \
    --bind "${WORK_DIR}":/work \
    "${VITAP_CONTAINER}" \
    diamond makedb \
        --in /work/uniref90.fasta \
        -d /work/uniref90.dmnd \
        --taxonmap /work/prot.accession2taxid.gz \
        --taxonnodes /work/nodes.dmp \
        --taxonnames /work/names.dmp \
        --threads "${THREADS}"

echo "========================================="
echo "  DONE"
echo "  Database: ${OUTPUT_DB}"
echo "========================================="
