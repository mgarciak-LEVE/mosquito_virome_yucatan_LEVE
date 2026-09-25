#!/bin/bash
# Author: Jorge Alberto Castro Rodríguez
# DIAMOND nr_tax database construction (protein, with taxonomy)
# Based on: https://github.com/JPbio/Metagenomics_dealing_with_blastn_diamond_results
# Ver. 1.1.0

####==================================####
####           CONFIGURATION          ####
####==================================####
set -euo pipefail

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- Storage ---
PERMANENT_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46"
DB_BASE="${PERMANENT_BASE}/databases"
DIAMOND_DB="${DB_BASE}/diamond_db"
SCRATCH_BASE="${PERMANENT_BASE}/projects"
TMP_DIR="${SCRATCH_BASE}/tmp/diamond_build"
DOWNLOAD_DIR="${SCRATCH_BASE}/downloads/diamond_db_build"

# --- Container ---
CONTAINERS="${PERMANENT_BASE}/containers"
DIAMOND_CONTAINER="${CONTAINERS}/diamond.sif"

# --- Runtime ---
THREADS=32
FORCE_REBUILD="${FORCE_REBUILD:-false}"

####==================================####
####          LOAD MODULES            ####
####==================================####
module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          MAIN EXECUTION          ####
####==================================####
echo "========================================="
echo "DIAMOND nr_tax Database Construction"
echo "Project: ${PROJECT_NAME}"
echo "Date: $(date)"
echo "========================================="

mkdir -p "${DIAMOND_DB}" "${TMP_DIR}" "${DOWNLOAD_DIR}"
cd "${DOWNLOAD_DIR}"

# --- Apptainer wrapper (bind /lustre + work dirs) ---
APPTAINER_BIND="/lustre:/lustre,${DOWNLOAD_DIR}:${DOWNLOAD_DIR},${TMP_DIR}:${TMP_DIR}"
DIAMOND_CMD="apptainer exec --bind ${APPTAINER_BIND} ${DIAMOND_CONTAINER} diamond"

# Verify DIAMOND
if [ ! -f "${DIAMOND_CONTAINER}" ]; then
    echo "Error: DIAMOND container not found at ${DIAMOND_CONTAINER}"
    exit 1
fi
if ! ${DIAMOND_CMD} --version >/dev/null 2>&1; then
    echo "Error: DIAMOND could not be executed inside the container."
    exit 1
fi
DIAMOND_VERSION="$(${DIAMOND_CMD} --version 2>&1 | head -n 1 || true)"
echo "DIAMOND version: ${DIAMOND_VERSION}"

####==================================####
####  STEP 1: DOWNLOAD nr (protein)   ####
####==================================####
echo "[1/5] Downloading NCBI nr (protein)..."
NR_GZ="${DOWNLOAD_DIR}/nr.gz"
NR_FASTA="${DOWNLOAD_DIR}/nr"

if [ ! -f "${NR_GZ}" ]; then
    wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz -O "${NR_GZ}"
else
    echo "nr.gz already present, skipping download."
fi
[ -s "${NR_GZ}" ] || { echo "Error: nr.gz is missing or empty"; exit 1; }
echo "nr.gz size: $(du -h "${NR_GZ}" | cut -f1)"

# Decompress (keep .gz so re-runs can skip download)
if [ ! -f "${NR_FASTA}" ]; then
    echo "Decompressing nr.gz (this takes a while)..."
    pigz -d -p ${THREADS} -k "${NR_GZ}"
fi
[ -s "${NR_FASTA}" ] || { echo "Error: nr FASTA is missing or empty"; exit 1; }
echo "nr FASTA size: $(du -h "${NR_FASTA}" | cut -f1)"

####==================================####
####  STEP 2: TAXONOMY FILES          ####
####==================================####
echo "[2/5] Downloading NCBI taxonomy files..."

# --- nodes.dmp / names.dmp (from new_taxdump) ---
if [ ! -f nodes.dmp ] || [ ! -f names.dmp ]; then
    wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
         -O new_taxdump.tar.gz
    tar -xzf new_taxdump.tar.gz nodes.dmp names.dmp
fi

# --- prot.accession2taxid (non-FULL is ~5GB vs ~100GB for FULL) ---
if [ ! -f prot.accession2taxid ]; then
    wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz \
         -O prot.accession2taxid.gz
    pigz -d -p ${THREADS} prot.accession2taxid.gz
fi

# --- Verify ---
for f in nodes.dmp names.dmp prot.accession2taxid; do
    [ -s "$f" ] || { echo "Error: required taxonomy file $f missing/empty"; exit 1; }
done
echo "Taxonomy files verified."

####==================================####
####  STEP 3: BUILD DIAMOND nr_tax    ####
####==================================####
echo "[3/5] Building DIAMOND nr_tax database..."
NR_TAX="${DIAMOND_DB}/nr_tax"

if [ -f "${NR_TAX}.dmnd" ] && [ "${FORCE_REBUILD}" != "true" ]; then
    echo "nr_tax.dmnd already exists. Set FORCE_REBUILD=true to rebuild."
else
    [ -f "${NR_TAX}.dmnd" ] && rm -f "${NR_TAX}.dmnd"
    echo "Building (this can take many hours and ~1TB RAM/disk)..."
    ${DIAMOND_CMD} makedb \
        --in "${NR_FASTA}" \
        -d "${NR_TAX}" \
        --taxonnodes nodes.dmp \
        --taxonnames names.dmp \
        --taxonmap prot.accession2taxid \
        --threads ${THREADS} \
        --tmpdir "${TMP_DIR}" \
        --verbose
fi

[ -f "${NR_TAX}.dmnd" ] || { echo "Error: DIAMOND DB build failed"; exit 1; }
echo "nr_tax.dmnd size: $(du -h "${NR_TAX}.dmnd" | cut -f1)"

####==================================####
####  STEP 4: TEST                    ####
####==================================####
echo "[4/5] Testing DIAMOND nr_tax with a protein query..."

cat > test_query.fasta << 'EOF'
>test_GFP
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGVACFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
EOF

# Protein query -> blastp (NOT blastx)
${DIAMOND_CMD} blastp \
    -q test_query.fasta \
    -d "${NR_TAX}.dmnd" \
    -k 5 \
    -p 4 \
    -e 0.001 \
    -f 6 qseqid sseqid pident length evalue bitscore staxids stitle \
    -o test_diamond.tab \
    --quiet || true

if [ -s test_diamond.tab ]; then
    echo "Test OK — top hits:"
    echo "----------------------------------------"
    head -n 3 test_diamond.tab | column -t
    echo "----------------------------------------"
else
    echo "WARNING: test returned no hits (check DB or query)."
fi

####==================================####
####  STEP 5: CLEANUP                 ####
####==================================####
echo "[5/5] Cleaning up..."

# Keep nr.gz (re-download avoidance) but drop the huge uncompressed FASTA
rm -f "${NR_FASTA}"
# Drop the compressed accession map (uncompressed version kept for reference)
rm -f prot.accession2taxid.gz
# Drop taxdump archive
rm -f new_taxdump.tar.gz

echo "Cleaning DIAMOND temp dir..."
rm -rf "${TMP_DIR:?}"/*

echo "========================================="
echo "DONE. Database: ${NR_TAX}.dmnd"
echo "========================================="
