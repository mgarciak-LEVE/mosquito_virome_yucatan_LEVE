#!/bin/bash
# Author: Jorge Alberto Castro Rodríguez
# DIAMOND nr_tax database construction (protein, with taxonomy)
# Based on: https://github.com/JPbio/Metagenomics_dealing_with_blastn_diamond_results
# Ver. 1.2.0

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
THREADS=40
FORCE_REBUILD="${FORCE_REBUILD:-true}"

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

# --- Apptainer wrapper ---
APPTAINER_BIND="/lustre:/lustre"
DIAMOND_CMD="apptainer exec --pwd / --bind ${APPTAINER_BIND} ${DIAMOND_CONTAINER} diamond"

# Verify container
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
echo "DIAMOND_DB:     ${DIAMOND_DB}"
echo "DOWNLOAD_DIR:   ${DOWNLOAD_DIR}"
echo "TMP_DIR:        ${TMP_DIR}"
echo "========================================="

####==================================####
####  STEP 1: DOWNLOAD nr (protein)   ####
####==================================####
echo "[1/5] Downloading NCBI nr (protein)..."
NR_GZ="${DOWNLOAD_DIR}/nr.gz"
NR_FASTA="${DOWNLOAD_DIR}/nr"

if [ ! -f "${NR_GZ}" ]; then
    wget -c -O "${NR_GZ}" "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
else
    echo "nr.gz already present, skipping download."
fi
[ -s "${NR_GZ}" ] || { echo "Error: nr.gz is missing or empty"; exit 1; }
echo "nr.gz size: $(du -h "${NR_GZ}" | cut -f1)"

if [ ! -f "${NR_FASTA}" ]; then
    echo "Decompressing nr.gz (this takes a while)..."
    pigz -d -p "${THREADS}" -k "${NR_GZ}"
fi
[ -s "${NR_FASTA}" ] || { echo "Error: nr FASTA is missing or empty"; exit 1; }
echo "nr FASTA size: $(du -h "${NR_FASTA}" | cut -f1)"

####==================================####
####  STEP 2: TAXONOMY FILES          ####
####==================================####
echo "[2/5] Downloading NCBI taxonomy files..."

if [ ! -f nodes.dmp ] || [ ! -f names.dmp ]; then
    wget -c -O new_taxdump.tar.gz \
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
    tar -xzf new_taxdump.tar.gz nodes.dmp names.dmp
fi

if [ ! -f prot.accession2taxid ]; then
    wget -c -O prot.accession2taxid.gz \
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
    pigz -d -p "${THREADS}" prot.accession2taxid.gz
fi

for f in nodes.dmp names.dmp prot.accession2taxid; do
    [ -s "$f" ] || { echo "Error: required taxonomy file $f missing/empty"; exit 1; }
done
echo "Taxonomy files verified."

####==================================####
####  STEP 3: BUILD DIAMOND nr_tax    ####
####==================================####
echo "[3/5] Building DIAMOND nr_tax database..."
TAXDIR="${DOWNLOAD_DIR}"
NR_TAX="${DIAMOND_DB}/nr_tax"

# Ensure the tmpdir exists inside the container path (it does because /lustre is bound 1:1)
mkdir -p "${TMP_DIR}"

if [ -f "${NR_TAX}.dmnd" ] && [ "${FORCE_REBUILD}" != "true" ]; then
    echo "nr_tax.dmnd already exists. Set FORCE_REBUILD=true to rebuild."
else
    # Move any existing DB aside so we don't get confused mid-build
    [ -f "${NR_TAX}.dmnd" ] && mv "${NR_TAX}.dmnd" "${NR_TAX}.dmnd.old.$(date +%s)"
    echo "Building (this can take 1-2 hours for nr; ~500 GB peak disk)..."
    ${DIAMOND_CMD} makedb \
        --in        "${NR_FASTA}" \
        -d          "${NR_TAX}" \
        --taxonnodes "${TAXDIR}/nodes.dmp" \
        --taxonnames "${TAXDIR}/names.dmp" \
        --taxonmap   "${TAXDIR}/prot.accession2taxid" \
        --threads   "${THREADS}" \
        --tmpdir    "${TMP_DIR}" \
        --log
fi

[ -f "${NR_TAX}.dmnd" ] || { echo "Error: DIAMOND DB build failed"; exit 1; }
echo "nr_tax.dmnd size: $(du -h "${NR_TAX}.dmnd" | cut -f1)"

####==================================####
####  STEP 4: TEST (real taxonomy check) ####
####==================================####
echo "[4/5] Testing DIAMOND nr_tax with a protein query..."

cat > test_query.fasta << 'EOF'
>test_GFP
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGVACFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
EOF

${DIAMOND_CMD} blastp \
    -q test_query.fasta \
    -d "${NR_TAX}.dmnd" \
    -k 5 \
    -p 4 \
    -e 0.001 \
    -f 6 qseqid sseqid pident length evalue bitscore staxids sscinames \
    -o test_diamond.tab \
    --quiet || true

if [ ! -s test_diamond.tab ]; then
    echo "ERROR: test returned no hits — DB build is broken."
    exit 1
fi

echo "Test hits:"
echo "----------------------------------------"
column -t test_diamond.tab | head -n 5
echo "----------------------------------------"

# Real taxonomy verification: check that at least one row has a non-empty taxid
N_WITH_TAXID=$(awk -F'\t' '$6 != "" && $6 != "N/A"' test_diamond.tab | wc -l)
N_TOTAL=$(wc -l < test_diamond.tab)
echo "Taxonomy check: ${N_WITH_TAXID}/${N_TOTAL} test hits have a taxid."

if [ "${N_WITH_TAXID}" -eq 0 ]; then
    echo "ERROR: DB has no taxonomy! Check that --taxonmap/--taxonnodes/--taxonnames were passed."
    echo "Inspect diamond_makedb.log for 'Loading taxonomy mapping file...'."
    exit 1
fi

echo "SUCCESS: DB contains taxonomy and produces staxids."

####==================================####
####  STEP 5: CLEANUP                 ####
####==================================####
echo "[5/5] Cleaning up..."

rm -f "${NR_FASTA}"
rm -f prot.accession2taxid.gz
rm -f new_taxdump.tar.gz
rm -f test_query.fasta

echo "Cleaning DIAMOND temp dir..."
rm -rf "${TMP_DIR:?}"/*

echo "========================================="
echo "DONE. Database: ${NR_TAX}.dmnd"
echo "========================================="
