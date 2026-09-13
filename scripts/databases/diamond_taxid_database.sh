#!/bin/bash
# Author: Jorge Alberto Castro Rodríguez

# Database construction script for DIAMOND with taxonomy support
# Adapted for HPC environment with Apptainer containers
# Based on: https://github.com/JPbio/Metagenomics_dealing_with_blastn_diamond_results

# 02/09/2026
# Ver. 1.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Output directory for DIAMOND database (permanent storage)
DIAMOND_DB="${PERMANENT_BASE}/databases/diamond_db"

# Temporary directory for DIAMOND build (use fast scratch space)
TMP_DIR="${SCRATCH_BASE}/tmp/diamond_build"

# Download directory (use scratch for downloads)
DOWNLOAD_DIR="${SCRATCH_BASE}/downloads/diamond_db_build"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
DIAMOND_CONTAINER="${CONTAINERS}/diamond.sif"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot 
# source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

# Set number of threads based on available resources
THREADS=32

####==================================####
####          LOAD MODULES            ####
####==================================####

# Load required modules
module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          MAIN EXECUTION          ####
####==================================####

# Set strict error handling
set -euo pipefail

# Send start notification
# tg_send "Starting DIAMOND database construction for ${PROJECT_NAME}" 2>/dev/null

echo "========================================="
echo "DIAMOND Database Construction"
echo "Project: ${PROJECT_NAME}"
echo "Date: $(date)"
echo "========================================="

# Create required directories
mkdir -p "${DIAMOND_DB}" "${TMP_DIR}" "${DOWNLOAD_DIR}"
cd "${DOWNLOAD_DIR}"

echo "Using temporary directory: ${TMP_DIR}"
echo "Download directory: ${DOWNLOAD_DIR}"
echo "Database output directory: ${DIAMOND_DB}"

####==================================####
####          STEP 1: DIAMOND         ####
####==================================####

echo "[1/6] Setting up DIAMOND..."

# Check if DIAMOND container exists
if [ ! -f "${DIAMOND_CONTAINER}" ]; then
    echo "Warning: DIAMOND container not found at ${DIAMOND_CONTAINER}"
    echo "Attempting to use system DIAMOND..."
    if ! command -v diamond &> /dev/null; then
        echo "Error: DIAMOND not available. Please ensure container exists or install DIAMOND."
        exit 1
    fi
    DIAMOND_CMD="diamond"
else
    echo "Using DIAMOND container: ${DIAMOND_CONTAINER}"
    DIAMOND_CMD="apptainer exec ${DIAMOND_CONTAINER} diamond"
fi

# Verify DIAMOND works
if ! ${DIAMOND_CMD} --version &> /dev/null; then
    echo "Error: DIAMOND command failed. Please check your setup."
    exit 1
fi

DIAMOND_VERSION=$(${DIAMOND_CMD} --version 2>&1 | head -n 1)
echo "DIAMOND version: ${DIAMOND_VERSION}"

####==================================####
####          STEP 2: PROTEIN DB      ####
####==================================####

echo "[2/7] Setting up protein database (nr)..."

DB_NAME="nr"
DB_FASTA="${DOWNLOAD_DIR}/${DB_NAME}.gz"

if [ ! -f "${DB_FASTA}" ]; then
    echo "Downloading ${DB_NAME} database..."
    wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz -O "${DB_FASTA}"
else
    echo "${DB_NAME}.gz already exists, skipping download."
fi

# Verify download
if [ ! -f "${DB_FASTA}" ] || [ ! -s "${DB_FASTA}" ]; then
    echo "Error: Failed to download ${DB_NAME} database"
    exit 1
fi

echo "Database file size: $(du -h "${DB_FASTA}" | cut -f1)"

####==================================####
####          STEP 3: TAXONOMY        ####
####==================================####

echo "[3/7] Downloading NCBI taxonomy files..."

# Download taxdump files
TAXDUMP_FILE="new_taxdump.tar.gz"
if [ ! -f "${TAXDUMP_FILE}" ]; then
    echo "Downloading taxonomy dump..."
    if ! wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump.tar.gz -O "${TAXDUMP_FILE}"; then
        echo "Error: Failed to download taxonomy dump"
        exit 1
    fi
    tar -xzvf "${TAXDUMP_FILE}" nodes.dmp names.dmp
else
    echo "Taxonomy dump already exists, extracting..."
    tar -xzvf "${TAXDUMP_FILE}" nodes.dmp names.dmp
fi

# Download accession2taxid mapping
ACCESSION_MAP="prot.accession2taxid.FULL.gz"
if [ ! -f "${ACCESSION_MAP}" ]; then
    echo "Downloading accession2taxid mapping..."
    if ! wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz -O "${ACCESSION_MAP}"; then
        echo "Error: Failed to download accession2taxid mapping"
        exit 1
    fi
    echo "Decompressing accession2taxid mapping..."
    pigz -d -p ${THREADS} -k "${ACCESSION_MAP}"
else
    echo "Accession2taxid file exists, decompressing if needed..."
    if [ -f "${ACCESSION_MAP}" ] && [ ! -f "prot.accession2taxid.FULL" ]; then
        pigz -d -p ${THREADS} -k "${ACCESSION_MAP}"
    fi
fi

# Verify required taxonomy files exist
required_files=("nodes.dmp" "names.dmp" "prot.accession2taxid.FULL")
for file in "${required_files[@]}"; do
    if [ ! -f "${file}" ] || [ ! -s "${file}" ]; then
        echo "Error: Required taxonomy file ${file} is missing or empty"
        exit 1
    fi
done

echo "Taxonomy files verified"

####==================================####
####          STEP 4: BUILD DB        ####
####==================================####

echo "[4/7] Building DIAMOND database with taxonomy..."
DB_OUTPUT="${DIAMOND_DB}/${DB_NAME}_tax"

if [ -f "${DB_OUTPUT}.dmnd" ]; then
    echo "Database already exists at ${DB_OUTPUT}.dmnd"
    read -p "Do you want to rebuild it? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping database build."
    else
        echo "Removing existing database..."
        rm -f "${DB_OUTPUT}.dmnd"
        echo "Building new database..."
        ${DIAMOND_CMD} makedb \
            --in "${DB_FASTA}" \
            -d "${DB_OUTPUT}" \
            --taxonnodes nodes.dmp \
            --taxonnames names.dmp \
            --taxonmap prot.accession2taxid.FULL \
            --threads ${THREADS} \
            --tmpdir "${TMP_DIR}" \
            --verbose
    fi
else
    echo "Building DIAMOND database (this may take several hours)..."
    ${DIAMOND_CMD} makedb \
        --in "${DB_FASTA}" \
        -d "${DB_OUTPUT}" \
        --taxonnodes nodes.dmp \
        --taxonnames names.dmp \
        --taxonmap prot.accession2taxid.FULL \
        --threads ${THREADS} \
        --tmpdir "${TMP_DIR}" \
        --verbose
fi

# Verify database was created
if [ ! -f "${DB_OUTPUT}.dmnd" ]; then
    echo "Error: Database creation failed"
    exit 1
fi

DB_SIZE=$(du -h "${DB_OUTPUT}.dmnd" | cut -f1)
echo "Protein database created successfully! Size: ${DB_SIZE}"

####==================================####
####          STEP 5: NUCLEOTIDE DB   ####
####==================================####

echo "[5/7] Building nucleotide database (nt)..."

NT_DB_NAME="nt"
NT_FASTA="${DOWNLOAD_DIR}/${NT_DB_NAME}.gz"
NT_DB_OUTPUT="${DIAMOND_DB}/${NT_DB_NAME}_tax"

# Check if nucleotide database already exists
if [ -f "${NT_DB_OUTPUT}.ndb" ] || [ -f "${NT_DB_OUTPUT}.nhr" ]; then
    echo "Nucleotide database already exists at ${NT_DB_OUTPUT}"
    read -p "Do you want to rebuild it? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping nucleotide database build."
    else
        echo "Removing existing nucleotide database..."
        rm -f "${NT_DB_OUTPUT}".*
        BUILD_NT=true
    fi
else
    BUILD_NT=true
fi

if [ "${BUILD_NT:-false}" = true ]; then
    # Download nt database
    if [ ! -f "${NT_FASTA}" ]; then
        echo "Downloading nt database..."
        wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz -O "${NT_FASTA}"
    else
        echo "nt.gz already exists, skipping download."
    fi
    
    # Verify download
    if [ ! -f "${NT_FASTA}" ] || [ ! -s "${NT_FASTA}" ]; then
        echo "Error: Failed to download nt database"
        exit 1
    fi
    
    echo "Nucleotide database file size: $(du -h "${NT_FASTA}" | cut -f1)"
    
    # Download nucleotide taxonomy mapping
    NT_TAX_MAP="nucl_gb.accession2taxid.gz"
    if [ ! -f "${NT_TAX_MAP}" ]; then
        echo "Downloading nucleotide accession2taxid mapping..."
        if ! wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz -O "${NT_TAX_MAP}"; then
            echo "Warning: Failed to download nucleotide taxonomy mapping. Building without taxonomy."
            NT_TAX_MAP=""
        else
            echo "Decompressing nucleotide taxonomy mapping..."
            pigz -d -p ${THREADS} -k "${NT_TAX_MAP}"
            NT_TAX_MAP="${DOWNLOAD_DIR}/nucl_gb.accession2taxid"
        fi
    else
        if [ -f "${NT_TAX_MAP}" ] && [ ! -f "${DOWNLOAD_DIR}/nucl_gb.accession2taxid" ]; then
            pigz -d -p ${THREADS} -k "${NT_TAX_MAP}"
            NT_TAX_MAP="${DOWNLOAD_DIR}/nucl_gb.accession2taxid"
        else
            NT_TAX_MAP="${DOWNLOAD_DIR}/nucl_gb.accession2taxid"
        fi
    fi
    
    echo "Building nucleotide database with makeblastdb..."
    module load blast/2.15.0 2>/dev/null || echo "BLAST module not available"
    
    if command -v makeblastdb &> /dev/null; then
        if [ -n "${NT_TAX_MAP}" ] && [ -f "${NT_TAX_MAP}" ]; then
            makeblastdb -in "${NT_FASTA}" \
                -dbtype nucl \
                -out "${NT_DB_OUTPUT}" \
                -parse_seqids \
                -taxid_map "${NT_TAX_MAP}" \
                -threads ${THREADS}
        else
            makeblastdb -in "${NT_FASTA}" \
                -dbtype nucl \
                -out "${NT_DB_OUTPUT}" \
                -parse_seqids \
                -threads ${THREADS}
        fi
        
        # Verify nucleotide database was created
        if [ -f "${NT_DB_OUTPUT}.ndb" ]; then
            NT_DB_SIZE=$(du -h "${NT_DB_OUTPUT}".ndb | cut -f1)
            echo "Nucleotide database created successfully! Size: ${NT_DB_SIZE}"
        else
            echo "Warning: Nucleotide database creation may have failed"
        fi
    else
        echo "Warning: makeblastdb not available. Skipping nucleotide database build."
    fi
fi

####==================================####
####          STEP 6: TEST DB         ####
####==================================####

echo "[6/7] Testing protein database with query..."

# Create a test query (GFP protein sequence)
cat > test_query.fasta << 'EOF'
>test_query_GFP
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGVACFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
EOF

echo "Running test search..."
${DIAMOND_CMD} blastx \
    -q test_query.fasta \
    -d "${DB_OUTPUT}.dmnd" \
    -k 5 \
    -p 4 \
    -e 0.001 \
    -f 6 qseqid qlen qstart qend qcovhsp length sseqid slen sstart send scovhsp pident evalue bitscore qstrand qframe stitle staxids \
    -o test_output.tab \
    --quiet

if [ -f "test_output.tab" ] && [ -s "test_output.tab" ]; then
    echo "Test successful! Database is working correctly."
    echo "Sample results:"
    echo "----------------------------------------"
    head -n 3 test_output.tab | column -t
    echo "----------------------------------------"
else
    echo "Warning: Test failed or produced no results. Please check your database."
fi

####==================================####
####          STEP 7: CLEANUP         ####
####==================================####

echo "[7/7] Cleaning up temporary files..."

# Remove large downloaded files but keep the databases
rm -f "${DB_FASTA}"
rm -f "${NT_FASTA}"
rm -f "prot.accession2taxid.FULL" "prot.accession2taxid.FULL.gz"
rm -f "nucl_gb.accession2taxid" "nucl_gb.accession2taxid.gz"
rm -f "new_taxdump.tar.gz"
# Keep taxonomy files for reference
echo "Downloaded files cleaned. Taxonomy files (nodes.dmp, names.dmp) kept for reference."

# Clean temporary directory
echo "Cleaning temporary DIAMOND files..."
rm -rf "${TMP_DIR}"/*

####=========================####
####          FINISH         ####
####=========================####

echo "========================================="
echo "Database construction complete!"
echo "========================================="
echo "Protein database: ${DB_OUTPUT}.dmnd (${DB_SIZE})"
if [ -f "${NT_DB_OUTPUT}.ndb" ]; then
    echo "Nucleotide database: ${NT_DB_OUTPUT} (${NT_DB_SIZE:-unknown})"
fi
echo "Taxonomy files: ${DOWNLOAD_DIR}/{nodes.dmp,names.dmp}"
echo ""
echo "Usage examples:"
echo ""
echo "Protein search (DIAMOND blastx):"
echo "${DIAMOND_CMD} blastx -q your_proteins.fasta -d ${DB_OUTPUT}.dmnd -k 10 -p ${THREADS} -o results.tab"
echo ""
if [ -f "${NT_DB_OUTPUT}.ndb" ]; then
    echo "Nucleotide search (BLASTN):"
    echo "blastn -query your_contigs.fasta -db ${NT_DB_OUTPUT} -out results_nt.tab -outfmt 6 -num_threads ${THREADS}"
fi
echo "========================================="

# Create metadata file
METADATA_FILE="${DIAMOND_DB}/database_metadata_$(date +%Y%m%d).txt"
cat > "${METADATA_FILE}" << EOF
Database Build Information
==========================
Build date: $(date)
Build host: $(hostname)
Project: ${PROJECT_NAME}
DIAMOND version: ${DIAMOND_VERSION}

Protein Database (nr_tax)
-------------------------
File: ${DB_OUTPUT}.dmnd
Size: ${DB_SIZE}
Source: NCBI nr (protein)
Build command: ${DIAMOND_CMD} makedb --in ${DB_FASTA} -d ${DB_OUTPUT} --taxonnodes nodes.dmp --taxonnames names.dmp --taxonmap prot.accession2taxid.FULL --threads ${THREADS} --tmpdir ${TMP_DIR}

Nucleotide Database (nt_tax)
---------------------------
File: ${NT_DB_OUTPUT}
Size: ${NT_DB_SIZE:-Not built}
Source: NCBI nt (nucleotide)
Build command: makeblastdb -in ${NT_FASTA} -dbtype nucl -out ${NT_DB_OUTPUT} -parse_seqids -taxid_map ${NT_TAX_MAP:-None} -threads ${THREADS}

Taxonomy Information
-------------------
Taxonomy version: $(date -r nodes.dmp +%Y-%m-%d)
Taxonomy files used:
  - nodes.dmp: $(wc -l nodes.dmp | cut -f1) entries
  - names.dmp: $(wc -l names.dmp | cut -f1) entries
  - prot.accession2taxid.FULL: $(wc -l prot.accession2taxid.FULL | cut -f1) entries

System information:
$(uname -a)

EOF

echo "Metadata saved to: ${METADATA_FILE}"

# Send completion notification
# tg_send "Database construction complete for ${PROJECT_NAME}! Protein DB: ${DB_SIZE}, Nucleotide DB: ${NT_DB_SIZE:-Not built}" 2>/dev/null

# Exit successfully
exit 0
