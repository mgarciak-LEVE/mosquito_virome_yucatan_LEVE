#!/bin/bash

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
# Input directory - where FASTQ files are
INPUT_DIR="${PROJECT_SCRATCH}/data/references/mosquito_genomes/genomic_files"
# Output directory for validation results (on Lustre)
OUTPUT_DIR="${PROJECT_SCRATCH}/data/references/mosquito_genomes/aedes_super_index"

mkdir -p "${INPUT_DIR}" "${OUTPUT_DIR}" "${STATS_DIR}"

# Genome statistics output directory
STATS_DIR="${PROJECT_SCRATCH}/docs/aedes_genomes_specs"

mkdir -p "${INPUT_DIR}" "${OUTPUT_DIR}" "${STATS_DIR}"

# Scripts directory (on NFS)
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot (source from NFS)
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load entrez-direct/21.6--he881be0_0 2>/dev/null || echo "entrez-direct not loaded"
module load samtools/1.20--h50ea8bc_0 2>/dev/null || echo "samtools not loaded"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Mosquito Genomes Concatenation Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Input: ${INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting genome concatenation for ${PROJECT_NAME}" 2>/dev/null || true

####==================================####
####          GENOME DOWNLOAD         ####
####==================================####

GENOMES=(
    "GCF_002204515.2" # Aedes aegypti (RefSeq, chromosome-level)
    "GCF_035046485.1" # Aedes albopictus (RefSeq, chromosome-level, 2024)
    "GCA_052575835.1" # Aedes mascarensis (GenBank, chromosome-level)
    "GCA_044231785.1" # Aedes sierrensis (Ochlerotatus)
    "GCA_052815935.1" # Aedes japonicus
    "GCA_040801935.1" # Aedes notoscriptus
    "GCA_024533555.2" # Aedes koreicus
)

for genome in "${GENOMES[@]}"; do
    echo "Downloading genome: $genome"
    tg_send "Downloading genome: $genome" 2>/dev/null || true
    # Use NCBI's efetch to download the genomic FASTA file
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../${genome}_*_genomic.fna.gz" -O "${INPUT_DIR}/${genome}.fna.gz"
    gunzip -f "${INPUT_DIR}/${genome}.fna.gz"
done

echo "Download completed!"
tg_send "Genome download completed" 2>/dev/null || true

####==================================####
####         RENAME GENOME FILES      ####
####==================================####

echo ""
echo "=== RENAMING GENOME FILES ==="

# Rename files to their actual genome names
# This uses the actual file names from NCBI

# Aedes aegypti (GCF_002204515.2) → GCF_002204515.2_AaegL5.0_genomic.fna
if [[ -f "${INPUT_DIR}/GCF_002204515.2.fna" ]]; then
    mv "${INPUT_DIR}/GCF_002204515.2.fna" "${INPUT_DIR}/GCF_002204515.2_AaegL5.0_genomic.fna"
fi

# Aedes albopictus (GCF_035046485.1) → GCF_035046485.1_AalbF5_genomic.fna
if [[ -f "${INPUT_DIR}/GCF_035046485.1.fna" ]]; then
    mv "${INPUT_DIR}/GCF_035046485.1.fna" "${INPUT_DIR}/GCF_035046485.1_AalbF5_genomic.fna"
fi 

# Aedes mascarensis (GCA_052575835.1) → GCA_052575835.1_Am_MascCH02_pri1.0_genomic.fna
if [[ -f "${INPUT_DIR}/GCA_052575835.1.fna" ]]; then
    mv "${INPUT_DIR}/GCA_052575835.1.fna" "${INPUT_DIR}/GCA_052575835.1_Am_MascCH02_pri1.0_genomic.fna"
fi

# Aedes sierrensis (GCA_044231785.1) → GCA_044231785.1_Aedes_sierrensis_1_genomic.fna
if [[ -f "${INPUT_DIR}/GCA_044231785.1.fna" ]]; then
    mv "${INPUT_DIR}/GCA_044231785.1.fna" "${INPUT_DIR}/GCA_044231785.1_Aedes_sierrensis_1_genomic.fna"
fi

# Aedes japonicus (GCA_052815935.1) → GCA_052815935.1_ASM5281593v1_genomic.fna
if [[ -f "${INPUT_DIR}/GCA_052815935.1.fna" ]]; then
    mv "${INPUT_DIR}/GCA_052815935.1.fna" "${INPUT_DIR}/GCA_052815935.1_ASM5281593v1_genomic.fna"
fi

# Aedes notoscriptus (GCA_040801935.1) → GCA_040801935.1_CSIRO_AGI_Anoto_v1_genomic.fna
if [[ -f "${INPUT_DIR}/GCA_040801935.1.fna" ]]; then
    mv "${INPUT_DIR}/GCA_040801935.1.fna" "${INPUT_DIR}/GCA_040801935.1_CSIRO_AGI_Anoto_v1_genomic.fna"
fi

# Aedes koreicus (GCA_024533555.2) → GCA_024533555.2_Akor_1.1_genomic.fna
if [[ -f "${INPUT_DIR}/GCA_024533555.2.fna" ]]; then
    mv "${INPUT_DIR}/GCA_024533555.2.fna" "${INPUT_DIR}/GCA_024533555.2_Akor_1.1_genomic.fna"
fi

echo "Renaming completed!"

####==================================####
####       GENOME CONCATENATION       ####
####==================================####

echo ""
echo "=== CONCATENATING GENOME FILES ==="

# Define files with their prefixes
declare -a FILES=(
    "GCF_002204515.2_AaegL5.0_genomic.fna|Aaeg|Aedes aegypti"
    "GCF_035046485.1_AalbF5_genomic.fna|Aalb|Aedes albopictus"
    "GCA_052575835.1_Am_MascCH02_pri1.0_genomic.fna|Amas|Aedes mascarensis"
    "GCA_044231785.1_Aedes_sierrensis_1_genomic.fna|Asie|Aedes sierrensis"
    "GCA_052815935.1_ASM5281593v1_genomic.fna|Ajap|Aedes japonicus"
    "GCA_040801935.1_CSIRO_AGI_Anoto_v1_genomic.fna|Anot|Aedes notoscriptus"
    "GCA_024533555.2_Akor_1.1_genomic.fna|Akor|Aedes koreicus"
)

# Create temp file for concatenation
TEMP_FILE="${INPUT_DIR}/superreference_temp.fna"
> "$TEMP_FILE"

# Check if files exist before processing
processed_count=0
for entry in "${FILES[@]}"; do
    IFS='|' read -r filename prefix name <<< "$entry"
    
    if [[ -f "${INPUT_DIR}/${filename}" ]]; then
        echo "Processing: ${name} (${filename})"
        tg_send "Processing: ${name}" 2>/dev/null || true
        
        awk -v prefix="$prefix" '{if ($0 ~ /^>/) print ">" prefix "_" substr($0, 2); else print}' "${INPUT_DIR}/${filename}" >> "$TEMP_FILE"
        
        ((processed_count++))
    else
        echo "  File not found: ${filename}"
        tg_send "File not found: ${filename}" 2>/dev/null || true
    fi
done

echo "Processed ${processed_count} of ${#FILES[@]} files"

# Move the superreference to output directory
if [[ -s "$TEMP_FILE" ]]; then
    mv "$TEMP_FILE" "${OUTPUT_DIR}/superreference.fna"
    echo "Superreference created: ${OUTPUT_DIR}/superreference.fna"
    tg_send "Superreference created" 2>/dev/null || true
else
    echo "Superreference creation failed (empty file)"
    tg_send "Superreference creation failed" 2>/dev/null || true
    exit 1
fi

####==================================####
####           GENERATE STATS         ####
####==================================####

echo ""
echo "=== GENERATING STATISTICS ==="

SUPERREF="${OUTPUT_DIR}/superreference.fna"

# Count sequences
SEQ_COUNT=$(grep -c "^>" "$SUPERREF" 2>/dev/null || echo "0")

# Count total bases
BASE_COUNT=$(awk '!/^>/ {sum += length($0)} END {print sum}' "$SUPERREF" 2>/dev/null || echo "0")

# Calculate average sequence length
AVG_LEN=$((BASE_COUNT / SEQ_COUNT))

echo "Total sequences: ${SEQ_COUNT}"
echo "Total bases: ${BASE_COUNT}"
echo "Average sequence length: ${AVG_LEN}"

# Save statistics to file
STATS_FILE="${STATS_DIR}/genome_stats.txt"
cat > "$STATS_FILE" << EOF
========================================
  Mosquito Superreference Statistics
========================================
Created: $(date)
Project: ${PROJECT_NAME}

Summary:
  - Total sequences: ${SEQ_COUNT}
  - Total bases: ${BASE_COUNT}
  - Average sequence length: ${AVG_LEN}
  - Files processed: ${processed_count}

Genomes included:
EOF

for entry in "${FILES[@]}"; do
    IFS='|' read -r filename prefix name <<< "$entry"
    if [[ -f "${INPUT_DIR}/${filename}" ]]; then
        echo "  ${name} (${filename})" >> "$STATS_FILE"
    else
        echo "  ${name} (${filename} - NOT FOUND)" >> "$STATS_FILE"
    fi
done

echo "Stats saved to: ${STATS_FILE}"
tg_send "Stats saved to: ${STATS_FILE}" 2>/dev/null || true

####==================================####
####              SUMMARY             ####
####==================================####

echo ""
echo "========================================="
echo "  Mosquito Genomes Concatenation Complete"
echo "  Superreference: ${OUTPUT_DIR}/superreference.fna"
echo "  Statistics: ${STATS_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Genome concatenation complete for ${PROJECT_NAME}" 2>/dev/null || true
