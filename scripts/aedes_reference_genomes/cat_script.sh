#!/bin/bash

eval "$(/Users/Parsimony/miniconda3/bin/conda shell.bash hook)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # We assume that the scripts are located in a subdirectory "scripts/"
# Works regardless of the current working directory.
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REFERENCE_DIR="${PROJECT_ROOT}/data/references/mosquito_genomes/genomic_files"
OUTDIR="${PROJECT_ROOT}/docs/aedes_genomes_specs"

mkdir -p "$REFERENCE_DIR" # Create the directory if it doesn't exist.

echo "=== BUILDING AEDES SUPERREFERENCE FOR HOST READ DEPLETION ==="

# Script for Aedes genome and chromosome annotation concatenation
# It is recommended to prefix chromosome names. This way we avoid naming conflicts

FILES=(
    # Aedes aegypti (RefSeq, chromosome-level)
    "GCF_002204515.2_AaegL5.0_genomic.fna|Aaeg|Aedes aegypti"
    # Aedes albopictus (RefSeq, chromosome-level, 2024)
    "GCF_035046485.1_AalbF5_genomic.fna|Aalb|Aedes albopictus"
    # Aedes mascarensis (GenBank, chromosome-level)
    "GCA_052575835.1_Am_MascCH02_pri1.0_genomic.fna|Amas|Aedes mascarensis"
    # Aedes sierrensis (Ochlerotatus)
    "GCA_044231785.1_Aedes_sierrensis_1_genomic.fna|Asie|Aedes sierrensis"
    # Aedes japonicus
    "GCA_052815935.1_ASM5281593v1_genomic.fna|Ajap|Aedes japonicus"
    # Aedes notoscriptus
    "GCA_040801935.1_CSIRO_AGI_Anoto_v1_genomic.fna|Anot|Aedes notoscriptus"
    # Aedes koreicus
    "GCA_024533555.2_Akor_1.1_genomic.fna|Akor|Aedes koreicus"
)

echo "=== CONCATENATING GENOME FILES ==="

TEMP_FILE="${REFERENCE_DIR}/superreference_temp.fna"
> "$TEMP_FILE"

# Process each file
for entry in "${FILES[@]}"; do
    IFS='|' read -r filename prefix name <<< "$entry"
    
    echo "Processing: $name ($filename)"

    # Add prefix to chromosome names and append to temp file
    sed "s/^>/>${prefix}_/g" "${REFERENCE_DIR}/${filename}" >> "$TEMP_FILE"
done

echo "Superreference files created successfully"

echo "=== SUMMARY REPORT ==="

# Generate basic statistics for the superreference file
SEQ_COUNT=$(grep -c "^>" "$TEMP_FILE" 2>/dev/null || echo "0")

BASE_COUNT=$(awk '!/^>/ {sum += length($0)} END {print sum}' "$TEMP_FILE" 2>/dev/null || echo "0")

echo "Total sequences: $SEQ_COUNT"
echo "Total bases: $BASE_COUNT"

# Save statistics to file
STATS_FILE="$OUTDIR/genome_stats.txt"
cat > "$STATS_FILE" << EOF
Created: $(date)
Files included: ${#FILES[@]}
Total sequences: $SEQ_COUNT
Total bases: $BASE_COUNT

Genomes:
EOF

# Add genome list to stats
for entry in "${FILES[@]}"; do
    IFS='|' read -r filename prefix name <<< "$entry"
    echo "  - $name ($filename)" >> "$STATS_FILE"
done



