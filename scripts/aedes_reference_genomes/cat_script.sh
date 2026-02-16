#!/bin/bash

eval "$(/Users/Parsimony/miniconda3/bin/conda shell.bash hook)"

OUTDIR="./data/Aedes_RNA_sequences/genomes"

echo "=== Building Aedes Superreference  for Host Reads Depletion ==="

# Script for Aedes genome and chromosome annotation concatenation
# It is recommended to prefix chromosome names. This way we avoid naming conflicts

FILES=(
    # Aedes aegypti (RefSeq, chromosome-level)
    "GCF_002204515.2_AaegL5.0_genomic.fna|Aaeg|Aedes aegypti"
    # Aedes albopictus (RefSeq, chromosome-level, 2024)
    "GCF_035046485.1_AalbF5_genomic.fna|Aalb|Aedes albopictus"
    # Aedes mascarensis (GenBank, chromosome-level)
    "GCA_052575835.1_genomic.fna|Amas|Aedes mascarensis"
    # Aedes sierrensis (Ochlerotatus)
    "GCA_044231785.1_genomic.fna|Asie|Aedes sierrensis"
    # Aedes japonicus
    "GCA_052815935.1_genomic.fna|Ajap|Aedes japonicus"
    # Aedes notoscriptus
    "GCA_040801935.1_genomic.fna|Anot|Aedes notoscriptus"
    # Aedes koreicus
    "GCA_024533555.2_genomic.fna|Akor|Aedes koreicus"
)

echo "=== Concatenating genome files ==="

TEMP_FILE="${OUTDIR}/superreference_temp.fna"
> "$TEMP_FILE"

# Process each file
for entry in "${FILES[@]}"; do
    IFS='|' read -r filename prefix name <<< "$entry"
    
    echo "Processing: $name ($filename)"

    # Add prefix to chromosome names and append to temp file
    sed "s/^>/>${prefix}_/g" "$filename" >> "$TEMP_FILE"
done

echo "Superreference files created successfully"

echo "=== Generating statistics ==="

# Generate basic statistics for the superreference file
SEQ_COUNT=$(grep -c "^>" "$TEMP_FILE" 2>/dev/null || echo "0")

BASE_COUNT=0
while IFS= read -r line; do
    if [[ ! "$line" =~ ^\> ]]; then
        BASE_COUNT=$((BASE_COUNT + ${#line}))
    fi
done < "$TEMP_FILE"

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



