#!/bin/bash
# The objetive of this script is to make a super reference based on nine Aedes species: 
# Ae. albopictus, Ae. aegypti, Ae. mascarensis, Ae. japonicus, Ae. notoscriptus, Ae. (Ochlerotatus) sierrensis, Ae. taeniorhynchus, Ae. koreicus, and Ae. polynesiensis.
# For details on the genomes used see .csv file in the folder

SECONDS=0

eval "$(/Users/Parsimony/miniconda3/bin/conda shell.bash hook)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # We assume that the scripts are located in a subdirectory "scripts/"
# Works regardless of the current working directory.
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Before running the script for genome indexing, we need to configure our script parameters 
THREADS=4
STAR_RAM=28000000000 
SUPER_REF="${PROJECT_ROOT}/data/references/mosquito_genomes/genomic_files/superreference_temp.fna"
INDEX_DIR="${PROJECT_ROOT}/data/references/aedes_super_index"
LOG_DIR="${PROJECT_ROOT}/logs"
GENOME_CHR_BIN_BITS=14   # For better memory distribution
SA_INDEX_BASES=13 # Parameters optimized for 12.1Gb genome (8 species)


# Create the index directory
mkdir -p "$INDEX_DIR"

# Remove any previous incomplete index
rm -rf "$INDEX_DIR"/*

echo "Building STAR index (this will take several hours)..."
echo "Parameters:"
echo "  • Threads: $THREADS"
echo "  • RAM limit: $(numfmt --to=iec $STAR_RAM)"
echo "  • GenomeChrBinNbits: $GENOME_CHR_BIN_BITS"
echo "  • SAindexNbases: $SA_INDEX_BASES"

conda activate star_env
# Build the STAR index
STAR \
  --runMode genomeGenerate \
  --genomeDir "$INDEX_DIR" \
  --genomeFastaFiles "$SUPER_REF" \
  --runThreadN "$THREADS" \
  --genomeChrBinNbits "$GENOME_CHR_BIN_BITS" \
  --limitGenomeGenerateRAM "$STAR_RAM" \
  --genomeSAindexNbases "$SA_INDEX_BASES"

conda deactivate

echo "=== SUPER REFERENCE CREATION COMPLETE ==="
echo "STAR index built in directory: $INDEX_DIR"
echo "Index contains these files:"
ls -la "$INDEX_DIR/"


echo "Time elapsed: $((SECONDS / 60)) minutes"