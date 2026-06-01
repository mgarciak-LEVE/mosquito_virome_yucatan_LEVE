#!/bin/bash
# Script: Ochlerotatus scapularis flavivirus (OSFV) phylogenetics.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.1.0
# 28/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Input: improved consensus from the completion script
COMPLETION_DIR="${PROJECT_ROOT}/results/czid/osfv_complete_genome"
ALL_CONSENSUS="${COMPLETION_DIR}/all_osfv_consensus.fasta"
PHYLO_OUT="${PROJECT_ROOT}/results/czid/osfv_complete_genome/phylogeny"

OSFV_DB_DIR="${PROJECT_ROOT}/data/databases"
OSFV_REFERENCES="${OSFV_DB_DIR}/ncbi_dataset/data/genomic.fna"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
MAFFT_SIF="${CONTAINER_DIR}/mafft_7.221--0.sif"
TRIMAL_SIF="${CONTAINER_DIR}/trimal_1.5.1--h9948957_0.sif"
IQTREE_SIF="${CONTAINER_DIR}/iqtree_3.1.2--h8471819_0.sif"
CDHIT_SIF="${CONTAINER_DIR}/cd-hit_4.8.1--h5ca1c30_13.sif"

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

THREADS=28
mkdir -p "${PHYLO_OUT}"

echo "=========================================="
echo "OSFV: Phylogenetics"
echo "=========================================="
tg_send "Starting OSFV phylogenetics..."

# Check if improved consensus exists
if [ ! -f "$ALL_CONSENSUS" ]; then
    echo "ERROR: ${ALL_CONSENSUS} not found. Run genome completion first."
    tg_send "ERROR: OSFV consensus not found."
    exit 1
fi

# Combine all improved consensus + references for phylogeny
echo "Combining improved genomes with references..."
cat "$ALL_CONSENSUS" "$OSFV_REFERENCES" > "${PHYLO_OUT}/osfv_all.fasta"

# Deduplicate
echo "Deduplicating sequences..."
apptainer exec "${CDHIT_SIF}" cd-hit-est \
    -i "${PHYLO_OUT}/osfv_all.fasta" \
    -o "${PHYLO_OUT}/osfv_nr.fasta" \
    -c 0.95 -n 8 -T "${THREADS}"

echo "Sequences for phylogeny: $(grep -c '^>' ${PHYLO_OUT}/osfv_nr.fasta)"

# Align
echo "Aligning with MAFFT..."
apptainer exec "${MAFFT_SIF}" mafft --auto --maxiterate 1000 --reorder --thread "${THREADS}" \
    "${PHYLO_OUT}/osfv_nr.fasta" > "${PHYLO_OUT}/osfv_aligned.fasta"

# Convert to uppercase
echo "Converting to uppercase..."
awk '!/^>/ {print toupper($0)} /^>/ {print}' \
    "${PHYLO_OUT}/osfv_aligned.fasta" \
    > "${PHYLO_OUT}/osfv_aligned_upper.fasta"

# Trim
echo "Trimming alignment..."
apptainer exec "${TRIMAL_SIF}" trimal \
    -in "${PHYLO_OUT}/osfv_aligned_upper.fasta" \
    -out "${PHYLO_OUT}/osfv_trimmed.fasta" \
    -gt 0.2 \
    -st 0.001

# Tree
echo "Building phylogenetic tree with IQ-TREE..."
apptainer exec "${IQTREE_SIF}" iqtree \
    -redo \
    -s "${PHYLO_OUT}/osfv_trimmed.fasta" \
    -m GTR+G+I -bb 1000 -alrt 1000 -nt AUTO \
    -pre "${PHYLO_OUT}/osfv_improved_tree"

echo ""
echo "Tree: ${PHYLO_OUT}/osfv_improved_tree.treefile"
echo "Complete."
tg_send "OSFV phylogeny complete."