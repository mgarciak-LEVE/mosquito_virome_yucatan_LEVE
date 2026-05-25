#!/bin/bash

# Script for viral phylogenetics analysis.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 24/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Input directories
CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes"
BLAST_DIR="${PROJECT_ROOT}/results/czid/BLAST"
PHYLO_OUTPUT="${PROJECT_ROOT}/results/czid/phylogenetics"

# OSFV database
OSFV_DB_DIR="${PROJECT_ROOT}/data/databases"
OSFV_REFERENCES="${OSFV_DB_DIR}/ncbi_dataset/data/genomic.fna"
OSFV_DB="${OSFV_DB_DIR}/osfv_references_db"

# Local viral protein database
VIRAL_DB="${PROJECT_ROOT}/data/databases/viral_refseq_prot"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
MAFFT_SIF="${CONTAINER_DIR}/mafft_7.221--0.sif"
TRIMAL_SIF="${CONTAINER_DIR}/trimal_1.5.1--h9948957_0.sif"
IQTREE_SIF="${CONTAINER_DIR}/iqtree_3.1.2--h8471819_0.sif"
BLAST_SIF="${CONTAINER_DIR}/blast_2.17.0--h66d330f_0.sif"
CDHIT_SIF="${CONTAINER_DIR}/cd-hit_4.8.1--h5ca1c30_13.sif"

THREADS=14

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

mkdir -p "${PHYLO_OUTPUT}"

echo "=========================================="
echo "Viral Phylogenetics Analysis"
echo "=========================================="
tg_send "Starting phylogenetics pipeline..."

####==================================####
####     OSFV: WHOLE GENOME (nt)      ####
####==================================####

echo "=========================================="
echo "OSFV: Whole Genome Phylogeny (nucleotide)"
echo "=========================================="

OSFV_OUT="${PHYLO_OUTPUT}/OSFV"
mkdir -p "${OSFV_OUT}"

echo "Collecting OSFV consensus genomes..."
cat "${CONSENSUS_DIR}/OSFV/"*_consensus.fa > "${OSFV_OUT}/osfv_consensus_all.fasta"

echo "Combining with reference genomes..."
cat "${OSFV_OUT}/osfv_consensus_all.fasta" \
    "${OSFV_REFERENCES}" \
    > "${OSFV_OUT}/osfv_all_sequences.fasta"

echo "Deduplicating sequences..."
apptainer exec "${CDHIT_SIF}" cd-hit-est \
    -i "${OSFV_OUT}/osfv_all_sequences.fasta" \
    -o "${OSFV_OUT}/osfv_nr.fasta" \
    -c 0.95 \
    -n 8 \
    -T "${THREADS}"

echo "Aligning with MAFFT..."
apptainer exec "${MAFFT_SIF}" mafft \
    --auto --reorder --thread "${THREADS}" \
    "${OSFV_OUT}/osfv_nr.fasta" \
    > "${OSFV_OUT}/osfv_aligned.fasta"

echo "Trimming alignment..."
apptainer exec "${TRIMAL_SIF}" trimal \
    -in "${OSFV_OUT}/osfv_aligned.fasta" \
    -out "${OSFV_OUT}/osfv_trimmed.fasta" \
    -gappyout

echo "Building phylogenetic tree with IQ-TREE..."
apptainer exec "${IQTREE_SIF}" iqtree \
    -s "${OSFV_OUT}/osfv_trimmed.fasta" \
    -m GTR+G+I -bb 1000 -alrt 1000 -nt AUTO \
    -pre "${OSFV_OUT}/osfv_tree"

echo "OSFV tree complete: ${OSFV_OUT}/osfv_tree.treefile"
tg_send "OSFV phylogeny complete."

####==================================####
####  TAV + GMMLV: RdRp PROTEIN (aa)  ####
####==================================####

echo "=========================================="
echo "TAV + GMMLV: Combined RdRp Phylogeny (aa)"
echo "=========================================="

COMBINED_OUT="${PHYLO_OUTPUT}/TAV_GMMLV"
mkdir -p "${COMBINED_OUT}"

# All RdRp ORF files
RDRP_FILES=(
    "${BLAST_DIR}/PM2469_S1_923325_split_orfs/PM2469_S1_923325_orf_6.fasta"
    "${BLAST_DIR}/PM2622_S32_923335_split_orfs/PM2622_S32_923335_orf_7.fasta"
    "${BLAST_DIR}/PM3183_S6_923339_split_orfs/PM3183_S6_923339_orf_1.fasta"
)

# Combine all RdRp sequences
echo "Combining all RdRp sequences..."
> "${COMBINED_OUT}/all_rdrp.fasta"
for rdrp_file in "${RDRP_FILES[@]}"; do
    if [ -f "$rdrp_file" ]; then
        cat "$rdrp_file" >> "${COMBINED_OUT}/all_rdrp.fasta"
    fi
done

if [ ! -s "${COMBINED_OUT}/all_rdrp.fasta" ]; then
    echo "ERROR: No RdRp files found."
    tg_send "ERROR: No RdRp files for combined phylogeny."
else
    # BLAST the longest RdRp against local viral DB
    echo "Finding reference sequences..."
    apptainer exec "${BLAST_SIF}" blastp \
        -query "${RDRP_FILES[0]}" \
        -db "$VIRAL_DB" \
        -out "${COMBINED_OUT}/combined_blast.tsv" \
        -outfmt "6 sseqid" \
        -max_target_seqs 100

    cut -f1 "${COMBINED_OUT}/combined_blast.tsv" | sort -u > "${COMBINED_OUT}/accessions.txt"

    echo "Downloading reference sequences..."
    > "${COMBINED_OUT}/references.fasta"
    while read acc; do
        apptainer exec "${BLAST_SIF}" blastdbcmd \
            -db "$VIRAL_DB" \
            -entry "$acc" \
            >> "${COMBINED_OUT}/references.fasta" 2>/dev/null
    done < "${COMBINED_OUT}/accessions.txt"

    # Combine all RdRp + references
    cat "${COMBINED_OUT}/all_rdrp.fasta" \
        "${COMBINED_OUT}/references.fasta" \
        > "${COMBINED_OUT}/combined_all.fasta"

    # Rename with clear labels
    awk '/^>/ {
        if ($0 ~ /Consensus_PM/) {
            match($0, /Consensus_(PM[^_]+_[^_]+_[^_]+).*_([0-9]+)/, arr)
            printf ">%s_orf%s\n", arr[1], arr[2]
        } else {
            gsub(/^>/, "")
            split($1, acc, " ")
            printf ">ref_%s\n", acc[1]
        }
        next
    }
    { print }' "${COMBINED_OUT}/combined_all.fasta" > "${COMBINED_OUT}/combined_renamed.fasta"

    # Deduplicate
    echo "Deduplicating sequences..."
    apptainer exec "${CDHIT_SIF}" cd-hit \
        -i "${COMBINED_OUT}/combined_renamed.fasta" \
        -o "${COMBINED_OUT}/combined_dedup.fasta" \
        -c 0.95 \
        -T "${THREADS}"

    echo "Sequences after deduplication: $(grep -c '^>' ${COMBINED_OUT}/combined_dedup.fasta)"

    # Align
    echo "Aligning with MAFFT..."
    apptainer exec "${MAFFT_SIF}" mafft \
        --auto --reorder \
        "${COMBINED_OUT}/combined_dedup.fasta" \
        > "${COMBINED_OUT}/combined_aligned.fasta"

    # Trim
    echo "Trimming alignment..."
    apptainer exec "${TRIMAL_SIF}" trimal \
        -in "${COMBINED_OUT}/combined_aligned.fasta" \
        -out "${COMBINED_OUT}/combined_trimmed.fasta" \
        -gappyout

    # Build tree
    echo "Building phylogenetic tree with IQ-TREE..."
    apptainer exec "${IQTREE_SIF}" iqtree \
        -s "${COMBINED_OUT}/combined_trimmed.fasta" \
        -m MFP -bb 1000 -alrt 1000 -nt AUTO \
        -pre "${COMBINED_OUT}/combined_tree"

    echo "Combined tree complete: ${COMBINED_OUT}/combined_tree.treefile"
    tg_send "TAV + GMMLV phylogeny complete."
fi

echo ""
echo "=========================================="
echo "Phylogenetics Pipeline Complete"
echo "=========================================="
echo "Trees: ${PHYLO_OUTPUT}/*/tree.treefile"
tg_send "Phylogenetics complete."