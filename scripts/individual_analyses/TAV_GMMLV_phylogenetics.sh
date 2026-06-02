#!/bin/bash
# Script: TAV + GMMLV combined RdRp phylogenetics.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 02/06/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Input directories
BLAST_DIR="${PROJECT_ROOT}/results/czid/BLAST"
PHYLO_OUT="${PROJECT_ROOT}/results/czid/phylogenetics/TAV_GMMLV"

# Viral protein database
VIRAL_DB="${PROJECT_ROOT}/data/databases/viral_refseq_prot"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
MAFFT_SIF="${CONTAINER_DIR}/mafft_7.221--0.sif"
TRIMAL_SIF="${CONTAINER_DIR}/trimal_1.5.1--h9948957_0.sif"
IQTREE_SIF="${CONTAINER_DIR}/iqtree_3.1.2--h8471819_0.sif"
CDHIT_SIF="${CONTAINER_DIR}/cd-hit_4.8.1--h5ca1c30_13.sif"
BLAST_SIF="${CONTAINER_DIR}/blast_2.17.0--h66d330f_0.sif"

THREADS=28
mkdir -p "${PHYLO_OUT}"

echo "=========================================="
echo "TAV + GMMLV: Combined RdRp Phylogenetics"
echo "=========================================="

####============================####
####   Collect RdRp Sequences    ####
####============================####

echo "Collecting RdRp sequences..."

> "${PHYLO_OUT}/all_rdrp.fasta"

# TAV RdRp (1 sequence)
if [ -f "${BLAST_DIR}/PM2469_S1_923325_split_orfs/PM2469_S1_923325_orf_6.fasta" ]; then
    awk -v label="TAV_PM2469" '/^>/{print ">"label; next} {print}' \
        "${BLAST_DIR}/PM2469_S1_923325_split_orfs/PM2469_S1_923325_orf_6.fasta" \
        >> "${PHYLO_OUT}/all_rdrp.fasta"
    echo "  Added TAV PM2469"
fi

# GMMLV RdRp (1 sequence — use the best one: PM2469 ORF6)
if [ -f "${BLAST_DIR}/PM2469_S1_923325_split_orfs/PM2469_S1_923325_orf_6.fasta" ]; then
    awk -v label="GMMLV_PM2469" '/^>/{print ">"label; next} {print}' \
        "${BLAST_DIR}/PM2469_S1_923325_split_orfs/PM2469_S1_923325_orf_6.fasta" \
        >> "${PHYLO_OUT}/all_rdrp.fasta"
    echo "  Added GMMLV PM2469"
fi

echo "RdRp sequences collected: $(grep -c '^>' ${PHYLO_OUT}/all_rdrp.fasta)"

####============================####
####   Get References (BLAST)    ####
####============================####

echo "Finding reference sequences..."

# Use the longest RdRp (TAV PM2469) to find references
apptainer exec "${BLAST_SIF}" blastp \
    -query "${BLAST_DIR}/PM2469_S1_923325_split_orfs/PM2469_S1_923325_orf_6.fasta" \
    -db "$VIRAL_DB" \
    -out "${PHYLO_OUT}/rdrp_blast.tsv" \
    -outfmt "6 sseqid" \
    -max_target_seqs 200

cut -f1 "${PHYLO_OUT}/rdrp_blast.tsv" | sort -u > "${PHYLO_OUT}/accessions.txt"

echo "Downloading reference sequences..."
> "${PHYLO_OUT}/references.fasta"
while read acc; do
    wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text" \
        >> "${PHYLO_OUT}/references.fasta"
    sleep 0.3
done < "${PHYLO_OUT}/accessions.txt"

echo "References downloaded: $(grep -c '^>' ${PHYLO_OUT}/references.fasta)"

####============================####
####      Combine + Rename      ####
####============================####

echo "Combining sequences..."
cat "${PHYLO_OUT}/all_rdrp.fasta" \
    "${PHYLO_OUT}/references.fasta" \
    > "${PHYLO_OUT}/combined_all.fasta"

# Rename reference headers to short format
awk '/^>/ {
    if ($0 ~ /^(TAV_|GMMLV_)/) {
        print
    } else {
        gsub(/^>/, "")
        split($1, acc, " ")
        printf ">ref_%s\n", acc[1]
    }
    next
}
{ print }' "${PHYLO_OUT}/combined_all.fasta" > "${PHYLO_OUT}/combined_renamed.fasta"

####============================####
####         Deduplicate        ####
####============================####

echo "Deduplicating sequences..."
apptainer exec "${CDHIT_SIF}" cd-hit \
    -i "${PHYLO_OUT}/combined_renamed.fasta" \
    -o "${PHYLO_OUT}/combined_dedup.fasta" \
    -c 0.85 \
    -T "${THREADS}"

# Ensure your 2 sequences are kept
cat "${PHYLO_OUT}/combined_dedup.fasta" \
    "${PHYLO_OUT}/all_rdrp.fasta" \
    > "${PHYLO_OUT}/combined_final.fasta"

# Remove any duplicate headers
awk '!seen[$0]++' "${PHYLO_OUT}/combined_final.fasta" > "${PHYLO_OUT}/combined_final_clean.fasta"
mv "${PHYLO_OUT}/combined_final_clean.fasta" "${PHYLO_OUT}/combined_final.fasta"

echo "Sequences after dedup: $(grep -c '^>' ${PHYLO_OUT}/combined_final.fasta)"

####============================####
####      Align, Trim, Tree     ####
####============================####

echo "Aligning with MAFFT..."
apptainer exec "${MAFFT_SIF}" mafft --auto --reorder \
    "${PHYLO_OUT}/combined_final.fasta" \
    > "${PHYLO_OUT}/combined_aligned.fasta"

echo "Trimming alignment..."
apptainer exec "${TRIMAL_SIF}" trimal \
    -in "${PHYLO_OUT}/combined_aligned.fasta" \
    -out "${PHYLO_OUT}/combined_trimmed.fasta" \
    -gt 0.1

echo "Building tree with IQ-TREE..."
apptainer exec "${IQTREE_SIF}" iqtree \
    -s "${PHYLO_OUT}/combined_trimmed.fasta" \
    -m MFP -bb 1000 -alrt 1000 -nt AUTO \
    -pre "${PHYLO_OUT}/combined_tree"

echo ""
echo "Tree: ${PHYLO_OUT}/combined_tree.treefile"
echo "Complete."