#!/bin/bash
# Script: Guadeloupe mosquito mononega-like virus (GMMLV) RdRp phylogenetics.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 07/06/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Input directories
COMPLETION_DIR="${PROJECT_ROOT}/results/czid/gmmlv_complete_genome"
ALL_CONSENSUS="${COMPLETION_DIR}/all_gmmlv_consensus.fasta"
PHYLO_OUT="${PROJECT_ROOT}/results/czid/phylogenetics/GMMLV"

# Viral protein database
VIRAL_DB="${PROJECT_ROOT}/data/databases/viral_refseq_prot"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
MAFFT_SIF="${CONTAINER_DIR}/mafft_7.221--0.sif"
TRIMAL_SIF="${CONTAINER_DIR}/trimal_1.5.1--h9948957_0.sif"
IQTREE_SIF="${CONTAINER_DIR}/iqtree_3.1.2--h8471819_0.sif"
CDHIT_SIF="${CONTAINER_DIR}/cd-hit_4.8.1--h5ca1c30_13.sif"
BLAST_SIF="${CONTAINER_DIR}/blast_2.17.0--h66d330f_0.sif"
EMBOSS_SIF="${CONTAINER_DIR}/emboss_6.6.0--h0f19ade_14.sif"

THREADS=28
mkdir -p "${PHYLO_OUT}"

echo "=========================================="
echo "GMMLV: RdRp Protein Phylogenetics"
echo "=========================================="

# Check if improved consensus exists
if [ ! -f "$ALL_CONSENSUS" ]; then
    echo "ERROR: ${ALL_CONSENSUS} not found. Run genome completion first."
    exit 1
fi

####============================####
####        Predict ORFs        ####
####============================####

echo "Predicting ORFs from GMMLV consensus..."
> "${PHYLO_OUT}/all_gmmlv_orfs.fasta"
for consensus in "${COMPLETION_DIR}"/*_improved.fa; do
    sample=$(basename "$consensus" _improved.fa)
    apptainer exec "${EMBOSS_SIF}" getorf \
        -sequence "$consensus" \
        -outseq "${PHYLO_OUT}/${sample}_orfs.fasta" \
        -minsize 300 -find 1
    cat "${PHYLO_OUT}/${sample}_orfs.fasta" >> "${PHYLO_OUT}/all_gmmlv_orfs.fasta"
done

####============================####
####     Find RdRp by BLAST     ####
####============================####

echo "BLASTing ORFs to find RdRp..."
apptainer exec "${BLAST_SIF}" blastp \
    -query "${PHYLO_OUT}/all_gmmlv_orfs.fasta" \
    -db "$VIRAL_DB" \
    -out "${PHYLO_OUT}/gmmlv_orfs_blast.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue stitle" \
    -max_target_seqs 10

echo "Top RdRp candidates:"
grep -i "polymerase\|RdRp\|rdrp\|L" "${PHYLO_OUT}/gmmlv_orfs_blast.tsv" | \
    awk '$5 < 1e-10' | head -10

####============================####
####      Extract RdRp ORFs     ####
####============================####

echo "Extracting RdRp ORFs..."
> "${PHYLO_OUT}/gmmlv_rdrp_queries.fasta"

for consensus in "${COMPLETION_DIR}"/*_improved.fa; do
    sample=$(basename "$consensus" _improved.fa)
    
    best_orf=$(grep "${sample}" "${PHYLO_OUT}/gmmlv_orfs_blast.tsv" | \
        grep -i "polymerase\|RdRp\|rdrp\|L" | \
        sort -k5,5g | head -1 | cut -f1)
    
    if [ -n "$best_orf" ]; then
        awk -v id="$best_orf" '
            /^>/ { keep = ($0 ~ id) }
            { if (keep) print }
        ' "${PHYLO_OUT}/${sample}_orfs.fasta" >> "${PHYLO_OUT}/gmmlv_rdrp_queries.fasta"
        echo "  ${sample}: $best_orf"
    fi
done

echo "RdRp queries: $(grep -c '^>' ${PHYLO_OUT}/gmmlv_rdrp_queries.fasta)"

####============================####
####      Get References        ####
####============================####

echo "Downloading curated mononegavirus RdRp references..."

> "${PHYLO_OUT}/references.fasta"

# GMMLV is a mononega-like virus
GMMLV_REFERENCES=(
    "YP_010799271.1"   # Culex tritaeniorhynchus Anphevirus RdRp
    "YP_009388622.1"   # Culex mononega-like virus 2
    "YP_010784507.1"   # Aedes aegypti anphevirus
    "YP_009302387.1"   # Xincheng Mosquito Virus
    "QEM39171.1"       # Guadeloupe mosquito mononega-like virus (REFERENCE)
    "XGU11768.1"       # Ferox mosquito mononega-like virus
    "XXK85565.1"       # Toxorhynchites amboinensis xinmovirus 1
    "WPR17627.1"       # Hemiptera mononega-like virus (polymerase, partial)
    "AYN07258.1"       # Tick mononegavirus (polymerase, partial)
    "QNM37838.1"       # Insect mononegavirales virus 2 (polymerase)
    "WPR16622.1"       # Millipede mononega-like virus (polymerase, partial)
    "DBA09009.1"       # Chalcocoris rutilans mononega-like virus 2 (polymerase)
    "DAZ89733.1"   # Chelostoma florisomne mononega-like virus (polymerase, partial)
)

for acc in "${GMMLV_REFERENCES[@]}"; do
    wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text" \
        >> "${PHYLO_OUT}/references.fasta"
    sleep 0.5
done

echo "References downloaded: $(grep -c '^>' ${PHYLO_OUT}/references.fasta)"

####============================####
####      Combine + Rename     ####
####============================####

echo "Combining sequences..."
cat "${PHYLO_OUT}/gmmlv_rdrp_queries.fasta" \
    "${PHYLO_OUT}/references.fasta" \
    > "${PHYLO_OUT}/combined_all.fasta"

awk '/^>/ {
    if ($0 ~ /Consensus_PM/) {
        match($0, /Consensus_(PM[^_]+_[^_]+_[^_]+).*/)
        printf ">GMMLV_%s\n", substr($1, 1)
    } else {
        gsub(/^>/, "")
        acc = $1
        desc = substr($0, index($0, " ") + 1)
        gsub(/ /, "_", desc)
        gsub(/[\[\]]/, "", desc)
        if (length(desc) > 3) {
            printf ">%s_%s\n", acc, desc
        } else {
            printf ">%s\n", acc
        }
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

# Ensure GMMLV sequences are kept
cat "${PHYLO_OUT}/combined_dedup.fasta" \
    "${PHYLO_OUT}/gmmlv_rdrp_queries.fasta" \
    > "${PHYLO_OUT}/combined_final.fasta"

awk '!seen[$0]++' "${PHYLO_OUT}/combined_final.fasta" > "${PHYLO_OUT}/combined_final_clean.fasta"
mv "${PHYLO_OUT}/combined_final_clean.fasta" "${PHYLO_OUT}/combined_final.fasta"

echo "Sequences after dedup: $(grep -c '^>' ${PHYLO_OUT}/combined_final.fasta)"

####============================####
####      Align, Trim, Tree     ####
####============================####

echo "Aligning with MAFFT..."
apptainer exec "${MAFFT_SIF}" mafft --auto --reorder \
    "${PHYLO_OUT}/combined_final.fasta" \
    > "${PHYLO_OUT}/gmmlv_rdrp_aligned.fasta"

echo "Trimming alignment..."
apptainer exec "${TRIMAL_SIF}" trimal \
    -in "${PHYLO_OUT}/gmmlv_rdrp_aligned.fasta" \
    -out "${PHYLO_OUT}/gmmlv_rdrp_trimmed.fasta" \
    -gt 0.1

echo "Building tree with IQ-TREE..."
apptainer exec "${IQTREE_SIF}" iqtree \
    -s "${PHYLO_OUT}/gmmlv_rdrp_trimmed.fasta" \
    -m MFP -bb 1000 -alrt 1000 -nt AUTO \
    -pre "${PHYLO_OUT}/gmmlv_rdrp_tree"

echo ""
echo "Tree: ${PHYLO_OUT}/gmmlv_rdrp_tree.treefile"
echo "Complete."