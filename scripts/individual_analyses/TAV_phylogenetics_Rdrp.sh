#!/bin/bash
# Script: Tesano Aedes virus (TAV) RdRp phylogenetics.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 07/06/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Input directories
COMPLETION_DIR="${PROJECT_ROOT}/results/czid/tav_complete_genome"
ALL_CONSENSUS="${COMPLETION_DIR}/all_tav_consensus.fasta"
PHYLO_OUT="${PROJECT_ROOT}/results/czid/phylogenetics/TAV"

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
echo "TAV: RdRp Protein Phylogenetics"
echo "=========================================="

# Check if improved consensus exists
if [ ! -f "$ALL_CONSENSUS" ]; then
    echo "ERROR: ${ALL_CONSENSUS} not found. Run genome completion first."
    exit 1
fi

####============================####
####        Predict ORFs        ####
####============================####

echo "Predicting ORFs from TAV consensus..."
> "${PHYLO_OUT}/all_tav_orfs.fasta"
for consensus in "${COMPLETION_DIR}"/*_improved.fa; do
    sample=$(basename "$consensus" _improved.fa)
    apptainer exec "${EMBOSS_SIF}" getorf \
        -sequence "$consensus" \
        -outseq "${PHYLO_OUT}/${sample}_orfs.fasta" \
        -minsize 300 -find 1
    cat "${PHYLO_OUT}/${sample}_orfs.fasta" >> "${PHYLO_OUT}/all_tav_orfs.fasta"
done

####============================####
####     Find RdRp by BLAST     ####
####============================####

echo "BLASTing ORFs to find RdRp..."
apptainer exec "${BLAST_SIF}" blastp \
    -query "${PHYLO_OUT}/all_tav_orfs.fasta" \
    -db "$VIRAL_DB" \
    -out "${PHYLO_OUT}/tav_orfs_blast.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue stitle" \
    -max_target_seqs 10

echo "Top RdRp candidates:"
grep -i "polymerase\|RdRp\|rdrp\|L\|NS5" "${PHYLO_OUT}/tav_orfs_blast.tsv" | \
    awk '$5 < 1e-10' | head -10

####============================####
####      Extract RdRp ORFs     ####
####============================####

echo "Extracting RdRp ORFs..."
> "${PHYLO_OUT}/tav_rdrp_queries.fasta"

for consensus in "${COMPLETION_DIR}"/*_improved.fa; do
    sample=$(basename "$consensus" _improved.fa)
    
    best_orf=$(grep "${sample}" "${PHYLO_OUT}/tav_orfs_blast.tsv" | \
        grep -i "polymerase\|RdRp\|rdrp\|L" | \
        sort -k5,5g | head -1 | cut -f1)
    
    if [ -n "$best_orf" ]; then
        awk -v id="$best_orf" '
            /^>/ { keep = ($0 ~ id) }
            { if (keep) print }
        ' "${PHYLO_OUT}/${sample}_orfs.fasta" >> "${PHYLO_OUT}/tav_rdrp_queries.fasta"
        echo "  ${sample}: $best_orf"
    fi
done

echo "RdRp queries: $(grep -c '^>' ${PHYLO_OUT}/tav_rdrp_queries.fasta)"

####============================####
####      Get References        ####
####============================####

echo "Downloading curated mononegavirus RdRp references..."

> "${PHYLO_OUT}/references.fasta"

# TAV is a mononegavirus (related to Anphevirus, Nyavirus, etc.)
TAV_REFERENCES=(
    "YP_010799271.1"   # Culex tritaeniorhynchus Anphevirus RdRp
    "YP_009342321.1"   # Wuhan insect virus 13 (Hypothetical protein)
    "YP_010784507.1"   # Aedes anphevirus (Rdrp)
    "YP_009302387.1"   # Xincheng Mosquito Virus (Rdrp)
    "YP_009388622.1"   # Culex mononega-like virus 2 (Rdrp)
    "QPF16707.1"       # Aedes aegypti To virus 1 (putative RdRp)
    "WNO14524.1"       # Aedes binegev-like virus 2 (Rdrp)
    "XGU11772.1"       # Albipes mosquito Gordis-like virus (Rdrp)
    "QOW03291.1"       # Kisumu mosquito virus (Rdrp)
    "XUU16279.1"       # Kwale mosquito virus (Rdrp)
    "XNP01474.1"       # Zheijang mosquito virus (hypothetical protein)
    "UNZ11823.1"       # Talaya insect virus 1 (putative Rdrp)
    "UNZ11794.1"       # Icha Creek insect virus (putative polyprotein)
    "QNQ73517.1"       # Insect tombusbipa-like virus 1 (Rdrp)
    "UNZ11827.1"       # Uzakla insect virus (putative polyrpotein)
    "BBN21000.1"       # Tesano Aedes virus (REFERENCE)
)

for acc in "${TAV_REFERENCES[@]}"; do
    wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text" \
        >> "${PHYLO_OUT}/references.fasta"
    sleep 0.5
done

echo "References downloaded: $(grep -c '^>' ${PHYLO_OUT}/references.fasta)"

####============================####
####      Combine + Rename     ####
####============================####

echo "Combining sequences..."
cat "${PHYLO_OUT}/tav_rdrp_queries.fasta" \
    "${PHYLO_OUT}/references.fasta" \
    > "${PHYLO_OUT}/combined_all.fasta"

awk '/^>/ {
    if ($0 ~ /Consensus_PM/) {
        match($0, /Consensus_(PM[^_]+_[^_]+_[^_]+).*/)
        printf ">TAV_%s\n", substr($1, 1)
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

# Ensure TAV sequences are kept
cat "${PHYLO_OUT}/combined_dedup.fasta" \
    "${PHYLO_OUT}/tav_rdrp_queries.fasta" \
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
    > "${PHYLO_OUT}/tav_rdrp_aligned.fasta"

echo "Trimming alignment..."
apptainer exec "${TRIMAL_SIF}" trimal \
    -in "${PHYLO_OUT}/tav_rdrp_aligned.fasta" \
    -out "${PHYLO_OUT}/tav_rdrp_trimmed.fasta" \
    -gt 0.1

echo "Building tree with IQ-TREE..."
apptainer exec "${IQTREE_SIF}" iqtree \
    -s "${PHYLO_OUT}/tav_rdrp_trimmed.fasta" \
    -m MFP -bb 1000 -alrt 1000 -nt AUTO \
    -pre "${PHYLO_OUT}/tav_rdrp_tree"

echo ""
echo "Tree: ${PHYLO_OUT}/tav_rdrp_tree.treefile"
echo "Complete."