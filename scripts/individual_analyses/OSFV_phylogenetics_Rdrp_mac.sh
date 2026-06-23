#!/bin/bash
# Script: Ochlerotatus scapularis flavivirus (OSFV) RdRp phylogenetics.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 02/06/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Input: improved consensus from the completion script
COMPLETION_DIR="${PROJECT_ROOT}/results/czid/osfv_complete_genome"
ALL_CONSENSUS="${COMPLETION_DIR}/all_osfv_consensus.fasta"
PHYLO_OUT="${PROJECT_ROOT}/results/czid/osfv_rdrp/phylogeny"

# Viral protein database
VIRAL_DB="${PROJECT_ROOT}/data/databases/viral_refseq_prot"

# Conda environment activation
source "$(conda info --base)/etc/profile.d/conda.sh"

THREADS=14
mkdir -p "${PHYLO_OUT}"

echo "=========================================="
echo "OSFV: RdRp Protein Phylogenetics"
echo "=========================================="

# Check if improved consensus exists
if [ ! -f "$ALL_CONSENSUS" ]; then
    echo "ERROR: ${ALL_CONSENSUS} not found. Run genome completion first."
    exit 1
fi

####============================####
####        Predict ORFs        ####
####============================####

echo "Predicting ORFs from OSFV consensus..."
conda activate emboss_env

> "${PHYLO_OUT}/all_osfv_orfs.fasta"
for consensus in "${COMPLETION_DIR}"/*_improved.fa; do
    sample=$(basename "$consensus" _improved.fa)
    getorf -sequence "$consensus" \
        -outseq "${PHYLO_OUT}/${sample}_orfs.fasta" \
        -minsize 300 -find 1
    cat "${PHYLO_OUT}/${sample}_orfs.fasta" >> "${PHYLO_OUT}/all_osfv_orfs.fasta"
done

conda deactivate

####============================####
####     Find RdRp by BLAST     ####
####============================####

echo "BLASTing ORFs to find RdRp..."
conda activate blast_env

blastp -query "${PHYLO_OUT}/all_osfv_orfs.fasta" \
    -db "$VIRAL_DB" \
    -out "${PHYLO_OUT}/osfv_orfs_blast.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue stitle" \
    -max_target_seqs 10

conda deactivate

echo "Top RdRp candidates:"
grep -i "polymerase\|NS5\|RdRp\|rdrp" "${PHYLO_OUT}/osfv_orfs_blast.tsv" | \
    awk '$5 < 1e-10' | head -10

####============================####
####      Extract RdRp ORFs     ####
####============================####

echo "Extracting RdRp ORFs..."
> "${PHYLO_OUT}/osfv_rdrp_queries.fasta"

for consensus in "${COMPLETION_DIR}"/*_improved.fa; do
    sample=$(basename "$consensus" _improved.fa)
    
    best_orf=$(grep "${sample}" "${PHYLO_OUT}/osfv_orfs_blast.tsv" | \
        grep -i "polymerase\|NS5\|RdRp\|rdrp" | \
        sort -k5,5g | head -1 | cut -f1)
    
    if [ -n "$best_orf" ]; then
        awk -v id="$best_orf" '
            /^>/ { keep = ($0 ~ id) }
            { if (keep) print }
        ' "${PHYLO_OUT}/${sample}_orfs.fasta" >> "${PHYLO_OUT}/osfv_rdrp_queries.fasta"
        echo "  ${sample}: $best_orf"
    fi
done

echo "RdRp queries: $(grep -c '^>' ${PHYLO_OUT}/osfv_rdrp_queries.fasta)"

####============================####
####      Get References        ####
####============================####

echo "Downloading curated flavivirus RdRp references..."

> "${PHYLO_OUT}/references.fasta"

# === VERIFIED FLAVIVIRUS ACCESSIONS ONLY ===
# All NP_ accessions were removed — they returned non-viral proteins.
# Only YP_ accessions and BCI56825.1 are verified viral proteins.

ALL_ACCESSIONS=(
    # Insect-specific flaviviruses
    "YP_009350102.1"   # Xishuangbanna aedes flavivirus
    "YP_009351861.1"   # Menghai flavivirus
    "YP_009345035.1"   # Shuangao insect-specific flavivirus
    "YP_009259488.1"   # Hanko flavivirus
    "YP_009164031.1"   # Aedes flavivirus
    "YP_009169331.1"   # Culex flavivirus
    "YP_009001464.1"   # Anopheles flavivirus
    "YP_009041466.1"   # Culex pipiens flavivirus
    "YP_009056847.1"   # Calbertado virus
    "YP_009388577.1"   # Hubei insect-specific flavivirus
    "YP_009552767.1"   # Hubei insect-specific flavivirus 2
    "BCI56825.1"       # Ochlerotatus scapularis flavivirus (REFERENCE)
    # Dual-host ISF
    "YP_009253929.1"   # Kamiti River virus
    "YP_009483324.1"   # Palm Creek virus
    "YP_009551951.1"   # Nakiwogo virus
    "YP_006846328.2"   # Ilomantsi virus
    # Mosquito-borne flaviviruses
    "YP_009164031.1"   # Usutu virus
    # No known vector
    "YP_009345019.1"   # Modoc virus
    # Pestiviruses
    "YP_009513190.1"   # Pestivirus K
    # Hepaciviruses
    "NP_671491.1"      # Hepatitis C virus NS5B (verified)
    "YP_001491523.1"   # GB virus B
    # Pegiviruses
    "YP_009272544.1"   # Pegivirus A
    "YP_009256192.1"   # Pegivirus B
)

for acc in "${ALL_ACCESSIONS[@]}"; do
    wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text" \
        >> "${PHYLO_OUT}/references.fasta"
    sleep 0.5
done

# Remove exact duplicates from references
conda activate blast_env
makeblastdb -in "${PHYLO_OUT}/references.fasta" -dbtype prot -out "${PHYLO_OUT}/refs_db" 2>/dev/null
conda deactivate

echo "References downloaded: $(grep -c '^>' ${PHYLO_OUT}/references.fasta)"

####============================####
####      Combine + Rename     ####
####============================####

echo "Combining sequences..."
cat "${PHYLO_OUT}/osfv_rdrp_queries.fasta" \
    "${PHYLO_OUT}/references.fasta" \
    > "${PHYLO_OUT}/combined_all.fasta"

# Rename to short headers
awk '/^>/ {
    if ($0 ~ /Consensus_PM/) {
        match($0, /Consensus_(PM[^_]+_[^_]+_[^_]+).*/)
        printf ">OSFV_%s\n", substr($1, 1)
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
conda activate blast_env
cd-hit -i "${PHYLO_OUT}/combined_renamed.fasta" \
    -o "${PHYLO_OUT}/combined_dedup.fasta" \
    -c 0.90 -T "${THREADS}"
conda deactivate

# Ensure your sequences are kept
cat "${PHYLO_OUT}/combined_dedup.fasta" \
    "${PHYLO_OUT}/osfv_rdrp_queries.fasta" \
    > "${PHYLO_OUT}/combined_final.fasta"

# Remove duplicate headers
awk '!seen[$0]++' "${PHYLO_OUT}/combined_final.fasta" > "${PHYLO_OUT}/combined_final_clean.fasta"
mv "${PHYLO_OUT}/combined_final_clean.fasta" "${PHYLO_OUT}/combined_final.fasta"

echo "Sequences after dedup: $(grep -c '^>' ${PHYLO_OUT}/combined_final.fasta)"

####============================####
####      Align, Trim, Tree     ####
####============================####

echo "Aligning with MAFFT..."
conda activate bioinformatics_base
mafft --auto --reorder --thread "${THREADS}" \
    "${PHYLO_OUT}/combined_final.fasta" \
    > "${PHYLO_OUT}/osfv_rdrp_aligned.fasta"

echo "Trimming alignment..."
trimal -in "${PHYLO_OUT}/osfv_rdrp_aligned.fasta" \
    -out "${PHYLO_OUT}/osfv_rdrp_trimmed.fasta" \
    -gt 0.1

echo "Building tree with IQ-TREE..."
iqtree -s "${PHYLO_OUT}/osfv_rdrp_trimmed.fasta" \
    -m MFP -bb 1000 -alrt 1000 -nt AUTO \
    -pre "${PHYLO_OUT}/osfv_rdrp_tree"
conda deactivate

echo ""
echo "Tree: ${PHYLO_OUT}/osfv_rdrp_tree.treefile"
echo "Complete."