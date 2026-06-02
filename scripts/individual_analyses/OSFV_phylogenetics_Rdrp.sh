#!/bin/bash
# Script: Ochlerotatus scapularis flavivirus (OSFV) RdRp phylogenetics.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.1.0
# 01/06/2026

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

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
MAFFT_SIF="${CONTAINER_DIR}/mafft_7.221--0.sif"
TRIMAL_SIF="${CONTAINER_DIR}/trimal_1.5.1--h9948957_0.sif"
IQTREE_SIF="${CONTAINER_DIR}/iqtree_3.1.2--h8471819_0.sif"
CDHIT_SIF="${CONTAINER_DIR}/cd-hit_4.8.1--h5ca1c30_13.sif"
BLAST_SIF="${CONTAINER_DIR}/blast_2.17.0--h66d330f_0.sif"
EMBOSS_SIF="${CONTAINER_DIR}/emboss_6.6.0--h0f19ade_14.sif"

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

THREADS=28
mkdir -p "${PHYLO_OUT}"

echo "=========================================="
echo "OSFV: RdRp Protein Phylogenetics"
echo "=========================================="
tg_send "Starting OSFV RdRp phylogenetics..."

# Check if improved consensus exists
if [ ! -f "$ALL_CONSENSUS" ]; then
    echo "ERROR: ${ALL_CONSENSUS} not found. Run genome completion first."
    tg_send "ERROR: OSFV consensus not found."
    exit 1
fi

####============================####
####        Predict ORFs        ####
####============================####

echo "Predicting ORFs from OSFV consensus..."
> "${PHYLO_OUT}/all_osfv_orfs.fasta"
for consensus in "${COMPLETION_DIR}"/*_improved.fa; do
    sample=$(basename "$consensus" _improved.fa)
    apptainer exec "${EMBOSS_SIF}" getorf \
        -sequence "$consensus" \
        -outseq "${PHYLO_OUT}/${sample}_orfs.fasta" \
        -minsize 300 -find 1
    cat "${PHYLO_OUT}/${sample}_orfs.fasta" >> "${PHYLO_OUT}/all_osfv_orfs.fasta"
done

####============================####
####     Find RdRp by BLAST     ####
####============================####

echo "BLASTing ORFs to find RdRp..."
apptainer exec "${BLAST_SIF}" blastp \
    -query "${PHYLO_OUT}/all_osfv_orfs.fasta" \
    -db "$VIRAL_DB" \
    -out "${PHYLO_OUT}/osfv_orfs_blast.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue stitle" \
    -max_target_seqs 10

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

# === INSECT-SPECIFIC FLAVIVIRUSES (cISF) ===
ISF_ACCESSIONS=(
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
)

# === DUAL-HOST INSECT-SPECIFIC FLAVIVIRUSES (dISF) ===
DISF_ACCESSIONS=(
    "YP_009253929.1"   # Kamiti River virus
    "YP_009483324.1"   # Palm Creek virus
    "YP_009551951.1"   # Nakiwogo virus
    "YP_006846328.2"   # Ilomantsi virus
)

# === MOSQUITO-BORNE FLAVIVIRUSES (MBF) ===
MBF_ACCESSIONS=(
    "NP_059433.1"      # Dengue virus 2
    "NP_041120.1"      # Yellow fever virus
    "NP_056776.1"      # West Nile virus
    "NP_619758.1"      # Zika virus
    "NP_872627.1"      # Japanese encephalitis virus
    "NP_051124.1"      # Murray Valley encephalitis virus
    "NP_041726.1"      # Saint Louis encephalitis virus
    "YP_009164031.1"   # Usutu virus (also in ISF, will be deduplicated)
)

# === TICK-BORNE FLAVIVIRUSES (TBF) ===
TBF_ACCESSIONS=(
    "NP_057949.1"      # Tick-borne encephalitis virus
    "NP_689391.1"      # Louping ill virus
    "NP_620044.1"      # Powassan virus
    "NP_658908.1"      # Omsk hemorrhagic fever virus
)

# === NO KNOWN VECTOR FLAVIVIRUSES (NKV) ===
NKV_ACCESSIONS=(
    "NP_891560.1"      # Rio Bravo virus
    "YP_009345019.1"   # Modoc virus
)

# === PESTIVIRUSES ===
PESTI_ACCESSIONS=(
    "NP_041913.1"      # Bovine viral diarrhea virus 1
    "NP_659601.1"      # Border disease virus
    "YP_009513190.1"   # Pestivirus K
)

# === HEPACIVIRUSES ===
HEPACI_ACCESSIONS=(
    "NP_671491.1"      # Hepatitis C virus NS5B
    "YP_001491523.1"   # GB virus B
)

# === PEGIVIRUSES ===
PEGI_ACCESSIONS=(
    "YP_009272544.1"   # Pegivirus A
    "YP_009256192.1"   # Pegivirus B
)

# Download all
for acc in "${ISF_ACCESSIONS[@]}" "${DISF_ACCESSIONS[@]}" "${MBF_ACCESSIONS[@]}" "${TBF_ACCESSIONS[@]}" "${NKV_ACCESSIONS[@]}" "${PESTI_ACCESSIONS[@]}" "${HEPACI_ACCESSIONS[@]}" "${PEGI_ACCESSIONS[@]}"; do
    wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text" \
        >> "${PHYLO_OUT}/references.fasta"
    sleep 0.5
done

# Remove exact duplicates from references
apptainer exec "${CDHIT_SIF}" cd-hit \
    -i "${PHYLO_OUT}/references.fasta" \
    -o "${PHYLO_OUT}/references_dedup.fasta" \
    -c 1.0 -T "${THREADS}"

mv "${PHYLO_OUT}/references_dedup.fasta" "${PHYLO_OUT}/references.fasta"

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
apptainer exec "${CDHIT_SIF}" cd-hit \
    -i "${PHYLO_OUT}/combined_renamed.fasta" \
    -o "${PHYLO_OUT}/combined_dedup.fasta" \
    -c 0.90 \
    -T "${THREADS}"

cat "${PHYLO_OUT}/combined_dedup.fasta" \
    "${PHYLO_OUT}/osfv_rdrp_queries.fasta" \
    > "${PHYLO_OUT}/combined_final.fasta"

echo "Sequences after dedup: $(grep -c '^>' ${PHYLO_OUT}/combined_dedup.fasta)"

####============================####
####      Align, Trim, Tree     ####
####============================####

echo "Aligning with MAFFT..."
apptainer exec "${MAFFT_SIF}" mafft --auto --reorder \
    "${PHYLO_OUT}/combined_final.fasta" \
    > "${PHYLO_OUT}/osfv_rdrp_aligned.fasta"

echo "Trimming alignment..."
apptainer exec "${TRIMAL_SIF}" trimal \
    -in "${PHYLO_OUT}/osfv_rdrp_aligned.fasta" \
    -out "${PHYLO_OUT}/osfv_rdrp_trimmed.fasta" \
    -gt 0.1

echo "Building tree with IQ-TREE..."
apptainer exec "${IQTREE_SIF}" iqtree \
    -s "${PHYLO_OUT}/osfv_rdrp_trimmed.fasta" \
    -m MFP -bb 1000 -alrt 1000 -nt AUTO \
    -pre "${PHYLO_OUT}/osfv_rdrp_tree"

echo ""
echo "Tree: ${PHYLO_OUT}/osfv_rdrp_tree.treefile"
echo "Complete."
tg_send "OSFV RdRp phylogeny complete."