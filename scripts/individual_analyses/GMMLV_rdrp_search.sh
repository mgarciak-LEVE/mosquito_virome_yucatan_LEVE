#!/bin/bash
# Script: Guadeloupe mosquito mononega-like virus (GMMLV) RdRp search.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 28/05/2026

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes/GMMLV"
READS_DIR="${PROJECT_ROOT}/data/raw/small_RNA"
OUTPUT_DIR="${PROJECT_ROOT}/results/czid/gmmlv_rdrp_search"

CONTAINER_DIR="${PROJECT_ROOT}/containers"
MINIMAP2_SIF="${CONTAINER_DIR}/minimap2_2.30--h577a1d6_0.sif"
SAMTOOLS_SIF="${CONTAINER_DIR}/samtools_1.20--h50ea8bc_1.sif"
EMBOSS_SIF="${CONTAINER_DIR}/emboss_6.6.0--h0f19ade_14.sif"
BLAST_SIF="${CONTAINER_DIR}/blast_2.17.0--h66d330f_0.sif"
VIRAL_DB="${PROJECT_ROOT}/data/databases/viral_refseq_prot"

THREADS=28
mkdir -p "${OUTPUT_DIR}"

REFERENCE=$(ls "${CONSENSUS_DIR}"/PM2469_S1_923325_consensus.fa | head -1)
echo "Reference: $(basename "$REFERENCE")"

# Clean the reference header
CLEAN_REF="${OUTPUT_DIR}/reference_clean.fa"
awk 'NR==1{print ">PM2469_GMMLV"; next} {print}' "$REFERENCE" > "$CLEAN_REF"
REFERENCE="$CLEAN_REF"

# Predict ORFs
apptainer exec "${EMBOSS_SIF}" getorf -sequence "$REFERENCE" \
    -outseq "${OUTPUT_DIR}/gmmlv_orfs.fasta" -minsize 300 -find 1

# BLAST ORFs
apptainer exec "${BLAST_SIF}" blastp \
    -query "${OUTPUT_DIR}/gmmlv_orfs.fasta" \
    -db "$VIRAL_DB" \
    -out "${OUTPUT_DIR}/orfs_blast.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue stitle"

echo "Top RdRp hits:"
grep -i "polymerase\|RdRp\|rdrp" "${OUTPUT_DIR}/orfs_blast.tsv" | head -10

# Map reads from all samples
for reads in "${READS_DIR}"/*.fa; do
    [ ! -f "$reads" ] && continue
    sample=$(basename "$reads" .fa)
    echo "Mapping: ${sample}"

    apptainer exec "${MINIMAP2_SIF}" minimap2 \
        -ax sr \
        -t "${THREADS}" \
        "$REFERENCE" \
        "$reads" \
        > "${OUTPUT_DIR}/${sample}.sam"

    apptainer exec "${SAMTOOLS_SIF}" samtools sort -@ "${THREADS}" \
        "${OUTPUT_DIR}/${sample}.sam" -o "${OUTPUT_DIR}/${sample}.bam"
    
    if [ $? -ne 0 ]; then
        echo "  ERROR: sorting failed for ${sample}, trying with header fix..."
        # If sort fails, add a minimal SAM header
        apptainer exec "${SAMTOOLS_SIF}" samtools view -h -S "${OUTPUT_DIR}/${sample}.sam" \
            | awk 'NR==1{print "@HD\tVN:1.6\tSO:queryname"} {print}' \
            | apptainer exec "${SAMTOOLS_SIF}" samtools sort -@ "${THREADS}" - \
            -o "${OUTPUT_DIR}/${sample}.bam"
    fi

    apptainer exec "${SAMTOOLS_SIF}" samtools index "${OUTPUT_DIR}/${sample}.bam"

    apptainer exec "${SAMTOOLS_SIF}" samtools depth -a "${OUTPUT_DIR}/${sample}.bam" \
        > "${OUTPUT_DIR}/${sample}_depth.txt"

    avg=$(awk '{sum+=$3; n++} END {printf "%.1f", sum/n}' "${OUTPUT_DIR}/${sample}_depth.txt")
    cov=$(awk '$3>0{c++} END {printf "%.1f", c/NR*100}' "${OUTPUT_DIR}/${sample}_depth.txt")
    echo "  ${sample}: ${avg}x depth, ${cov}% coverage"

    rm -f "${OUTPUT_DIR}/${sample}.sam"
done

echo "Complete."
tg_send "GMMLV RdRp search complete."