#!/bin/bash
# Script: Tesano Aedes virus (TAV) genome completion.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 27/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes/TAV"
READS_DIR="${PROJECT_ROOT}/data/raw/small_RNA"
OUTPUT_DIR="${PROJECT_ROOT}/results/czid/tav_complete_genome"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
BWA_SIF="${CONTAINER_DIR}/bwa_0.7.19--h577a1d6_1.sif"
SAMTOOLS_SIF="${CONTAINER_DIR}/samtools_1.20--h50ea8bc_1.sif"
IVAR_SIF="${CONTAINER_DIR}/ivar_1.4.4--h077b44d_0.sif"

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

THREADS=28
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "TAV: Genome Completion"
echo "=========================================="
tg_send "Starting TAV genome completion..."

ALL_CONSENSUS="${OUTPUT_DIR}/all_tav_consensus.fasta"
> "$ALL_CONSENSUS"

for consensus in "${CONSENSUS_DIR}"/*_consensus.fa; do
    sample=$(basename "$consensus" _consensus.fa)
    base_id=$(echo "$sample" | grep -oP '^[A-Z]+\d+')
    
    echo "Processing: ${sample} (base ID: ${base_id})"
    
    reads=$(ls "${READS_DIR}/${base_id}"*.fa 2>/dev/null | head -1)
    [ -z "$reads" ] && { echo "  No reads for ${base_id}, skipping."; continue; }
    echo "  Reads: $(basename "$reads")"
    
    clean_reads="${OUTPUT_DIR}/${sample}_clean_reads.fa"
    awk '/^>/ { print ">" NR; next } { print }' "$reads" > "$clean_reads"
    
    # Copy consensus locally and index
    local_consensus="${OUTPUT_DIR}/${sample}_consensus.fa"
    cp "$consensus" "$local_consensus"
    
    echo "  Building bwa index..."
    apptainer exec "${BWA_SIF}" bwa index "$local_consensus"
    
    echo "  Aligning reads..."
    apptainer exec "${BWA_SIF}" bwa aln -t "${THREADS}" -n 2 -o 0 \
        "$local_consensus" "$clean_reads" \
        > "${OUTPUT_DIR}/${sample}.sai"
    
    apptainer exec "${BWA_SIF}" bwa samse \
        "$local_consensus" "${OUTPUT_DIR}/${sample}.sai" "$clean_reads" \
        > "${OUTPUT_DIR}/${sample}.sam"
    
    # Convert to BAM
    echo "  Converting to BAM..."
    apptainer exec "${SAMTOOLS_SIF}" samtools sort -@ "${THREADS}" \
        "${OUTPUT_DIR}/${sample}.sam" -o "${OUTPUT_DIR}/${sample}.bam"
    apptainer exec "${SAMTOOLS_SIF}" samtools index "${OUTPUT_DIR}/${sample}.bam"
    
    # Generate improved consensus
    echo "  Generating improved consensus..."
    apptainer exec "${SAMTOOLS_SIF}" samtools mpileup -aa -A -d 0 \
        -f "$local_consensus" "${OUTPUT_DIR}/${sample}.bam" \
        | apptainer exec "${IVAR_SIF}" ivar consensus \
            -p "${OUTPUT_DIR}/${sample}_improved" \
            -m 3 -t 0.5 -n N
    
    # Coverage stats
    depth_file="${OUTPUT_DIR}/${sample}_depth.txt"
    apptainer exec "${SAMTOOLS_SIF}" samtools depth -a "${OUTPUT_DIR}/${sample}.bam" > "$depth_file"
    
    avg_depth=$(awk '{sum+=$3; n++} END {printf "%.1f", sum/n}' "$depth_file")
    pct_cov=$(awk '$3>0{c++} END {printf "%.1f", c/NR*100}' "$depth_file")
    gaps=$(awk '$3==0{c++} END {print c}' "$depth_file")
    
    echo "  ${sample}: ${avg_depth}x depth, ${pct_cov}% coverage, ${gaps} gaps"
    
    cat "${OUTPUT_DIR}/${sample}_improved.fa" >> "$ALL_CONSENSUS"
    
    # Cleanup
    rm -f "${OUTPUT_DIR}/${sample}.sam" "${OUTPUT_DIR}/${sample}.sai" "$clean_reads"
done

echo ""
echo "Genome completion complete."
echo "Improved consensus: ${OUTPUT_DIR}/*_improved.fa"
echo "All consensus: ${ALL_CONSENSUS}"
tg_send "TAV genome completion complete."