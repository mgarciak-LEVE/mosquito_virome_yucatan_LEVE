#!/bin/bash
# Script: Ochlerotatus scapularis flavivirus (OSFV) genome completion.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 02/06/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes/OSFV"
READS_DIR="${PROJECT_ROOT}/data/raw/small_RNA"
OUTPUT_DIR="${PROJECT_ROOT}/results/czid/osfv_complete_genome"

# Conda environment activation
source "$(conda info --base)/etc/profile.d/conda.sh"

THREADS=14
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "OSFV: Genome Completion"
echo "=========================================="

ALL_CONSENSUS="${OUTPUT_DIR}/all_osfv_consensus.fasta"
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
    
    # Activate environment with bwa
    conda activate bwa_env
    
    echo "  Building bwa index..."
    bwa index "$local_consensus"
    
    echo "  Aligning reads..."
    bwa aln -t "${THREADS}" -n 2 -o 0 \
        "$local_consensus" "$clean_reads" \
        > "${OUTPUT_DIR}/${sample}.sai"
    
    bwa samse \
        "$local_consensus" "${OUTPUT_DIR}/${sample}.sai" "$clean_reads" \
        > "${OUTPUT_DIR}/${sample}.sam"
    
    conda deactivate
    
    # Convert to BAM
    echo "  Converting to BAM..."
    conda activate samtools_env
    samtools sort -@ "${THREADS}" \
        "${OUTPUT_DIR}/${sample}.sam" -o "${OUTPUT_DIR}/${sample}.bam"
    samtools index "${OUTPUT_DIR}/${sample}.bam"
    
    # Generate mpileup
    samtools mpileup -aa -A -d 0 \
        -f "$local_consensus" "${OUTPUT_DIR}/${sample}.bam" \
        > "${OUTPUT_DIR}/${sample}.mpileup"
    
    conda deactivate
    
    # Generate improved consensus
    echo "  Generating improved consensus..."
    conda activate ivar_env
    cat "${OUTPUT_DIR}/${sample}.mpileup" | ivar consensus \
        -p "${OUTPUT_DIR}/${sample}_improved" \
        -m 3 -t 0.5 -n N
    
    conda deactivate
    
    # Coverage stats
    conda activate samtools_env
    depth_file="${OUTPUT_DIR}/${sample}_depth.txt"
    samtools depth -a "${OUTPUT_DIR}/${sample}.bam" > "$depth_file"
    
    conda deactivate
    
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