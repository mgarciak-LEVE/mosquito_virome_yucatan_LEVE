#!/bin/bash
# Script: Guadeloupe mosquito mononega-like virus (GMMLV) genome completion.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 2.0.0
# 01/06/2026

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes/GMMLV"
READS_DIR="${PROJECT_ROOT}/data/raw/small_RNA"
OUTPUT_DIR="${PROJECT_ROOT}/results/czid/gmmlv_complete_genome"

CONTAINER_DIR="${PROJECT_ROOT}/containers"
BWA_SIF="${CONTAINER_DIR}/bwa_0.7.19--h577a1d6_1.sif"
SAMTOOLS_SIF="${CONTAINER_DIR}/samtools_1.20--h50ea8bc_1.sif"
IVAR_SIF="${CONTAINER_DIR}/ivar_1.4.4--h077b44d_0.sif"
EMBOSS_SIF="${CONTAINER_DIR}/emboss_6.6.0--h0f19ade_14.sif"
BLAST_SIF="${CONTAINER_DIR}/blast_2.17.0--h66d330f_0.sif"
VIRAL_DB="${PROJECT_ROOT}/data/databases/viral_refseq_prot"

THREADS=28
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "GMMLV: Genome Completion"
echo "=========================================="

ALL_CONSENSUS="${OUTPUT_DIR}/all_gmmlv_consensus.fasta"
> "$ALL_CONSENSUS"

# Use the best GMMLV consensus as reference for read mapping
REFERENCE=$(ls "${CONSENSUS_DIR}"/PM2469_S1_923325_consensus.fa | head -1)
echo "Reference: $(basename "$REFERENCE")"

# Sample-to-library mapping
declare -A SAMPLE_MAP
SAMPLE_MAP["PM2469"]="PM2469"
SAMPLE_MAP["PM2622"]="PM2622"
SAMPLE_MAP["PM3183"]="PM3183"
SAMPLE_MAP["PM3050"]="PM3050"
SAMPLE_MAP["PM2486"]="PM2486"
SAMPLE_MAP["PM3153"]="PM3153"
SAMPLE_MAP["PM2616"]="PM2616"
SAMPLE_MAP["PM2601"]="PM2601"
SAMPLE_MAP["PM2577"]="PM2577"
SAMPLE_MAP["PM2539"]="PM2539"
SAMPLE_MAP["PM2528"]="PM2528"
SAMPLE_MAP["PM2524"]="PM2524"
SAMPLE_MAP["PM2494"]="PM2494"


# Add more mappings as needed

for consensus in "${CONSENSUS_DIR}"/*_consensus.fa; do
    sample=$(basename "$consensus" _consensus.fa)
    base_id=$(echo "$sample" | grep -oP '^[A-Z]+\d+')
    
    echo "Processing: ${sample} (base ID: ${base_id})"
    
    # Check for mapped library
    mapped_id="${SAMPLE_MAP[$base_id]}"
    if [ -n "$mapped_id" ]; then
        reads=$(ls "${READS_DIR}/${mapped_id}"*.fa 2>/dev/null | head -1)
        echo "  Using mapped library: ${mapped_id}"
    else
        reads=$(ls "${READS_DIR}/${base_id}"*.fa 2>/dev/null | head -1)
    fi
    
    [ -z "$reads" ] && { echo "  No reads for ${base_id}, skipping."; continue; }
    echo "  Reads: $(basename "$reads")"
    
    clean_reads="${OUTPUT_DIR}/${sample}_clean_reads.fa"
    awk '/^>/ { print ">" NR; next } { print }' "$reads" > "$clean_reads"
    
    # Copy the REFERENCE consensus locally
    local_consensus="${OUTPUT_DIR}/${sample}_consensus.fa"
    cp "$REFERENCE" "$local_consensus"
    
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