#!/bin/bash
# Script: Map PM2528 small RNA reads to PM3153 consensus.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 09/06/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes/OSFV"
READS_DIR="${PROJECT_ROOT}/data/raw/small_RNA"
OUTPUT_DIR="${PROJECT_ROOT}/results/czid/osfv_complete_genome"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
BWA_SIF="${CONTAINER_DIR}/bwa_0.7.19--h577a1d6_1.sif"
SAMTOOLS_SIF="${CONTAINER_DIR}/samtools_1.20--h50ea8bc_1.sif"
IVAR_SIF="${CONTAINER_DIR}/ivar_1.4.4--h077b44d_0.sif"

THREADS=28
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "OSFV: PM2528 reads → PM3153 consensus"
echo "=========================================="

# Reference: PM3153 consensus (most complete, 10.2 kb)
REFERENCE="${CONSENSUS_DIR}/PM3153_S5_923338_consensus.fa"
echo "Reference: $(basename "$REFERENCE")"

# Reads: PM2528 small RNA
READS="${READS_DIR}/PM2528.clean.fa"

if [ ! -f "$READS" ]; then
    echo "ERROR: PM2528 reads not found: ${READS}"
    exit 1
fi
echo "Reads: $(basename "$READS")"

# Clean read headers
CLEAN_READS="${OUTPUT_DIR}/PM2528_clean.fa"
awk '/^>/ { print ">" NR; next } { print }' "$READS" > "$CLEAN_READS"

# Index reference
echo "Building bwa index..."
apptainer exec "${BWA_SIF}" bwa index "$REFERENCE"

# Align
echo "Aligning reads..."
apptainer exec "${BWA_SIF}" bwa aln -t "${THREADS}" -n 2 -o 0 \
    "$REFERENCE" "$CLEAN_READS" \
    > "${OUTPUT_DIR}/PM2528_on_PM3153.sai"

apptainer exec "${BWA_SIF}" bwa samse \
    "$REFERENCE" "${OUTPUT_DIR}/PM2528_on_PM3153.sai" "$CLEAN_READS" \
    > "${OUTPUT_DIR}/PM2528_on_PM3153.sam"

# Convert to BAM
echo "Converting to BAM..."
apptainer exec "${SAMTOOLS_SIF}" samtools sort -@ "${THREADS}" \
    "${OUTPUT_DIR}/PM2528_on_PM3153.sam" -o "${OUTPUT_DIR}/PM2528_on_PM3153.bam"
apptainer exec "${SAMTOOLS_SIF}" samtools index "${OUTPUT_DIR}/PM2528_on_PM3153.bam"

# Generate consensus
echo "Generating consensus..."
apptainer exec "${SAMTOOLS_SIF}" samtools mpileup -aa -A -d 0 \
    -f "$REFERENCE" "${OUTPUT_DIR}/PM2528_on_PM3153.bam" \
    | apptainer exec "${IVAR_SIF}" ivar consensus \
        -p "${OUTPUT_DIR}/PM2528_on_PM3153_improved" \
        -m 3 -t 0.5 -n N

# Coverage stats
DEPTH_FILE="${OUTPUT_DIR}/PM2528_on_PM3153_depth.txt"
apptainer exec "${SAMTOOLS_SIF}" samtools depth -a "${OUTPUT_DIR}/PM2528_on_PM3153.bam" > "$DEPTH_FILE"

AVG=$(awk '{sum+=$3; n++} END {printf "%.1f", sum/n}' "$DEPTH_FILE")
COV=$(awk '$3>0{c++} END {printf "%.1f", c/NR*100}' "$DEPTH_FILE")
GAPS=$(awk '$3==0{c++} END {print c}' "$DEPTH_FILE")
MAX=$(awk '$3>max{max=$3} END {print max+0}' "$DEPTH_FILE")

echo ""
echo "=========================================="
echo "Results: PM2528 reads on PM3153 consensus"
echo "=========================================="
echo "  Average depth:  ${AVG}x"
echo "  Max depth:      ${MAX}x"
echo "  Coverage:       ${COV}%"
echo "  Gaps:           ${GAPS}"
echo ""
echo "Improved consensus: ${OUTPUT_DIR}/PM2528_on_PM3153_improved.fa"

# Cleanup
rm -f "${OUTPUT_DIR}/PM2528_on_PM3153.sam" \
      "${OUTPUT_DIR}/PM2528_on_PM3153.sai" \
      "$CLEAN_READS"

echo "Complete."