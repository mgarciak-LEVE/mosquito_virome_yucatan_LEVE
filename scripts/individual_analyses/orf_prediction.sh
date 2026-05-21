#!/bin/bash

# Script for ORF prediction from viral consensus genomes and statistics related to them.
# Ver. 1.0.0
# 20/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes"
ORF_DIR="${PROJECT_ROOT}/data/raw/czid_raw/orfs"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
EMBOSS_SIF="${CONTAINER_DIR}/samtools_1.20--h50ea8bc_1.sif" #### CHECKOUT IT OUT

# Virus names
VIRUS_NAMES=(
    "OSFV"
    "TAV"
    "GMMLV"
)

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

# Output directory for coverage results
stats_csv=${PROJECT_ROOT}/results/orf_predicition/orf_coverage_stats.csv
# Crear directorio si no existe
mkdir -p "$(dirname "$stats_csv")"

# CSV file for summary statistics
echo "sample,virus,orf_id,start,end,strand,frame,length_nt,length_aa" > "$stats_csv"


echo "=========================================="
echo "ORF prediction from viral consensus genomes"
echo "=========================================="
echo "Consensus directory: ${CONSENSUS_DIR}"
tg_send "Starting ORF prediction..."

####==================================####
####     PROCESS EACH VIRUS TYPE      ####
####==================================####

for i in "${!VIRUS_NAMES[@]}"; do

    VIRUS_NAME="${VIRUS_NAMES[$i]}"
    
    VIRUS_INPUT="${CONSENSUS_DIR}/${VIRUS_NAME}"
    VIRUS_OUTPUT="${ORF_DIR}/${VIRUS_NAME}"

    echo "=========================================="
    echo "Virus: ${VIRUS_NAME}"
    echo "Input:       ${VIRUS_INPUT}"
    echo "Output:      ${VIRUS_OUTPUT}"
    echo "=========================================="

    # Check if virus directory exists
    if [ ! -d "$VIRUS_INPUT" ]; then
        echo "WARNING: Directory not found for ${VIRUS_NAME}"
        tg_send "WARNING: No  directory for ${VIRUS_NAME}. Skipping."
        continue
    fi

    # Create output directory
    mkdir -p "${VIRUS_OUTPUT}"

####==================================####
####   PROCESS EACH SAMPLE DIRECTORY  ####
####==================================####

    for sample_dir in "${VIRUS_INPUT}"/*/; do

        sample_name=$(basename "$sample_dir")
        fasta_file="${sample_dir}/${sample_name}_consensus.fa"

        # Skip if BAM file doesn't exist
        if [ ! -f "$fasta_file" ]; then
            echo "WARNING: No fasta file for ${sample_name}. Skipping."
            continue
        fi

        # Find the original consensus file
        fasta_file=$(ls "${VIRUS_INPUT}/${sample_name}"* 2>/dev/null | head -1)

        if [ -z "$fasta_file" ] || [ ! -f "$fasta_file" ]; then
            echo "  WARNING: No FASTA found for ${sample_name}. Skipping."
            tg_send "WARNING: No FASTA for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        if [ ! -s "$fasta_file" ]; then
            echo "  WARNING: FASTA file is empty for ${sample_name}. Skipping."
            tg_send "WARNING: Empty FASTA file for ${VIRUS_NAME}/${sample_name}"
            continue
        fi

        echo "----------------------------------------"
        echo "Sample: ${sample_name}"
        echo "FASTA:    $(basename "$fasta_file")"
        echo "----------------------------------------"

        apptainer exec "${EMBOSS_SIF}" getorf -sequence "$fasta_file" -outseq "${VIRUS_OUTPUT}/${sample_name}_orfs.fasta" -minsize 300 -find 1

        ####============================####
        ####        ORF STATISTICS      ####
        ####============================####

        # ORF output header
        header=$(grep "^>" "${VIRUS_OUTPUT}/${sample_name}_orfs.fasta" | head -1)
        
        # ORF id 
        orf_id=$(grep "^>" "${VIRUS_OUTPUT}/${sample_name}_orfs.fasta" |  head -1 | awk '{print $1}' | rev | cut -d'_' -f1 | rev) # Gets the consensus number

        # Start position
        start=$(echo "$header" | grep -oP '\[\K\d+')

        # End position
        end=$(echo "$header" | grep -oP '\d+(?=\])')

        # Strand
        if echo "$header" | grep -q "REVERSE"; then
            strand="-"
        else
            strand="+"
        fi

        # Frame calculation
        frame=$(( (start - 1) % 3 + 1 ))

        # Length in nucleotides
        length_nt=$((end - start + 1))

        # Length in amino acids (approximate)
        length_aa=$((length_nt / 3))

        # CSV summary for all samples
        echo "${sample_name},${VIRUS_NAME},${orf_id},${start},${end},${strand},${frame},${length_nt},${length_aa}" >> "$stats_csv"
    
    done

done