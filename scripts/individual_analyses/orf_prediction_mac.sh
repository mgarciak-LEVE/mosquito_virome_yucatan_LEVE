#!/bin/bash

# Script for ORF prediction from viral consensus genomes and statistics related to them for mac environment.
# Ver. 1.0.0
# 22/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONSENSUS_DIR="${PROJECT_ROOT}/data/raw/czid_raw/consensus_genomes"
ORF_DIR="${PROJECT_ROOT}/data/raw/czid_raw/orfs"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
EMBOSS_SIF="${CONTAINER_DIR}/emboss_6.6.0--h0f19ade_14.sif"

# Virus names
VIRUS_NAMES=(
    "OSFV"
    "TAV"
    "GMMLV"
)

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"
# Conda activation
source "$(conda info --base)/etc/profile.d/conda.sh" 

# Output directory for coverage results
stats_csv=${PROJECT_ROOT}/results/orf_predicition/czid/orf_coverage_stats.csv
# Create directory for stats if it doesn't exist
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

     for fasta_file in "${VIRUS_INPUT}"/*_consensus.fa; do

        if [ ! -f "$fasta_file" ]; then
            continue
        fi

        if [ ! -s "$fasta_file" ]; then
            echo "  WARNING: FASTA file for ${VIRUS_NAME} is empty. Skipping."
            tg_send "WARNING: FASTA file for ${VIRUS_NAME} is empty. Skipping."
            continue
        fi

        sample_name=$(basename "$fasta_file" _consensus.fa)

        echo "----------------------------------------"
        echo "Sample: ${sample_name}"
        echo "FASTA:  $(basename "$fasta_file")"
        echo "----------------------------------------"

        conda activate emboss_env
        
        getorf -sequence "$fasta_file" -outseq "${VIRUS_OUTPUT}/${sample_name}_orfs.fasta" -minsize 300 -find 1

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
        tg_send "ORF prediction for ${sample_name} (${VIRUS_NAME}): ORF ${orf_id}, Start: ${start}, End: ${end}, Strand: ${strand}, Frame: ${frame}, Length (nt): ${length_nt}, Length (aa): ${length_aa}"

    done

done