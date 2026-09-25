#!/bin/bash

# Script for ORF BLAST analysis. 
# Ver. 1.1.0
# 20/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ORF_DIR="${PROJECT_ROOT}/data/raw/czid_raw/orfs"
BLAST_OUTPUT_DIR="${PROJECT_ROOT}/results/czid/BLAST"

# OSFV reference sequence for BLAST database
OSFV_DB_DIR="${PROJECT_ROOT}/data/databases"
OSFV_REFERENCES="${OSFV_DB_DIR}/ncbi_dataset/data/genomic.fna"
OSFV_DB="${OSFV_DB_DIR}/osfv_references_db"

# Container paths
CONTAINER_DIR="${PROJECT_ROOT}/containers"
BLAST_SIF="${CONTAINER_DIR}/blast_2.17.0--h66d330f_0.sif"

# Virus names
VIRUS_NAMES=(
    "OSFV"
    "TAV"
    "GMMLV"
)

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

# Output directory for coverage results
stats_csv=${BLAST_OUTPUT_DIR}/orf_BLAST_stats.csv
# Create directory for stats if it doesn't exist
mkdir -p "$(dirname "$stats_csv")"

# CSV file for summary statistics
echo "sample,virus,orf_id,qseqid,sseqid,pident,length,qlen,slen,evalue,stitle" > "$stats_csv"

echo "=========================================="
echo "ORF BLAST for viral consensus genomes"
echo "=========================================="
echo "ORFdirectory: ${ORF_DIR}"
tg_send "Starting ORF BLAST..."

####==================================####
####     PROCESS EACH VIRUS TYPE      ####
####==================================####

for i in "${!VIRUS_NAMES[@]}"; do

    VIRUS_NAME="${VIRUS_NAMES[$i]}"
    
    VIRUS_INPUT="${ORF_DIR}/${VIRUS_NAME}"
    VIRUS_OUTPUT="${BLAST_OUTPUT_DIR}"

    echo "=========================================="
    echo "Virus: ${VIRUS_NAME}"
    echo "Input:       ${VIRUS_INPUT}"
    echo "Output:      ${BLAST_OUTPUT_DIR}"
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

     for orf_file in "${VIRUS_INPUT}"/*_orfs.fasta; do

        if [ ! -f "$orf_file" ]; then
            continue
        fi

        if [ ! -s "$orf_file" ]; then
            echo "  WARNING: ORF FASTA file for ${VIRUS_NAME} is empty. Skipping."
            tg_send "WARNING: ORF FASTA file for ${VIRUS_NAME} is empty. Skipping."
            continue
        fi

        sample_name=$(basename "$orf_file" _orfs.fasta)

        echo "----------------------------------------"
        echo "Sample: ${sample_name}"
        echo "FASTA:  $(basename "$orf_file")"
        echo "----------------------------------------"

        ####============================####
        ####    SPLIT MULTI-FASTA ORFS  ####
        ####============================####

        orf_split_dir="${VIRUS_OUTPUT}/${sample_name}_split_orfs"
        mkdir -p "$orf_split_dir"

        # Split multi-FASTA into individual ORF files
        awk -v virus="$VIRUS_NAME" -v sample="$sample_name" '/^>/ {
            match($0, /_([0-9]+) /, arr)
            orf_num = arr[1]
            # Create distinct header: >OSFV_PM2528_S2_923330_orf1
            new_header = sprintf(">%s_%s_orf%s", virus, sample, orf_num)
            print new_header
            next
        }
        {
            print
        }'  "$orf_file" > "${orf_split_dir}/${sample_name}_orfs.fasta"

        echo "  ORFs split into: $(basename "$orf_split_dir")"

        ####============================####
        ####   BLAST EACH INDIVIDUAL ORF ####
        ####============================####

        if [ "$VIRUS_NAME" = "OSFV" ]; then
            ####  OSFV: Nucleotide BLAST against custom DB  ####

            # Build the database once (only if not already built)
            if [ ! -f "${OSFV_DB}.nsq" ] && [ ! -f "${OSFV_DB}.nin" ]; then
                echo "  Building OSFV nucleotide database..."
                apptainer exec "${BLAST_SIF}" makeblastdb \
                    -in "$OSFV_REFERENCES" \
                    -dbtype nucl \
                    -parse_seqids \
                    -title "Flavivirus References" \
                    -out "$OSFV_DB"
            fi

            for single_orf in "${orf_split_dir}"/*_orf_*.fasta; do
                if [ ! -f "$single_orf" ]; then
                    continue
                fi
                orf_num=$(basename "$single_orf" .fasta | rev | cut -d'_' -f1 | rev)
                echo "  BLASTing ORF ${orf_num} (nucleotide, local DB)..."

                blast_out="${VIRUS_OUTPUT}/${sample_name}_orf${orf_num}_blastn.tsv"

                apptainer exec "${BLAST_SIF}" tblastn \
                    -query "$single_orf" \
                    -db "$OSFV_DB" \
                    -out "$blast_out" \
                    -outfmt "6 qseqid sseqid pident length qlen slen evalue stitle"

                if [ $? -eq 0 ]; then
                    # Parse top hit for CSV for OSFV
                    top_hit=$(head -1 "$blast_out" 2>/dev/null)

                    if [ -n "$top_hit" ]; then
                        qseqid=$(echo "$top_hit" | cut -f1)
                        sseqid=$(echo "$top_hit" | cut -f2)
                        pident=$(echo "$top_hit" | cut -f3)
                        length=$(echo "$top_hit" | cut -f4)
                        qlen=$(echo "$top_hit" | cut -f5)
                        slen=$(echo "$top_hit" | cut -f6)
                        evalue=$(echo "$top_hit" | cut -f7)
                        stitle=$(echo "$top_hit" | cut -f8-)
                    else
                        qseqid=""; sseqid=""; pident=""; length=""; qlen=""
                        slen=""; evalue=""; stitle="No hits found"
                    fi

                    echo "${sample_name},${VIRUS_NAME},${orf_num},${qseqid},${sseqid},${pident},${length},${qlen},${slen},${evalue},${stitle}" >> "$stats_csv"
                else
                    echo "    ERROR: BLAST failed for ORF ${orf_num}"
                    tg_send "ERROR: BLAST failed for ${VIRUS_NAME}/${sample_name} ORF ${orf_num}"
                fi
            done

        else
            ####  TAV & GMMLV: Protein BLAST remote against nr  ####

            for single_orf in "${orf_split_dir}"/*_orf_*.fasta; do
                if [ ! -f "$single_orf" ]; then
                    continue
                fi
                orf_num=$(basename "$single_orf" .fasta | rev | cut -d'_' -f1 | rev)
                echo "  BLASTing ORF ${orf_num} (protein, remote)..."

                blast_out="${VIRUS_OUTPUT}/${sample_name}_orf${orf_num}_blastp.tsv"

                apptainer exec "${BLAST_SIF}" blastp \
                    -query "$single_orf" \
                    -db "${PROJECT_ROOT}/data/databases/viral_refseq_prot" \
                    -out "$blast_out" \
                    -outfmt "6 qseqid sseqid pident length qlen slen evalue stitle"

                if [ $? -eq 0 ]; then
                    
                    # Parse top hit for CSV for TAV and GMMLV
                    top_hit=$(head -1 "$blast_out" 2>/dev/null)

                    if [ -n "$top_hit" ]; then
                        qseqid=$(echo "$top_hit" | cut -f1)
                        sseqid=$(echo "$top_hit" | cut -f2)
                        pident=$(echo "$top_hit" | cut -f3)
                        length=$(echo "$top_hit" | cut -f4)
                        qlen=$(echo "$top_hit" | cut -f5)
                        slen=$(echo "$top_hit" | cut -f6)
                        evalue=$(echo "$top_hit" | cut -f7)
                        stitle=$(echo "$top_hit" | cut -f8-)
                    else
                        qseqid=""; sseqid=""; pident=""; length=""; qlen=""
                        slen=""; evalue=""; stitle="No hits found"
                    fi

                    echo "${sample_name},${VIRUS_NAME},${orf_num},${qseqid},${sseqid},${pident},${length},${qlen},${slen},${evalue},${stitle}" >> "$stats_csv"
                else
                    echo "    ERROR: BLAST failed for ORF ${orf_num}"
                    tg_send "ERROR: BLAST failed for ${VIRUS_NAME}/${sample_name} ORF ${orf_num}"
                fi
            done

        fi

    done

done

echo "=========================================="
echo "ORF BLAST Complete"
echo "=========================================="
tg_send "ORF BLAST complete."
