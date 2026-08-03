#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to assemble unmapped reads from STAR and Bowtie2
# 03/08/2026
# Ver. 1.0.0 (farm-ready)

####==================================####
####           CONFIGURATION          ####
####==================================####

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directories
RNASPADES_INPUT_DIR="${PROJECT_SCRATCH}/results/assembly/rnaSPAdes"
METASPADES_INPUT_DIR="${PROJECT_SCRATCH}/results/assembly/metaSPAdes"
MEGAHIT_INPUT_DIR="${PROJECT_SCRATCH}/results/assembly/MEGAhit"


# Output directory
OUTPUT_DIR="${PROJECT_SCRATCH}/results/statistics"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Assembly Statistics Recovery Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Inputs:  ${RNASPADES_INPUT_DIR} , ${METASPADES_INPUT_DIR} , ${MEGAHIT_INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting assembly statistics recovery for ${PROJECT_NAME}" 2>/dev/null || true

# Create output directory
mkdir -p "${OUTPUT_DIR}/results/statistics"

echo "Generating assembly statistics with seqkit..."
tg_send "Generating assembly statistics with seqkit..." 2>/dev/null || true

#Creates CSV file
stats_file="${OUTDIR}/assembly/statistics/assembly_summary.csv"

#Specifies the file column headers
echo "Sample,Alignment,File,Num_Contigs,Total_Length,Max_Length,Min_Length,Avg_Length,N50,GC_Percent" > "$stats_file"

conda activate seqkit_env

# rnaSPAdes statistics
for sample_dir in "$OUTDIR/assembly/rnaSPAdes"/*/; do
    if [[ -d "$sample_dir" ]]; then # If its a directory...
        sample=$(basename "$sample_dir")
        
        # RNAspAdes typically produces these files:
        for contig_file in "$sample_dir/transcripts.fasta" \
                            "$sample_dir/soft_filtered_transcripts.fasta" \
                            "$sample_dir/hard_filtered_transcripts.fasta"; do
            if [[ -f "$contig_file" ]]; then # If the file exists...
                echo "Processing $(basename "$contig_file") for $sample - rnaSPAdes"
                
                # Get statistics with seqkit
                stats=$(seqkit stats "$contig_file" -T | awk 'NR==2')  # -T seqkit data output in tabular format. NR==2 means that it only gets the data row
                
                # Parse statistics
                num_contigs=$(echo "$stats" | awk '{print $4}')
                total_length=$(echo "$stats" | awk '{print $5}')
                min_length=$(echo "$stats" | awk '{print $6}')
                max_length=$(echo "$stats" | awk '{print $7}')
                avg_length=$(echo "$stats" | awk '{print $8}')

                # N50 calculation
                # seqkit fx2tab -n -l converts FASTA to table with names and lengths.
                n50=$(seqkit fx2tab -n -l "$contig_file" | awk '
                {len[NR]=$2; sum+=$2} 
                END {
                    half=sum/2;
                    for(i=1;i<=NR;i++) {
                        running_sum+=len[i];
                        if(running_sum>=half) {
                            print len[i];
                            exit;
                        }
                    }
                }')
                #len[NR]=$2; sum+=$2: Stores lengths in array, calculates total sum
                #half=sum/2: Calculates halfway point
                #running_sum+=len[i]: Adds lengths sequentially
                #if(running_sum>=half): Checks when we cross the halfway point
                #print len[i]: Outputs the N50 value


                # GC content
                #seqkit fx2tab -n -g gets GC percentage for each sequence
                gc_percent=$(seqkit fx2tab -n -g "$contig_file" | awk '
                {sum_gc+=$2; count++} 
                END {if(count>0) printf "%.2f", sum_gc/count; else print "0"}')

                filename=$(basename "$contig_file")
                echo "$sample,rnaSPAdes,$filename,$num_contigs,$total_length,$max_length,$min_length,$avg_length,$n50,$gc_percent" >> "$stats_file"
            fi
        done
    fi
done


# metaSPAdes statistics

for sample_dir in "$OUTDIR/assembly/metaSPAdes"/*/; do
    if [[ -d "$sample_dir" ]]; then
        sample=$(basename "$sample_dir")
        
        # metaSPAdes typically has these files
        for contig_file in "$sample_dir/contigs.fasta" "$sample_dir/scaffolds.fasta"; do
            if [[ -f "$contig_file" ]]; then
                file_desc=$(basename "$contig_file" .fasta)
                echo "Processing $file_desc for $sample - metaSPAdes"
                
                # Get statistics with seqkit
                stats=$(seqkit stats "$contig_file" -T | awk 'NR==2')
                
                # Parse statistics
                num_contigs=$(echo "$stats" | awk '{print $4}')
                total_length=$(echo "$stats" | awk '{print $5}')
                min_length=$(echo "$stats" | awk '{print $6}')
                max_length=$(echo "$stats" | awk '{print $7}')
                avg_length=$(echo "$stats" | awk '{print $8}')


                # N50 calculation
                n50=$(seqkit fx2tab -n -l "$contig_file" | awk '
                {len[NR]=$2; sum+=$2} 
                END {
                    half=sum/2;
                    for(i=1;i<=NR;i++) {
                        running_sum+=len[i];
                        if(running_sum>=half) {
                            print len[i];
                            exit;
                        }
                    }
                }')
                
                # GC content
                gc_percent=$(seqkit fx2tab -n -g "$contig_file" | awk '
                {sum_gc+=$2; count++} 
                END {if(count>0) printf "%.2f", sum_gc/count; else print "0"}')

                
                echo "$sample,metaSPAdes,$file_desc,$num_contigs,$total_length,$max_length,$min_length,$avg_length,$n50,$gc_percent" >> "$stats_file"
            fi
        done
    fi
done

echo "Statistics saved to: $stats_file"
echo "=== Assembly Statistics ==="
column -t -s, "$stats_file"

conda deactivate   
