#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to recover data associated with mapping statistics.
# 30/07/2026
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

# Input: Mapped statistics from STAR
INPUT_DIR="${PROJECT_SCRATCH}/results/"

# Output: alignment results
OUTPUT_DIR="${PROJECT_SCRATCH}/docs/mapping"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            
########==================================####

echo "========================================="
echo "  Mapping Statistics Recovery Module"
echo "  Project: ${PROJECT_NAME}"
echo "  Mapping logs: ${INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting mapping statistics recovery for ${PROJECT_NAME}" 2>/dev/null || true

# ASSEMBLY INFORMATION
echo "Generating mapping statistics..."
mkdir -p "${ALIGNED_DIR}/statistics"
    
#Creates CSV file
stats_file="${OUTDIR}/statistics/mapping_summary.csv"
#Specifies the file column headers
echo "Sample,Input_Reads,Average_Length,Uniquely_Mapped_Reads,Uniquely_Percent,Total_Multimapped,Total_Multimapped_Percent,Unmapped_Reads,Unmapped_Percent,Unmapped_Too_Short,Unmapped_Other" > "$stats_file"


for star_log in "$OUTDIR/aligned"/*Log.final.out; do
	if [[ -f "$star_log" ]]; then  # If the file exists... 
            filename=$(basename "$star_log" "_Log.final.out")
            
            # Parse sample and read type based on naming convention
            if [[ "$filename" == *"_paired_"* ]]; then
                sample=$(echo "$filename" | sed 's/_paired_//')
                read_type="paired"
            elif [[ "$filename" == *"_R1_unpaired_"* ]]; then
                sample=$(echo "$filename" | sed 's/_R1_unpaired_//')
                read_type="R1_unpaired"
            elif [[ "$filename" == *"_R2_unpaired_"* ]]; then
                sample=$(echo "$filename" | sed 's/_R2_unpaired_//')
                read_type="R2_unpaired"
            else
                sample="$filename"
                read_type="unknown"
            fi


                    # Parse statistics
                    input_reads=$(grep "Number of input reads" $star_log | awk '{print $6}')
                    avg_length=$(grep "Average input read length" $star_log | awk '{print $6}')
                    
                    #Uniquely mapped
                    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $star_log  | awk '{print $6}')
                    uniquely_percent=$(grep "Uniquely mapped reads %" $star_log  | awk '{print $6}' | sed 's/%//')
                    
                    #Multimapped
                    multi_mapped1=$(grep "Number of reads mapped to multiple loci" $star_log  | awk '{print $9}')
                    multi_mapped1_percent=$(grep "% of reads mapped to multiple loci" $star_log  | awk '{print $9}' | sed 's/%//')
                    multi_mapped2=$(grep "Number of reads mapped to too many loci" $star_log  | awk '{print $10}')
                    multi_mapped2_percent=$(grep "% of reads mapped to too many loci" $star_log  | awk '{print $10}' | sed 's/%//')
                    total_multi_mapped=$(($multi_mapped1 + $multi_mapped2))
                    total_multi_percent=$(echo "scale=2; $total_multi_mapped * 100 / $input_reads" | bc)

                    #Unmapped
                    unmapped_mismatches=$(grep "Number of reads unmapped: too many mismatches" $star_log  | awk '{print $9}')
                    unmapped_short=$(grep "Number of reads unmapped: too short" $star_log  | awk '{print $8}')
                    unmapped_other=$(grep "Number of reads unmapped: other" $star_log  | awk '{print $7}')
                    total_unmapped=$(($unmapped_mismatches + $unmapped_short + $unmapped_other))
                    unmapped_percent=$(echo "scale=2; $total_unmapped * 100 / $input_reads" | bc)

                    echo "$sample,$input_reads,$avg_length,$uniquely_mapped_reads,$uniquely_percent,$total_multi_mapped,$total_multi_percent,$total_unmapped,$unmapped_percent,$unmapped_short,$unmapped_other" >> "$stats_file"
                fi
    done
    
    echo "Statistics saved to: $stats_file"
    echo "=== Mapping Statistics ==="
    column -t -s, "$stats_file" 


