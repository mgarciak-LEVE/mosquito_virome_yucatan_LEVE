#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to recover data associated with mapping statistics for STAR and Bowtie2.
# 30/07/2026
# Ver. 1.1.0 (farm-ready)

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
STAR_INPUT_DIR="${PROJECT_SCRATCH}/results/aligned/STAR_alignment"

LSF_LOG_DATE="04_08_2026"  # Change this to your date if needed
BOWTIE_LSF_LOG_DIR="${HOME}/lsf_logs/${LSF_LOG_DATE}"

# Output directory
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
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Mapping Statistics Recovery Module"
echo "  Project: ${PROJECT_NAME}"
echo "  STAR logs: ${STAR_INPUT_DIR}"
echo "  Bowtie2 logs: ${BOWTIE_INPUT_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting mapping statistics recovery for ${PROJECT_NAME}" 2>/dev/null || true

# Create output directory
mkdir -p "${OUTPUT_DIR}/statistics"

####=====================================####
####     STAR MAPPING STATISTICS         ####
####=====================================####

echo ""
echo "=== STAR Mapping Statistics ==="

STAR_STATS_FILE="${OUTPUT_DIR}/statistics/star_mapping_summary.csv"
echo "Sample,Read_Type,Input_Reads,Average_Length,Uniquely_Mapped_Reads,Uniquely_Percent,Total_Multimapped,Total_Multimapped_Percent,Unmapped_Reads,Unmapped_Percent,Unmapped_Too_Short,Unmapped_Other" > "$STAR_STATS_FILE"

for star_log in "${STAR_INPUT_DIR}"/*_Log.final.out; do
    if [[ -f "$star_log" ]]; then
        filename=$(basename "$star_log" "_Log.final.out")
        
        # Parse sample and read type
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
        input_reads=$(grep "Number of input reads" "$star_log" | awk '{print $6}' 2>/dev/null || echo "0")
        avg_length=$(grep "Average input read length" "$star_log" | awk '{print $6}' 2>/dev/null || echo "0")
        
        uniquely_mapped_reads=$(grep "Uniquely mapped reads number" "$star_log" | awk '{print $6}' 2>/dev/null || echo "0")
        uniquely_percent=$(grep "Uniquely mapped reads %" "$star_log" | awk '{print $6}' | sed 's/%//' 2>/dev/null || echo "0")
        
        multi_mapped1=$(grep "Number of reads mapped to multiple loci" "$star_log" | awk '{print $9}' 2>/dev/null || echo "0")
        multi_mapped1_percent=$(grep "% of reads mapped to multiple loci" "$star_log" | awk '{print $9}' | sed 's/%//' 2>/dev/null || echo "0")
        multi_mapped2=$(grep "Number of reads mapped to too many loci" "$star_log" | awk '{print $10}' 2>/dev/null || echo "0")
        multi_mapped2_percent=$(grep "% of reads mapped to too many loci" "$star_log" | awk '{print $10}' | sed 's/%//' 2>/dev/null || echo "0")
        total_multi_mapped=$((multi_mapped1 + multi_mapped2))
        
        if [[ "$input_reads" -gt 0 ]]; then
            total_multi_percent=$(echo "scale=2; $total_multi_mapped * 100 / $input_reads" | bc 2>/dev/null || echo "0")
        else
            total_multi_percent="0"
        fi

        unmapped_mismatches=$(grep "Number of reads unmapped: too many mismatches" "$star_log" | awk '{print $9}' 2>/dev/null || echo "0")
        unmapped_short=$(grep "Number of reads unmapped: too short" "$star_log" | awk '{print $8}' 2>/dev/null || echo "0")
        unmapped_other=$(grep "Number of reads unmapped: other" "$star_log" | awk '{print $7}' 2>/dev/null || echo "0")
        total_unmapped=$((unmapped_mismatches + unmapped_short + unmapped_other))
        
        if [[ "$input_reads" -gt 0 ]]; then
            unmapped_percent=$(echo "scale=2; $total_unmapped * 100 / $input_reads" | bc 2>/dev/null || echo "0")
        else
            unmapped_percent="0"
        fi

        echo "$sample,$read_type,$input_reads,$avg_length,$uniquely_mapped_reads,$uniquely_percent,$total_multi_mapped,$total_multi_percent,$total_unmapped,$unmapped_percent,$unmapped_short,$unmapped_other" >> "$STAR_STATS_FILE"
    fi
done

echo "STAR statistics saved to: $STAR_STATS_FILE"

####=====================================####
####     BOWTIE2 MAPPING STATISTICS      ####
####=====================================####

echo ""
echo "=== Bowtie2 Mapping Statistics ==="

BOWTIE_STATS_FILE="${OUTPUT_DIR}/statistics/bowtie2_mapping_summary.csv"
echo "Sample,Read_Type,Total_Reads,Aligned_Concordantly_0,Aligned_Concordantly_1,Aligned_Concordantly_1_Percent,Aligned_Concordantly_GT1,Overall_Alignment_Rate" > "$BOWTIE_STATS_FILE"

for bowtie_log in "${BOWTIE_LSF_LOG_DIR}"/mapping_*_*.out; do
    if [[ ! -f "$bowtie_log" ]]; then
        continue
    fi
    
    # Extract sample name
    sample=$(grep "Sample:" "$bowtie_log" | head -1 | sed 's/.*Sample: //' || echo "unknown")
    
    # Skip if no sample found
    if [[ "$sample" == "unknown" ]]; then
        continue
    fi
    
    echo "Processing: $bowtie_log (Sample: $sample)"
    
    # Parse each stats block using a state machine approach
    total_reads=""
    aligned_0=""
    aligned_1=""
    aligned_gt1=""
    overall_rate=""
    read_type=""
    stats_block_started=0
    
    while IFS= read -r line; do
        # Detect start of a stats block
        if echo "$line" | grep -q "reads; of these:"; then
            # If we have data from previous block, save it
            if [[ -n "$total_reads" ]] && [[ -n "$read_type" ]]; then
                echo "$sample,$read_type,$total_reads,$aligned_0,$aligned_1,$aligned_gt1,$overall_rate" >> "$BOWTIE_STATS_FILE"
            fi
            
            # Reset for new block
            total_reads=""
            aligned_0=""
            aligned_1=""
            aligned_gt1=""
            overall_rate=""
            read_type=""
            stats_block_started=1
            
            # Get total reads
            total_reads=$(echo "$line" | awk '{print $1}' | sed 's/,//g')
        fi
        
        # Detect read type
        if echo "$line" | grep -q "were paired"; then
            read_type="paired"
        elif echo "$line" | grep -q "were unpaired"; then
            # Determine if R1 or R2 unpaired by looking at the surrounding context
            if echo "$line" | grep -q "R1_unpaired"; then
                read_type="R1_unpaired"
            elif echo "$line" | grep -q "R2_unpaired"; then
                read_type="R2_unpaired"
            else
                read_type="unpaired"
            fi
        fi
        
        # Parse alignment statistics
        if echo "$line" | grep -q "aligned 0 times"; then
            aligned_0=$(echo "$line" | awk '{print $1}' | sed 's/,//g')
        fi
        if echo "$line" | grep -q "aligned exactly 1 time"; then
            aligned_1=$(echo "$line" | awk '{print $1}' | sed 's/,//g')
        fi
        if echo "$line" | grep -q "aligned >1 times"; then
            aligned_gt1=$(echo "$line" | awk '{print $1}' | sed 's/,//g')
        fi
        if echo "$line" | grep -q "overall alignment rate"; then
            overall_rate=$(echo "$line" | awk '{print $1}' | sed 's/%//')
        fi
        
    done < "$bowtie_log"
    
    # Save the last block
    if [[ -n "$total_reads" ]] && [[ -n "$read_type" ]]; then
        echo "$sample,$read_type,$total_reads,$aligned_0,$aligned_1,$aligned_gt1,$overall_rate" >> "$BOWTIE_STATS_FILE"
    fi
    
done

echo "Bowtie2 statistics saved to: $BOWTIE_STATS_FILE"

####================================####
####        DISPLAY SUMMARY          ####
####================================####

echo ""
echo "========================================="
echo "  Mapping Statistics Summary"
echo "========================================="

if [[ -f "$STAR_STATS_FILE" ]]; then
    echo ""
    echo "=== STAR Mapping ==="
    column -t -s, "$STAR_STATS_FILE" 2>/dev/null || cat "$STAR_STATS_FILE"
fi

if [[ -f "$BOWTIE_STATS_FILE" ]]; then
    echo ""
    echo "=== Bowtie2 Mapping ==="
    column -t -s, "$BOWTIE_STATS_FILE" 2>/dev/null || cat "$BOWTIE_STATS_FILE"
fi

echo ""
echo "========================================="

####================================####
####     COPY RESULTS TO NFS         ####
####================================####

echo ""
echo "Copying results to permanent storage..."

PERMANENT_RESULTS="/nfs/users/nfs_j/jr46/git_repos/${PROJECT_NAME}/docs/mapping"
mkdir -p "${PERMANENT_RESULTS}"

if [[ -f "$STAR_STATS_FILE" ]]; then
    cp "$STAR_STATS_FILE" "${PERMANENT_RESULTS}/"
    echo "STAR stats copied to: ${PERMANENT_RESULTS}/star_mapping_summary.csv"
fi

if [[ -f "$BOWTIE_STATS_FILE" ]]; then
    cp "$BOWTIE_STATS_FILE" "${PERMANENT_RESULTS}/"
    echo "Bowtie2 stats copied to: ${PERMANENT_RESULTS}/bowtie2_mapping_summary.csv"
fi

tg_send "Mapping statistics recovery completed for ${PROJECT_NAME}" 2>/dev/null || true

echo ""
echo "========================================="
echo "  Mapping Statistics Recovery Complete"
echo "  STAR stats: ${STAR_STATS_FILE}"
echo "  Bowtie2 stats: ${BOWTIE_STATS_FILE}"
echo "  Permanent storage: ${PERMANENT_RESULTS}"
echo "========================================="
