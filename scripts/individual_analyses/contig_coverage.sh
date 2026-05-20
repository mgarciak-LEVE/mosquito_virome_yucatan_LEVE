#!/bin/bash

# Script for contig coverage analysis after read mapping.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 18/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw/czid_raw"
MAPPING_DIR="${RAW_DATA_DIR}/mapping"

# Virus names
VIRUS_NAMES=(
    "OSFV"
    "TAV"
    "GMMLV"
)

# Telegram bot
source "${SCRIPT_DIR}/bot_telegram.sh"

source "$(conda info --base)/etc/profile.d/conda.sh"

echo "=========================================="
echo "  Coverage Analysis"
echo "=========================================="
echo "Mapping directory: ${MAPPING_DIR}"
tg_send "Starting coverage analysis"

####==================================####
####     PROCESS EACH VIRUS TYPE      ####
####==================================####

for VIRUS_NAME in "${VIRUS_NAMES[@]}"; do

    VIRUS_OUTPUT="${MAPPING_DIR}/${VIRUS_NAME}"

    echo "=========================================="
    echo "Virus: ${VIRUS_NAME}"
    echo "Directory: ${VIRUS_OUTPUT}"
    echo "=========================================="

    # Check if virus directory exists
    if [ ! -d "$VIRUS_OUTPUT" ]; then
        echo "WARNING: Directory not found for ${VIRUS_NAME}"
        tg_send "WARNING: No mapping directory for ${VIRUS_NAME}. Skipping."
        continue
    fi

####==================================####
####   PROCESS EACH SAMPLE DIRECTORY  ####
####==================================####

    for sample_dir in "${VIRUS_OUTPUT}"/*/; do

        sample_name=$(basename "$sample_dir")
        bam_file="${sample_dir}/${sample_name}_sorted.bam"

        # Skip if BAM file doesn't exist
        if [ ! -f "$bam_file" ]; then
            echo "WARNING: No BAM file for ${sample_name}. Skipping."
            continue
        fi

        echo "----------------------------------------"
        echo "Sample: ${sample_name}"
        echo "BAM: $(basename "$bam_file")"
        echo "----------------------------------------"

        ####============================####
        ####       COVERAGE DEPTH        ####
        ####============================####

        echo "Calculating per-base depth..."
        
        conda activate samtools_env

        # Generate per-base depth file
        samtools depth \
            -a \
            "$bam_file" \
            > "${sample_dir}/${sample_name}_depth.txt"

        if [ $? -ne 0 ]; then
            echo "  ERROR: Depth calculation failed for ${sample_name}"
            tg_send "ERROR: Depth failed for ${VIRUS_NAME}/${sample_name}"
            conda deactivate
            continue
        fi

        echo "  Depth file: ${sample_name}_depth.txt"

        ####============================####
        ####     COVERAGE STATISTICS     ####
        ####============================####

        echo "Calculating coverage statistics..."

        # Total positions in the contig
        total_positions=$(wc -l < "${sample_dir}/${sample_name}_depth.txt" | tr -d ' ')

        # Average depth
        avg_depth=$(awk '{sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}' \
            "${sample_dir}/${sample_name}_depth.txt")

        # Positions with depth >= 3
        covered_positions=$(awk '$3 >= 3 {count++} END {print count+0}' \
            "${sample_dir}/${sample_name}_depth.txt")

        # Positions with depth == 0 (gaps)
        gap_positions=$(awk '$3 == 0 {count++} END {print count+0}' \
            "${sample_dir}/${sample_name}_depth.txt")

        # Percentage covered at >= 3x
        if [ "$total_positions" -gt 0 ]; then
            pct_covered=$(awk -v c="$covered_positions" -v t="$total_positions" \
                'BEGIN {printf "%.2f", (c/t)*100}')
            pct_gaps=$(awk -v g="$gap_positions" -v t="$total_positions" \
                'BEGIN {printf "%.2f", (g/t)*100}')
        else
            pct_covered="0.00"
            pct_gaps="0.00"
        fi

        # Min and max depth
        min_depth=$(awk 'NR==1 {min=$3} $3<min {min=$3} END {print min+0}' \
            "${sample_dir}/${sample_name}_depth.txt")
        max_depth=$(awk 'NR==1 {max=$3} $3>max {max=$3} END {print max+0}' \
            "${sample_dir}/${sample_name}_depth.txt")

        conda deactivate

        ####============================####
        ####        SAVE STATISTICS      ####
        ####============================####

        # Write statistics to a summary file
        stats_file="${sample_dir}/${sample_name}_coverage_stats.txt"

        cat > "$stats_file" << EOF
        
Coverage Statistics for ${sample_name}
======================================
Total positions:      ${total_positions}
Average depth:        ${avg_depth}x
Minimum depth:        ${min_depth}x
Maximum depth:        ${max_depth}x
Positions >= 3x:      ${covered_positions}
Positions with gaps:  ${gap_positions}
Percent covered:      ${pct_covered}%
Percent gaps:         ${pct_gaps}%
======================================
EOF

        ####============================####
        ####        DISPLAY RESULTS      ####
        ####============================####

        echo ""
        echo "  ┌─────────────────────────────────┐"
        printf  "  │ %-31s │\n" "Coverage Summary"
        echo "  ├─────────────────────────────────┤"
        printf  "  │ Total positions:    %-11s │\n" "$total_positions"
        printf  "  │ Average depth:      %-11s │\n" "${avg_depth}x"
        printf  "  │ Min depth:          %-11s │\n" "${min_depth}x"
        printf  "  │ Max depth:          %-11s │\n" "${max_depth}x"
        printf  "  │ Positions >= 3x:    %-11s │\n" "$covered_positions"
        printf  "  │ Positions with 0x:  %-11s │\n" "$gap_positions"
        printf  "  │ Percent >= 3x:      %-10s%% │\n" "$pct_covered"
        printf  "  │ Percent gaps:       %-10s%% │\n" "$pct_gaps"
        echo "  └─────────────────────────────────┘"
        echo ""

        # Determine quality flag
        if [ "$(echo "$pct_covered >= 80" | bc -l 2>/dev/null || echo 0)" -eq 1 ] && [ "$(echo "$avg_depth >= 10" | bc -l 2>/dev/null || echo 0)" -eq 1 ]; then
            echo "  Quality: GOOD (>=80% coverage, >=10x depth)"
        elif [ "$(echo "$pct_covered >= 50" | bc -l 2>/dev/null || echo 0)" -eq 1 ] && [ "$(echo "$avg_depth >= 3" | bc -l 2>/dev/null || echo 0)" -eq 1 ]; then
            echo "  Quality: ACCEPTABLE (>=50% coverage, >=3x depth)"
        else
            echo "  Quality: LOW (<50% coverage or <3x depth)"
        fi

        # Notify
        tg_send "${VIRUS_NAME}/${sample_name}: Depth=${avg_depth}x, Coverage=${pct_covered}%, Gaps=${pct_gaps}%"

        echo ""

    done

    echo "Completed coverage analysis for ${VIRUS_NAME}"
    tg_send "Completed coverage for ${VIRUS_NAME}"

done

####==================================####
####            SUMMARY               ####
####==================================####

echo ""
echo "=========================================="
echo "  Coverage Analysis Complete"
echo "=========================================="
echo ""
echo "Output files per sample:"
echo "  - <sample>_depth.txt          (per-base depth)"
echo "  - <sample>_coverage_stats.txt (summary statistics)"
echo ""
echo "Location: ${MAPPING_DIR}/<virus>/<sample>/"

tg_send "Coverage analysis complete. Results in ${MAPPING_DIR}"