#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to recover assembly statistics from rnaSPAdes, metaSPAdes, metaviralSPAdes, and MEGAhit
# 08/08/2026
# Ver. 2.2.0 (farm-ready - with parent directory tracking)

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

# STAR assembly directory (parent indicates mapping tool)
STAR_ASSEMBLY_DIR="${PROJECT_SCRATCH}/results/assembly/star_assembly"

# Bowtie assembly directory (parent indicates mapping tool)
BOWTIE_ASSEMBLY_DIR="${PROJECT_SCRATCH}/results/assembly/bowtie_assembly"

# Output directory
OUTPUT_DIR="${PROJECT_SCRATCH}/results/assembly/statistics"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
SEQKIT_CONTAINER="${CONTAINERS}/seqkit_2.8.2.sif"

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
echo "  STAR assemblies: ${STAR_ASSEMBLY_DIR}"
echo "  Bowtie assemblies: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting assembly statistics recovery for ${PROJECT_NAME}" 2>/dev/null || true

# Create output directory
mkdir -p "${OUTPUT_DIR}"

####================================####
####       STATISTICS FUNCTION       ####
####================================####

get_assembly_stats() {
    local contig_file="$1"
    local sample="$2"
    local assembler="$3"
    local file_desc="$4"
    local mapping_tool="$5"
    local stats_file="$6"
    
    # Check if file exists and is not empty
    if [[ ! -f "$contig_file" ]] || [[ ! -s "$contig_file" ]]; then
        echo "  WARNING: File not found or empty: $contig_file"
        return 1
    fi
    
    echo "  Processing $file_desc for $sample - $assembler (mapping: $mapping_tool)"
    
    # Get basic statistics with seqkit stats
    stats=$(apptainer exec \
        --bind "$(dirname "$contig_file"):/input:ro" \
        "$SEQKIT_CONTAINER" \
        seqkit stats "/input/$(basename "$contig_file")" -T 2>/dev/null | awk 'NR==2')
    
    if [[ -z "$stats" ]]; then
        echo "    WARNING: seqkit stats failed for $contig_file"
        return 1
    fi
    
    # Parse statistics
    num_contigs=$(echo "$stats" | awk '{print $4}')
    total_length=$(echo "$stats" | awk '{print $5}')
    min_length=$(echo "$stats" | awk '{print $6}')
    max_length=$(echo "$stats" | awk '{print $7}')
    avg_length=$(echo "$stats" | awk '{print $8}')
    
    # N50 calculation
    n50=$(apptainer exec \
        --bind "$(dirname "$contig_file"):/input:ro" \
        "$SEQKIT_CONTAINER" \
        seqkit fx2tab -n -l "/input/$(basename "$contig_file")" 2>/dev/null | \
        awk '
        {len[NR]=$2; sum+=$2} 
        END {
            if (NR==0) {print 0; exit}
            for(i=1;i<=NR;i++) {
                for(j=i+1;j<=NR;j++) {
                    if(len[i] < len[j]) {
                        temp=len[i]; len[i]=len[j]; len[j]=temp
                    }
                }
            }
            half=sum/2;
            running=0;
            for(i=1;i<=NR;i++) {
                running+=len[i];
                if(running>=half) {
                    print len[i];
                    exit;
                }
            }
        }')
    
    [[ -z "$n50" ]] && n50=0
    
    # GC content
    gc_percent=$(apptainer exec \
        --bind "$(dirname "$contig_file"):/input:ro" \
        "$SEQKIT_CONTAINER" \
        seqkit fx2tab -n -g "/input/$(basename "$contig_file")" 2>/dev/null | \
        awk '{sum_gc+=$2; count++} END {if(count>0) printf "%.2f", sum_gc/count; else print "0"}')
    
    [[ -z "$gc_percent" ]] && gc_percent=0
    
    # Write to CSV
    echo "$sample,$mapping_tool,$assembler,$file_desc,$num_contigs,$total_length,$max_length,$min_length,$avg_length,$n50,$gc_percent" >> "$stats_file"
    
    echo "    $file_desc: $num_contigs contigs, N50=$n50, GC=$gc_percent%"
}

####================================####
####          PROCESS ASSEMBLIES     ####
####================================####

process_assembly_directory() {
    local base_dir="$1"
    local mapping_tool="$2"
    local stats_file="$3"
    
    if [[ ! -d "$base_dir" ]]; then
        echo "  Directory not found: $base_dir"
        return 1
    fi
    
    local files_processed=0
    
    echo ""
    echo "=== Processing $mapping_tool assemblies ==="
    
    # rnaSPAdes
    local rnaspades_dir="${base_dir}/rnaSPAdes"
    if [[ -d "$rnaspades_dir" ]]; then
        for sample_dir in "$rnaspades_dir"/*/; do
            if [[ -d "$sample_dir" ]]; then
                sample=$(basename "$sample_dir")
                
                for contig_file in "$sample_dir/transcripts.fasta" \
                                   "$sample_dir/soft_filtered_transcripts.fasta" \
                                   "$sample_dir/hard_filtered_transcripts.fasta"; do
                    if [[ -f "$contig_file" ]] && [[ -s "$contig_file" ]]; then
                        file_desc=$(basename "$contig_file" .fasta)
                        get_assembly_stats "$contig_file" "$sample" "rnaSPAdes" "$file_desc" "$mapping_tool" "$stats_file"
                        ((files_processed++))
                    fi
                done
            fi
        done
    fi
    
    # metaSPAdes
    local metaspades_dir="${base_dir}/metaSPAdes"
    if [[ -d "$metaspades_dir" ]]; then
        for sample_dir in "$metaspades_dir"/*/; do
            if [[ -d "$sample_dir" ]]; then
                sample=$(basename "$sample_dir")
                
                for contig_file in "$sample_dir/contigs.fasta" "$sample_dir/scaffolds.fasta"; do
                    if [[ -f "$contig_file" ]] && [[ -s "$contig_file" ]]; then
                        file_desc=$(basename "$contig_file" .fasta)
                        get_assembly_stats "$contig_file" "$sample" "metaSPAdes" "$file_desc" "$mapping_tool" "$stats_file"
                        ((files_processed++))
                    fi
                done
            fi
        done
    fi
    
    # metaviralSPAdes
    local metaviralspades_dir="${base_dir}/metaviralSPAdes"
    if [[ -d "$metaviralspades_dir" ]]; then
        for sample_dir in "$metaviralspades_dir"/*/; do
            if [[ -d "$sample_dir" ]]; then
                sample=$(basename "$sample_dir")
                
                for contig_file in "$sample_dir/contigs.fasta" "$sample_dir/scaffolds.fasta"; do
                    if [[ -f "$contig_file" ]] && [[ -s "$contig_file" ]]; then
                        file_desc=$(basename "$contig_file" .fasta)
                        get_assembly_stats "$contig_file" "$sample" "metaviralSPAdes" "$file_desc" "$mapping_tool" "$stats_file"
                        ((files_processed++))
                    fi
                done
            fi
        done
    fi
    
    # MEGAhit
    local megahit_dir="${base_dir}/MEGAhit"
    if [[ -d "$megahit_dir" ]]; then
        for sample_dir in "$megahit_dir"/*/; do
            if [[ -d "$sample_dir" ]]; then
                sample=$(basename "$sample_dir")
                
                contig_file="$sample_dir/final.contigs.fa"
                if [[ -f "$contig_file" ]] && [[ -s "$contig_file" ]]; then
                    get_assembly_stats "$contig_file" "$sample" "MEGAhit" "final.contigs" "$mapping_tool" "$stats_file"
                    ((files_processed++))
                else
                    contig_file="$sample_dir/contigs.fa"
                    if [[ -f "$contig_file" ]] && [[ -s "$contig_file" ]]; then
                        get_assembly_stats "$contig_file" "$sample" "MEGAhit" "contigs" "$mapping_tool" "$stats_file"
                        ((files_processed++))
                    fi
                fi
                
                # Extract MEGAhit log stats
                local log_file="$sample_dir/megahit.log"
                if [[ -f "$log_file" ]]; then
                    local total_reads=$(grep -E "total reads:" "$log_file" | tail -1 | awk '{print $NF}' | sed 's/,//g' 2>/dev/null)
                    local k_list=$(grep -E "k list:" "$log_file" | tail -1 | sed 's/.*k list: //' 2>/dev/null)
                    local runtime=$(grep -E "ALL DONE. Time elapsed:" "$log_file" | tail -1 | sed 's/.*Time elapsed: //' 2>/dev/null)
                    
                    LOG_STATS_FILE="${OUTPUT_DIR}/megahit_log_stats.csv"
                    if [[ ! -f "$LOG_STATS_FILE" ]]; then
                        echo "Sample,Mapping_Tool,Total_Reads,Kmer_List,Runtime" > "$LOG_STATS_FILE"
                    fi
                    echo "$sample,$mapping_tool,${total_reads:-0},${k_list:-N/A},${runtime:-N/A}" >> "$LOG_STATS_FILE"
                fi
            fi
        done
    fi
    
    echo "  Total files processed for $mapping_tool: $files_processed"
    return $files_processed
}

####================================####
####          GENERATE STATS         ####
####================================####

echo ""
echo "Generating assembly statistics with seqkit..."
tg_send "Generating assembly statistics with seqkit..." 2>/dev/null || true

# Create CSV file with headers
STATS_FILE="${OUTPUT_DIR}/assembly_summary.csv"
echo "Sample,Mapping_Tool,Assembler,File,Num_Contigs,Total_Length,Max_Length,Min_Length,Avg_Length,N50,GC_Percent" > "$STATS_FILE"

total_processed=0

# Process STAR assemblies
process_assembly_directory "$STAR_ASSEMBLY_DIR" "STAR" "$STATS_FILE"
total_processed=$((total_processed + $?))

# Process Bowtie assemblies
process_assembly_directory "$BOWTIE_ASSEMBLY_DIR" "Bowtie" "$STATS_FILE"
total_processed=$((total_processed + $?))

####================================####
####        DISPLAY SUMMARY         ####
####================================####

echo ""
echo "========================================="
echo "  Assembly Statistics Summary"
echo "========================================="

if [[ $total_processed -eq 0 ]]; then
    echo "WARNING: No assembly files were processed!"
    echo "Checked directories:"
    echo "  - STAR: $STAR_ASSEMBLY_DIR"
    echo "  - Bowtie: $BOWTIE_ASSEMBLY_DIR"
else
    echo "Total files processed: $total_processed"
    echo ""
    echo "Statistics saved to: $STATS_FILE"
    echo ""
    echo "=== Assembly Statistics ==="
    column -t -s, "$STATS_FILE" 2>/dev/null || cat "$STATS_FILE"
fi

echo ""
echo "========================================="

####================================####
####     COPY RESULTS TO NFS         ####
####================================####

echo ""
echo "Copying results to permanent storage..."

PERMANENT_RESULTS="/nfs/users/nfs_j/jr46/git_repos/${PROJECT_NAME}/docs/assembly_stats"
mkdir -p "${PERMANENT_RESULTS}"

if [[ -f "$STATS_FILE" ]]; then
    cp "$STATS_FILE" "${PERMANENT_RESULTS}/"
    echo "Stats copied to: ${PERMANENT_RESULTS}/assembly_summary.csv"
    tg_send "Stats copied to NFS: ${PERMANENT_RESULTS}" 2>/dev/null || true
fi

if [[ -f "$LOG_STATS_FILE" ]]; then
    cp "$LOG_STATS_FILE" "${PERMANENT_RESULTS}/"
    echo "MEGAhit log stats copied to: ${PERMANENT_RESULTS}/megahit_log_stats.csv"
fi

echo ""
echo "========================================="
echo "  Assembly Statistics Recovery Complete"
echo "  Stats file: ${STATS_FILE}"
echo "  Permanent storage: ${PERMANENT_RESULTS}"
echo "========================================="

tg_send "Assembly statistics recovery completed for ${PROJECT_NAME}" 2>/dev/null || true