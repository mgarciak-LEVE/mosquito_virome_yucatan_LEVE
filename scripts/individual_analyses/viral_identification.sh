#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to identify assembled contigs
# 18/08/2026
# Ver. 1.3.1 (farm-ready with conda env fix)

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
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/assembly"

# STAR assembly directory (parent indicates mapping tool)
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/star_assembly"

# Bowtie assembly directory (parent indicates mapping tool)
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/bowtie_assembly"

# Output directory
OUTPUT_DIR="${PROJECT_SCRATCH}/results/identification"

# Database paths
VIRSORTER2_DB="/nfs/users/nfs_j/jr46/databases/virsorter2_db/db"
CHECKV_DB="/nfs/users/nfs_j/jr46/databases/checkv_db/checkv-db-v1.5"
DIAMOND_DB="/nfs/users/nfs_j/jr46/databases/diamond_db"
GRAVITY_DB="/nfs/users/nfs_j/jr46/databases/gravity_db"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
VIRSORTER2_CONTAINER="${CONTAINERS}/virsorter2.sif"
DEEPVIRFINDER_CONTAINER="${CONTAINERS}/deepvirfinder.sif"
CHECKV_CONTAINER="${CONTAINERS}/checkv.sif"
DIAMOND_CONTAINER="${CONTAINERS}/diamond.sif"
GRAVITY_CONTAINER="${CONTAINERS}/gravityv2_v2.2.sif"

# Conda environment path (created on the host)
# Create with: conda create -n virsorter2 -c conda-forge -c bioconda virsorter snakemake mamba
CONDA_ENV_VIRSORTER2="/software/treeoflife/conda/users/envs/team222/jr46/virsorter2"

# Gravity configuration
GRAVITY_EMAIL="jacr@iibiomedicas.unam.mx"  # Required for NCBI downloads
GRAVITY_TAXO_GROUP="Genus"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"
module load conda 2>/dev/null || echo "Conda module not available"

THREADS=12
CORES=6

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Identification Pipeline"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  Conda env: ${CONDA_ENV_VIRSORTER2}"
echo "  Input directories (mapping_tool/assembler/sample):"
echo "    - STAR: ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting viral identification pipeline" 2>/dev/null || true

####==================================####
####     BUILD SAMPLE-ASSEMBLER LIST  ####
####==================================####

# Create arrays to store sample-mapping-assembler combinations
declare -a SAMPLE_ARRAY
declare -a MAPPING_ARRAY
declare -a ASSEMBLER_ARRAY
declare -a ASSEMBLY_FILE_ARRAY

echo "Building sample-assembler list..."

# Function to add sample-assembler combination
add_combination() {
    local mapping_tool=$1
    local assembler=$2
    local sample=$3
    local assembly_file=$4
    
    if [[ -f "$assembly_file" ]] && [[ -s "$assembly_file" ]]; then
        SAMPLE_ARRAY+=("$sample")
        MAPPING_ARRAY+=("$mapping_tool")
        ASSEMBLER_ARRAY+=("$assembler")
        ASSEMBLY_FILE_ARRAY+=("$assembly_file")
        echo "  Added: $mapping_tool/$assembler/$sample - $(basename "$assembly_file")"
    else
        echo "  WARNING: Assembly file not found or empty for $mapping_tool/$assembler/$sample: $assembly_file"
    fi
}

# Function to process a mapping tool directory
process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2
    
    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi
    
    echo "  Processing mapping tool: $mapping_tool"
    
    # For each assembler under this mapping tool
    for assembler_dir in "$base_dir"/*/; do
        if [[ -d "$assembler_dir" ]]; then
            assembler=$(basename "$assembler_dir")
            echo "    Processing assembler: $assembler"
            
            # For each sample under this assembler
            for sample_dir in "$assembler_dir"/*/; do
                if [[ -d "$sample_dir" ]]; then
                    sample=$(basename "$sample_dir")
                    
                    # Determine which assembly file to use based on assembler
                    case $assembler in
                        "MEGAhit")
                            assembly_file="${sample_dir}/final.contigs.fa"
                            ;;
                        "metaSPAdes"|"metaviralSPAdes")
                            assembly_file="${sample_dir}/contigs.fasta"
                            ;;
                        "rnaSPAdes")
                            if [[ -f "${sample_dir}/hard_filtered_transcripts.fasta" ]]; then
                                assembly_file="${sample_dir}/hard_filtered_transcripts.fasta"
                            else
                                assembly_file="${sample_dir}/transcripts.fasta"
                            fi
                            ;;
                        *)
                            echo "      WARNING: Unknown assembler: $assembler"
                            continue
                            ;;
                    esac
                    
                    add_combination "$mapping_tool" "$assembler" "$sample" "$assembly_file"
                fi
            done
        fi
    done
}

# Process STAR assemblies
process_mapping_tool "STAR" "$STAR_ASSEMBLY_DIR"

# Process Bowtie assemblies
process_mapping_tool "Bowtie" "$BOWTIE_ASSEMBLY_DIR"

# Total number of jobs
TOTAL_JOBS=${#SAMPLE_ARRAY[@]}

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No valid assembly files found"
    tg_send "ERROR: No valid assembly files found for identification" 2>/dev/null || true
    exit 1
fi

echo ""
echo "========================================="
echo "  Total jobs: ${TOTAL_JOBS}"
echo "  Unique samples: $(printf '%s\n' "${SAMPLE_ARRAY[@]}" | sort -u | wc -l)"
echo "  Summary by mapping tool:"
echo "    - STAR: $(grep -c "STAR" <<< "${MAPPING_ARRAY[@]}" || echo "0")"
echo "    - Bowtie: $(grep -c "Bowtie" <<< "${MAPPING_ARRAY[@]}" || echo "0")"
echo "========================================="

####==================================####
####          GET ARRAY TASK          ####
####==================================####

# LSB_JOBINDEX is the array task ID (1-based)
JOB_INDEX=$((LSB_JOBINDEX - 1))

if [[ $JOB_INDEX -ge ${#SAMPLE_ARRAY[@]} ]]; then
    echo "ERROR: Invalid job index ${LSB_JOBINDEX} (max ${#SAMPLE_ARRAY[@]})"
    exit 1
fi

# Get sample, mapping tool, and assembler for this job
SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"
MAPPING_TOOL="${MAPPING_ARRAY[$JOB_INDEX]}"
ASSEMBLER="${ASSEMBLER_ARRAY[$JOB_INDEX]}"
ASSEMBLY_FILE="${ASSEMBLY_FILE_ARRAY[$JOB_INDEX]}"

echo "========================================="
echo "  Processing job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Sample: ${SAMPLE}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Assembly file: ${ASSEMBLY_FILE}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Processing: ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true

####==================================####
####          CREATE OUTPUT DIR       ####
####==================================####

# Output structure: mapping_tool/assembler/sample
SAMPLE_OUTPUT="${OUTPUT_DIR}/${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
mkdir -p "${SAMPLE_OUTPUT}"

# Create writable directories for VirSorter2
VS2_CACHE="${SAMPLE_OUTPUT}/virsorter2_cache"
VS2_TMP="${SAMPLE_OUTPUT}/virsorter2_tmp"
CONDA_PKGS="${SAMPLE_OUTPUT}/conda_pkgs"
mkdir -p "$VS2_CACHE" "$VS2_TMP" "$CONDA_PKGS"

# Change to output directory
cd "$SAMPLE_OUTPUT" || exit 1

run_virsorter2() {
    echo ""
    echo "--- Running VirSorter2 ---"
    
    local output_dir="${SAMPLE_OUTPUT}/virsorter2"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/final-viral-combined.fa" ]] && [[ -s "${output_dir}/final-viral-combined.fa" ]]; then
        echo "  VirSorter2 already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running VirSorter2..."
    
    # Ensure CORES is set
    if [[ -z "${CORES}" ]]; then
        CORES=4  # Default to 4 cores if not set
        echo "  WARNING: CORES not set, defaulting to ${CORES}"
    fi

    # Check if database exists
    if [[ ! -d "${VIRSORTER2_DB}" ]] || [[ ! -f "${VIRSORTER2_DB}/hmm/viral/combined.hmm" ]]; then
        echo "  WARNING: VirSorter2 database not found or incomplete at ${VIRSORTER2_DB}"
        echo "  Expected to find: ${VIRSORTER2_DB}/hmm/viral/combined.hmm"
        echo "  Skipping VirSorter2 analysis"
        return 0
    fi

    local CONDA_ENV_PATH="/software/treeoflife/conda/users/envs/team222/jr46/virsorter2"

    # Verify the conda environment exists
    if [[ ! -d "${CONDA_ENV_PATH}" ]]; then
        echo "  ERROR: Conda environment not found at ${CONDA_ENV_PATH}"
        echo "  Please create it with: conda create -n virsorter2_env -c conda-forge -c bioconda virsorter snakemake mamba"
        return 1
    fi
    
    # Get the basename of the assembly file
    local assembly_basename=$(basename "$ASSEMBLY_FILE")
    
    # Create config file for VirSorter2
    cat > "${output_dir}/config.yaml" << EOF
input: /input/${assembly_basename}
output: /output/virsorter2
db_dir: /virsorter2_db
threads: ${CORES}
EOF
    
    echo "  Using conda environment: ${CONDA_ENV_PATH}"
    
    # Run VirSorter2 using the conda environment
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        --bind "${VIRSORTER2_DB}":/virsorter2_db:ro \
        --bind "${CONDA_ENV_PATH}":/conda_env \
        --env "TMPDIR=/tmp" \
        --env "CONDA_PKGS_DIRS=/opt/conda/pkgs" \
        "${VIRSORTER2_CONTAINER}" \
        /bin/bash -c "
            # Initialize conda inside the container
            source /opt/conda/etc/profile.d/conda.sh
            
            # Activate the mounted conda environment
            conda activate /conda_env
            
            # Verify VirSorter2 is available
            echo '  VirSorter2 version:'
            virsorter --version
            
            # Run VirSorter2 using Snakemake
            virsorter run \
                -i /input/${assembly_basename} \
                -w /output/virsorter2 \
                --db-dir /virsorter2_db \
                --include-groups dsDNAphage,ssDNA,ssRNA,dsRNA \
                --min-length 350 \
                --min-score 0.5 \
                -j ${CORES} \
                all
        " 2>&1 | tee "${output_dir}/virsorter2.log"

    if [[ $? -eq 0 ]] && [[ -f "${output_dir}/final-viral-combined.fa" ]]; then
        echo "VirSorter2 completed successfully"
        viral_count=$(grep -c '^>' "${output_dir}/final-viral-combined.fa" 2>/dev/null || echo "0")
        echo "  Found ${viral_count} viral contigs"
        return 0
    else
        echo "VirSorter2 failed or produced no output"
        return 1
    fi
}

####==================================####
####    FUNCTION: RUN DEEPVIRFINDER   ####
####==================================####

run_deepvirfinder() {
    echo ""
    echo "--- Running DeepVirFinder ---"
    
    local output_dir="${SAMPLE_OUTPUT}/deepvirfinder"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/${SAMPLE}_prediction.txt" ]] && [[ -s "${output_dir}/${SAMPLE}_prediction.txt" ]]; then
        echo "  DeepVirFinder already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running DeepVirFinder..."
    
    # Ensure CORES is set
    if [[ -z "${CORES}" ]]; then
        CORES=4
        echo "  WARNING: CORES not set, defaulting to ${CORES}"
    fi
    
    # Create writable directory for DeepVirFinder cache
    local dvf_cache="${SAMPLE_OUTPUT}/dvf_cache"
    mkdir -p "$dvf_cache"
    
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        --bind "${VS2_TMP}":/tmp \
        --bind "${dvf_cache}":/root/.dvf \
        --env "TMPDIR=/tmp" \
        "${DEEPVIRFINDER_CONTAINER}" \
        python /DeepVirFinder/dvf.py \
            -i "/input/$(basename "$ASSEMBLY_FILE")" \
            -o "/output/deepvirfinder" \
            -l 350 \
            -c "${CORES}" \
            2>&1 | tee "${output_dir}/deepvirfinder.log"
    
    if [[ $? -eq 0 ]]; then
        echo "  DeepVirFinder completed successfully"
        if [[ -f "${output_dir}/${SAMPLE}_prediction.txt" ]]; then
            awk '$3 <= 0.05' "${output_dir}/${SAMPLE}_prediction.txt" > \
                "${output_dir}/viral_contigs.txt"
            viral_count=$(wc -l < "${output_dir}/viral_contigs.txt")
            echo "  Found ${viral_count} viral contigs (p < 0.05)"
        fi
        return 0
    else
        echo "  DeepVirFinder failed"
        return 1
    fi
}

####==================================####
####      FUNCTION: RUN CHECKV       ####
####==================================####

run_checkv() {
    echo ""
    echo "--- Running CheckV ---"
    
    local output_dir="${SAMPLE_OUTPUT}/checkv"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/quality_summary.tsv" ]] && [[ -s "${output_dir}/quality_summary.tsv" ]]; then
        echo "  CheckV already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    echo "  Running CheckV..."
    
    # Ensure THREADS is set
    if [[ -z "${THREADS}" ]]; then
        THREADS=4
        echo "  WARNING: THREADS not set, defaulting to ${THREADS}"
    fi
    
    # Check if database exists
    if [[ ! -d "${CHECKV_DB}" ]]; then
        echo "  ERROR: CheckV database not found at ${CHECKV_DB}"
        echo "  Please download CheckV database or set the correct path"
        return 1
    fi
    
    # Determine which viral contigs to use for CheckV
    local viral_input="${ASSEMBLY_FILE}"
    if [[ -f "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" ]] && [[ -s "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" ]]; then
        viral_input="${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa"
        echo "  Using VirSorter2 output for CheckV analysis"
    elif [[ -f "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]] && [[ -s "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]]; then
        # Extract viral sequences from assembly based on DeepVirFinder results
        python3 -c "
import sys
with open('${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt') as f:
    contigs = set(line.split()[0] for line in f)
with open('${ASSEMBLY_FILE}') as f:
    write_seq = False
    for line in f:
        if line.startswith('>'):
            contig_id = line.strip()[1:].split()[0]
            write_seq = contig_id in contigs
        if write_seq:
            sys.stdout.write(line)
" > "${output_dir}/viral_contigs.fa"
        viral_input="${output_dir}/viral_contigs.fa"
        echo "  Using DeepVirFinder output for CheckV analysis"
    fi
    
    # Create CheckV temp directory
    local checkv_temp="${output_dir}/tmp"
    mkdir -p "$checkv_temp"
    
    apptainer exec \
        --bind "$(dirname "$viral_input")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        --bind "${CHECKV_DB}":/checkv_db:ro \
        --bind "${checkv_temp}":/tmp \
        --env "TMPDIR=/tmp" \
        "${CHECKV_CONTAINER}" \
        checkv end_to_end \
            "/input/$(basename "$viral_input")" \
            "/output/checkv" \
            -d /checkv_db \
            -t "${THREADS}" \
            2>&1 | tee "${output_dir}/checkv.log"
    
    if [[ $? -eq 0 ]]; then
        echo "  CheckV completed successfully"
        return 0
    else
        echo "  CheckV failed"
        return 1
    fi
}

####==================================####
####      FUNCTION: RUN DIAMOND      ####
####==================================####

run_diamond() {
    echo ""
    echo "--- Running DIAMOND ---"
    
    local output_dir="${SAMPLE_OUTPUT}/diamond"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/diamond_results.tsv" ]] && [[ -s "${output_dir}/diamond_results.tsv" ]]; then
        echo "  DIAMOND already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    # Check if database exists
    if [[ ! -d "${DIAMOND_DB}" ]] || [[ ! -f "${DIAMOND_DB}/viral_nr.dmnd" ]]; then
        echo "  WARNING: DIAMOND database not found at ${DIAMOND_DB}/viral_nr.dmnd"
        echo "  Skipping DIAMOND analysis - database needs to be created"
        echo "  To create: diamond makedb --in viral_proteins.faa -d viral_nr"
        return 0
    fi
    
    echo "  Running DIAMOND..."
    
    # Ensure THREADS is set
    if [[ -z "${THREADS}" ]]; then
        THREADS=4
        echo "  WARNING: THREADS not set, defaulting to ${THREADS}"
    fi
    
    # Create diamond temp directory
    local diamond_temp="${output_dir}/tmp"
    mkdir -p "$diamond_temp"
    
    apptainer exec \
        --bind "$(dirname "$ASSEMBLY_FILE")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        --bind "${DIAMOND_DB}":/db:ro \
        --bind "${diamond_temp}":/tmp \
        --env "TMPDIR=/tmp" \
        "${DIAMOND_CONTAINER}" \
        diamond blastx \
            -q "/input/$(basename "$ASSEMBLY_FILE")" \
            -d "/db/viral_nr.dmnd" \
            -o "/output/diamond/diamond_results.tsv" \
            --ultra-sensitive \
            --evalue 1e-5 \
            --threads "${THREADS}" \
            --outfmt 6 \
            2>&1 | tee "${output_dir}/diamond.log"
    
    if [[ $? -eq 0 ]]; then
        echo "  DIAMOND completed successfully"
        return 0
    else
        echo "  DIAMOND failed"
        return 1
    fi
}

####==================================####
####      FUNCTION: RUN GRAVITY       ####
####==================================####

run_gravity() {
    echo ""
    echo "--- Running GRAViTy-V2 ---"
    
    local output_dir="${SAMPLE_OUTPUT}/gravity"
    mkdir -p "$output_dir"
    
    # Check if already completed
    if [[ -f "${output_dir}/taxonomy_results.csv" ]] && [[ -s "${output_dir}/taxonomy_results.csv" ]]; then
        echo "  GRAViTy-V2 already completed for ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
        return 0
    fi
    
    # Get viral contigs (prioritize VirSorter2, then DeepVirFinder)
    local viral_contigs=""
    if [[ -f "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" ]] && [[ -s "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" ]]; then
        viral_contigs="${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa"
        echo "  Using VirSorter2 viral contigs for GRAViTy-V2 analysis"
    elif [[ -f "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]] && [[ -s "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]]; then
        # Extract sequences from assembly based on DeepVirFinder contig IDs
        local dvf_contigs="${output_dir}/dvf_viral_contigs.fa"
        python3 -c "
import sys
with open('${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt') as f:
    contigs = set(line.split()[0] for line in f)
with open('${ASSEMBLY_FILE}') as f:
    write_seq = False
    for line in f:
        if line.startswith('>'):
            contig_id = line.strip()[1:].split()[0]
            write_seq = contig_id in contigs
        if write_seq:
            sys.stdout.write(line)
" > "$dvf_contigs"
        viral_contigs="$dvf_contigs"
        echo "  Using DeepVirFinder viral contigs for GRAViTy-V2 analysis"
    else
        echo "  No viral contigs found for GRAViTy-V2 analysis"
        return 0
    fi
    
    if [[ ! -s "$viral_contigs" ]]; then
        echo "  No viral contigs found for GRAViTy-V2 analysis"
        return 0
    fi
    
    # Count viral contigs
    local viral_count=$(grep -c '^>' "$viral_contigs")
    echo "  Found ${viral_count} viral contigs for GRAViTy-V2 analysis"
    
    if [[ ${viral_count} -eq 0 ]]; then
        echo "  No viral contigs found"
        return 0
    fi
    
    echo "  Running GRAViTy-V2..."
    
    # Ensure THREADS is set
    if [[ -z "${THREADS}" ]]; then
        THREADS=4
        echo "  WARNING: THREADS not set, defaulting to ${THREADS}"
    fi
    
    # Determine which workflow to use based on viral contigs
    if [[ ${viral_count} -lt 50 ]]; then
        WORKFLOW="similar_viruses"
        echo "  Using 'similar_viruses' workflow (small dataset)"
    else
        WORKFLOW="new_classification_full"
        echo "  Using 'new_classification_full' workflow (large dataset)"
    fi
    
    # Create gravity temp directory
    local gravity_temp="${output_dir}/tmp"
    mkdir -p "$gravity_temp"
    
    # Run GRAViTy-V2
    apptainer exec \
        --bind "$(dirname "$viral_contigs")":/input:ro \
        --bind "${SAMPLE_OUTPUT}":/output \
        --bind "${GRAVITY_DB}":/gravity_db:ro \
        --bind "${gravity_temp}":/tmp \
        --env "TMPDIR=/tmp" \
        "${GRAVITY_CONTAINER}" \
        gravy_cli \
            --GenomeDescTableFile "/output/gravity/input_vmr.csv" \
            --GenomeSeqFile "/output/gravity/input_sequences.gb" \
            --ExpDir "/output/gravity" \
            --TaxoGrouping_Header "${GRAVITY_TAXO_GROUP}" \
            --genbank_email "${GRAVITY_EMAIL}" \
            --ProteinLength_Cutoff 70 \
            --UseBlast true \
            --NThreads "${THREADS}" \
            --Bootstrap_method "sumtrees" \
            --N_Bootstrap 10 \
            --workflow "${WORKFLOW}" \
            2>&1 | tee "${output_dir}/gravity.log"
    
    local exit_code=$?
    
    if [[ ${exit_code} -eq 0 ]]; then
        echo "  GRAViTy-V2 completed successfully"
        if [[ -f "${output_dir}/classification_results.csv" ]]; then
            echo "  Classification results generated"
            local classified=$(tail -n +2 "${output_dir}/classification_results.csv" | wc -l)
            echo "  Classified ${classified} viral contigs"
        fi
        return 0
    else
        echo "  GRAViTy-V2 failed with exit code ${exit_code}"
        return 1
    fi
}

####==================================####
####      FUNCTION: CONSOLIDATE       ####
####==================================####

consolidate_results() {
    echo ""
    echo "--- Consolidating Results ---"
    
    local report="${SAMPLE_OUTPUT}/consolidated_report.txt"
    
    cat > "$report" << EOF
===== VIRAL IDENTIFICATION REPORT =====
Mapping Tool: ${MAPPING_TOOL}
Assembler: ${ASSEMBLER}
Sample: ${SAMPLE}
Assembly file: $(basename "$ASSEMBLY_FILE")
Date: $(date)
Job: ${LSB_JOBINDEX}/${TOTAL_JOBS}

--- ASSEMBLY STATISTICS ---
$(if [[ -f "$ASSEMBLY_FILE" ]]; then
    contigs=$(grep -c '^>' "$ASSEMBLY_FILE" 2>/dev/null || echo "0")
    total_len=$(grep -v '^>' "$ASSEMBLY_FILE" | tr -d '\n' | wc -c 2>/dev/null || echo "0")
    echo "Total contigs: $contigs"
    echo "Total length: $total_len bp"
else
    echo "Assembly file not found"
fi)

--- VirSorter2 Results ---
$(if [[ -f "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" ]]; then
    viral_count=$(grep -c '^>' "${SAMPLE_OUTPUT}/virsorter2/final-viral-combined.fa" 2>/dev/null || echo "0")
    echo "Viral contigs found: $viral_count"
else
    echo "VirSorter2 results not available"
fi)

--- DeepVirFinder Results ---
$(if [[ -f "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" ]]; then
    viral_count=$(wc -l < "${SAMPLE_OUTPUT}/deepvirfinder/viral_contigs.txt" 2>/dev/null || echo "0")
    echo "Viral contigs found (p < 0.05): $viral_count"
else
    echo "DeepVirFinder results not available"
fi)

--- CheckV Completeness ---
$(if [[ -f "${SAMPLE_OUTPUT}/checkv/quality_summary.tsv" ]]; then
    echo "Completeness summary available"
    head -5 "${SAMPLE_OUTPUT}/checkv/quality_summary.tsv"
else
    echo "CheckV results not available"
fi)

=======================================
EOF
    
    echo "  Report created: $report"
}

####==================================####
####          RUN ALL TOOLS          ####
####==================================####

echo ""
echo "========================================="
echo "  Starting viral identification tools"
echo "  ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE}"
echo "========================================="

# Run each tool
run_virsorter2
run_deepvirfinder
run_checkv
run_diamond
run_gravity

# Consolidate results
consolidate_results

####==================================####
####          COMPLETION              ####
####==================================####

echo ""
echo "========================================="
echo "  Completed job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  Mapping Tool: ${MAPPING_TOOL}"
echo "  Assembler: ${ASSEMBLER}"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${SAMPLE_OUTPUT}"
echo "  Report: ${SAMPLE_OUTPUT}/consolidated_report.txt"
echo "========================================="

tg_send "Completed: ${MAPPING_TOOL}/${ASSEMBLER}/${SAMPLE} [${LSB_JOBINDEX}/${TOTAL_JOBS}]" 2>/dev/null || true
