#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/individual_analyses/directory_structure.sh

# Author: Jorge Alberto Castro Rodríguez
# Script to generate directories and verify their existence.
# 29/06/2026
# Version 2.1.1 (farm-ready)

####==================================####
####          CONFIGURATION           ####
####==================================####

# Directory where scripts are located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_pipeline}"  # Default project name if not provided as an argument

# --------- STORAGE LOCATIONS ---------

# Permanent storage
PERMANENT_BASE="/nfs/team222/projects"

# Scratch storage
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46"

# --------- PROJECT DIRECTORIES ---------

# Working directory on Lustre 
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Permanent results directory
PERMANENT_RESULTS="${PERMANENT_BASE}/${PROJECT_NAME}/results"

# LSF output directory
LSF_LOGS="${HOME}/lsf_logs/${DATE}"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

####==================================####
####  FUNCTION FOR DIRECTORY CREATION ####
####==================================####

create_dir_validation() {
    local dir_path="$1"
    local dir_description="$2"
    
    if [[ -d "$dir_path" ]]; then
        echo "Directory already exists: ${dir_description}"
    else
        if mkdir -p "$dir_path"; then
            echo "New directory created: ${dir_description}"
        else
            echo "ERROR: Failed to create ${dir_description} at ${dir_path}"
            exit 1
        fi
    fi
}

####==================================####
####       DIRECTORY CREATION         ####
####==================================####

create_directory_structure() {
    local base="$1"
    
    echo "========================================="
    echo "  Directory structure is being created"
    echo "  Project: ${PROJECT_NAME}"
    echo "  Location: ${base}"
    echo "========================================="
    
    # === Create LSF logs directory ===
    # This ensures LSF can write to it in future jobs
    create_dir_validation "${LSF_LOGS}" "LSF job logs (NFS)"
    
    # === Base directory on Lustre ===
    create_dir_validation "$base" "Project base directory (Lustre scratch)"
    
    # === Data directories ===
    create_dir_validation "${base}/data/raw/total_RNA" "data/raw/total_RNA"
    create_dir_validation "${base}/data/raw/small_RNA" "data/raw/small_RNA"
    create_dir_validation "${base}/data/metadata" "data/metadata"
    
    # === References directories ===
    create_dir_validation "${base}/data/references/mosquito_genomes/aedes_super_index" "References: Aedes supergenome"
    create_dir_validation "${base}/data/references/databases/BLAST" "References: BLAST DB"
    create_dir_validation "${base}/data/references/databases/DIAMOND" "References: DIAMOND DB"
    
    # === Permanent results directory ===
    create_dir_validation "${PERMANENT_RESULTS}" "Permanent results storage"

    # === Results directories ===
    create_dir_validation "${base}/results/untrimmed_qc/fastqc" "Results: FastQC (before trimming)"
    create_dir_validation "${base}/results/untrimmed_qc/multiqc" "Results: MultiQC (before trimming)"
    create_dir_validation "${base}/results/untrimmed_qc/stats" "Results: FASTQ validation stats (before trimming)"
    
    create_dir_validation "${base}/results/trimmed_qc/fastqc" "Results: FastQC (after trimming)"
    create_dir_validation "${base}/results/trimmed_qc/multiqc" "Results: MultiQC (after trimming)"
    create_dir_validation "${base}/results/trimmed_qc/stats" "Results: FASTQ validation stats (after trimming)"
    
    create_dir_validation "${base}/results/trimmed" "Results: Trimmed sequences"
    create_dir_validation "${base}/results/aligned" "Results: Alignments"
    create_dir_validation "${base}/results/assembly/statistics" "Results: Assembly statistics"
    create_dir_validation "${base}/results/assembly/rnaSPAdes" "Results: rnaSPAdes assembly"
    create_dir_validation "${base}/results/assembly/metaSPAdes" "Results: metaSPAdes assembly"
    create_dir_validation "${base}/results/assembly/MEGAhit" "Results: MEGAhit assembly"
    
    # === Logs directories ===
    create_dir_validation "${base}/logs/trimming" "Logs: Trimmomatic"
    create_dir_validation "${base}/logs/mapping" "Logs: STAR alignment"
    create_dir_validation "${base}/logs/assembly" "Logs: Assembly"
    create_dir_validation "${base}/logs/blast" "Logs: BLAST/DIAMOND"
    
    
    echo "========================================="
    echo "Directory structure completed"
    echo "Project location (Lustre): ${base}"
    echo "Permanent results: ${PERMANENT_RESULTS}"
    echo "LSF logs: ${LSF_LOGS}"
    echo "========================================="
}

# Create the main directory structure
create_directory_structure "$PROJECT_SCRATCH"

echo ""
echo "   - Work: ${PROJECT_SCRATCH}"
echo "   - Copy final results to: ${PERMANENT_RESULTS}"
echo "   - LSF logs: ${LSF_LOGS}"
echo "   - Symlink: ${HOME}/git_repos/mosquito_virome_pipeline/project_data"
echo ""
