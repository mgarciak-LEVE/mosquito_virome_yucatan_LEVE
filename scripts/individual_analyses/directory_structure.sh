#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to generate direcotries and verify their existence.
# 23/06/2026
# Version 2.0.0

####==================================####
####          CONFIGURATION           ####
####==================================####

# Directory where scripts are located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Directorio 

# Telegram bot function
source "${SCRIPT_DIR}/bot_telegram.sh" 

# Proyect configuration 
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"
BASE_DIR="$(pwd)/${PROJECT_NAME}"

####==================================####
####  FUNCTION FOR DIRECTORY CREATION  ####
####==================================####

create_dir_validation() {
    local dir_path="$1"
    local dir_description="$2"
    
    if [[ -d "$dir_path" ]]; then
        echo "Directory already exists: ${dir_description}"
        tg_send "Directory already exists: ${dir_description}"

    else
        if mkdir -p "$dir_path"; then
            echo "New directory created: ${dir_description}"
            tg_send "New directory created: ${dir_description}"

        else
            echo "Error: ${dir_description}"
            tg_send "Error: ${dir_description}"
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
    echo "  Proyect: ${PROJECT_NAME}"
    echo "  Location: ${base}"
    echo "========================================="
    
    tg_send "Initiating directory structure creation: ${PROJECT_NAME}"
    
    # Base directory validation
    create_dir_validation "$base" "Base directory of the proyect"
    
    # Container directory
    create_dir_validation "${base}/containers" "/containers"

    # Data directories
    create_dir_validation "${base}/data/raw/total_RNA" "data/raw/total_RNA"
    create_dir_validation "${base}/data/raw/small_RNA" "data/raw/small_RNA"
    create_dir_validation "${base}/data/metadata" "data/metadata"
    
    # References directories
    create_dir_validation "${base}/data/references/mosquito_genomes/aedes_super_index" "References: Aedes supergenome"
    create_dir_validation "${base}/data/references/databases/BLAST" "References: BLAST DB"
    create_dir_validation "${base}/data/references/databases/DIAMOND" "References: DIAMOND DB"
    
    # Results directories
    create_dir_validation "${base}/results/untrimmed_qc/fastqc" "Results: FastQC (before trimming)"
    create_dir_validation "${base}/results/untrimmed_qc/multiqc" "Results: MultiQC (before trimming)"
    create_dir_validation "${base}/results/untrimmed_qc/stats" ".fastq file validation statistics (before trimming)"

    create_dir_validation "${base}/results/trimmed_qc/fastqc" "Results: FastQC (after trimming)"
    create_dir_validation "${base}/results/trimmed_qc/multiqc" "Results: MultiQC (after trimming)"
    create_dir_validation "${base}/results/trimmed_qc/stats" ".fastq file validation statistics (after trimming)"
    
    create_dir_validation "${base}/results/trimmed" "Results: Trimmed sequences"
    create_dir_validation "${base}/results/aligned" "Results: Alignments"
    create_dir_validation "${base}/results/assembly/statistics" "Results: Assembly statistics"
    create_dir_validation "${base}/results/assembly/rnaSPAdes" "Results: rnaSPAdes assembly"
    create_dir_validation "${base}/results/assembly/metaSPAdes" "Results: metaSPAdes assembly"
    create_dir_validation "${base}/results/assembly/MEGAhit" "Results: MEGAhit assembly"
    
    # Logs direcories
    create_dir_validation "${base}/logs/trimming" "Logs: Trimmomatic"
    create_dir_validation "${base}/logs/mapping" "Logs: STAR alignment"
    create_dir_validation "${base}/logs/assembly" "Logs: Assembly"
    create_dir_validation "${base}/logs/blast" "Logs: BLAST/DIAMOND"
    
    # Scripts and documents directories
    create_dir_validation "${base}/docs/aedes_genomes_specs" "Documentation: Genomes specs"
    create_dir_validation "${base}/scripts/aedes_reference_genomes" "Scripts: Reference genomes"
    create_dir_validation "${base}/scripts/databases" "Scripts: Databases"
    create_dir_validation "${base}/scripts/pipeline_whole" "Scripts: Complete pipeline"
    create_dir_validation "${base}/scripts/individual_analyses" "Scripts: Individual analyses"

}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    create_directory_structure "$BASE_DIR"
    
    echo ""
    echo "========================================="
    echo "Directory structure completed"
    echo "Proyect locatio : ${BASE_DIR}"
    echo "========================================="

fi