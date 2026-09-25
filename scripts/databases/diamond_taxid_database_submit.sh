#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/databases/diamond_taxid_database_submit.sh

# Set the date
DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"

# Create dated log directory
mkdir -p "${DATE_DIR}"

# Submit the job
bsub -o "${DATE_DIR}/diamond_taxid_database_%J.out" \
     -e "${DATE_DIR}/diamond_taxid_database_%J.err" \
     -q week \
     -n 32 \
     -M 128000 \
     -R "select[mem>128000] rusage[mem=128000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/databases/diamond_taxid_database.sh"
