#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/databases/viral_refseq_database_submit.sh

# Set the date
DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"

# Create dated log directory
mkdir -p "${DATE_DIR}"

# Submit the job
bsub -o "${DATE_DIR}/viral_refseq_database_%J.out" \
     -e "${DATE_DIR}/viral_refseq_database_%J.err" \
     -q normal \
     -n 8 \
     -M 16000 \
     -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/databases/viral_refseq_database.sh"




