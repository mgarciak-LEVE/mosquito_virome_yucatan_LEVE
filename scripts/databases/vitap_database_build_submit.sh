#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/databases/vitap_database_build_submit.sh

# Set the date
DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"

# Create dated log directory
mkdir -p "${DATE_DIR}"

# Submit the job
bsub -o "${DATE_DIR}/vitap_database_build_%J.out" \
     -e "${DATE_DIR}/vitap_database_build_%J.err" \
     -q long \
     -n 32 \
     -M 64000 \
     -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/databases/vitap_database_build.sh"
