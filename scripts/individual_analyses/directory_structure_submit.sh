#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/individual_analyses/directory_structure.sh

# Set the date
DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"

# Create dated log directory
mkdir -p "${DATE_DIR}"

# Submit the job
bsub -o "${DATE_DIR}/directory_structure_%J.out" \
     -e "${DATE_DIR}/directory_structure_%J.err" \
     -q normal \
     -n 1 \
     -M 1000 \
     -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/directory_structure.sh"
