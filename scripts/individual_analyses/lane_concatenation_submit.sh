#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/individual_analyses/lane_concatenation_submit.sh

# Set the date
DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"

# Create dated log directory
mkdir -p "${DATE_DIR}"

# Submit the job
bsub -o "${DATE_DIR}/lane_concatenation_%J.out" \
     -e "${DATE_DIR}/lane_concatenation_%J.err" \
     -q normal \
     -n 16 \
     -M 32000 \
     -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/lane_concatenation.sh"
