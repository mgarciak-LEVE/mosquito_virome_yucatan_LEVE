#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/mapping_stats_submit.sh

# Set the date
DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"

# Create dated log directory
mkdir -p "${DATE_DIR}"

bsub -o "${DATE_DIR}/mapping_stats_%J.out" \
     -e "${DATE_DIR}/mapping_stats_%J.err" \
     -q normal \
     -n 1 \
     -M 4000 \
     -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/mapping_stats.sh"
