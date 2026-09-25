#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/databases/uniref_build.sh

# Set the date
DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"

# Create dated log directory
mkdir -p "${DATE_DIR}"

bsub -q long -n 12 -M 64000 \
     -o "${DATE_DIR}/uniref_build.%J.out" \
     -e "${DATE_DIR}/uniref_build.%J.err" \
     -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/databases/uniref_build.sh"

