#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/multiqc_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

bsub -o "${DATE_DIR}/multiqc_%J.out" \
     -e "${DATE_DIR}/multiqc_%J.err" \
     -q normal \
     -n 1 \
     -M 2000 \
     -R "select[mem>2000] rusage[mem=2000]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/multiqc_analysis.sh"

echo "MultiQC job submitted!"
