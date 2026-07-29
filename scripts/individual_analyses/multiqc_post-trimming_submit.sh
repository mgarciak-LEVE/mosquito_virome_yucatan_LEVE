#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/multiqc_post-trimming_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

bsub -o "${DATE_DIR}/multiqc_post-trimming_%J.out" \
     -e "${DATE_DIR}/multiqc_post-trimming_%J.err" \
     -q normal \
     -n 1 \
     -M 4000 \
     -R "select[mem>4000] rusage[mem=4000]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/multiqc_post-trimming.sh"

echo "Post-trimming MultiQC job submitted!"
