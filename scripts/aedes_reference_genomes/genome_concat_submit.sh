#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/aedes_reference_genomes/genome_concat_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

bsub -o "${DATE_DIR}/genome_concat_%J.out" \
     -e "${DATE_DIR}/genome_concat_%J.err" \
     -q normal \
     -n 1 \
     -M 8000 \
     -R "select[mem>8000] rusage[mem=8000]" \
     -G team222 \
     < ${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/aedes_reference_genomes/genome_concat.sh

echo "Concatenation job submitted!"
