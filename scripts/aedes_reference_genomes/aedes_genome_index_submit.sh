#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/aedes_reference_genomes/aedes_genome_index_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

bsub -o "${DATE_DIR}/aedes_genome_index_%J.out" \
     -e "${DATE_DIR}/aedes_genome_index_%J.err" \
     -q long \
     -n 8 \
     -M 64000 \
     -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" \
     -G team222 \
     < "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/aedes_reference_genomes/aedes_genome_index.sh"
     
echo "STAR indexing job submitted!"
