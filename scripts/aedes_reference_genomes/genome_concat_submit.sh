#!/bin/bash
# ~/git_repos/mosquito_virome_pipeline/scripts/aedes_reference_genomes/genome_concat_submit.sh

bsub -o ${HOME}/lsf_logs/%D/genome_concat_%J.out \
     -e ${HOME}/lsf_logs/%D/genome_concat_%J.err \
     -q normal -n 1 -M 8000 \
     -R "select[mem>8000] rusage[mem=8000]" \
     -G team222 \
     < ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/aedes_reference_genomes/genome_concat.sh
