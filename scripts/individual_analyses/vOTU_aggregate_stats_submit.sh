#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/vOTU_aggregate_stats_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

bsub -o "${DATE_DIR}/vOTU_aggregate_stats_%J.out" \
     -e "${DATE_DIR}/vOTU_aggregate_stats_%J.err" \
     -q normal \
     -n 1 \
     -M 2000 \
     -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" \
     -G team222 \
     -J "vOTU_aggregate_stats" \
     bash ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/vOTU_aggregate_stats.sh

