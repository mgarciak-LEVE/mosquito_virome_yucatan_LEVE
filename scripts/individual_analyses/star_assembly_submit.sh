#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/star_assembly_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

STAR_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/aligned/STAR_alignment"
cd "$STAR_DIR" || exit 1

sample_count=$(ls -1 *_paired_Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)

if [[ $sample_count -eq 0 ]]; then
    echo "No STAR paired samples found"
    exit 1
fi

echo "Found ${sample_count} STAR paired samples"

bsub -o "${DATE_DIR}/star_assembly_%J_%I.out" \
     -e "${DATE_DIR}/star_assembly_%J_%I.err" \
     -q long \
     -n 32 \
     -M 128000 \
     -R "select[mem>128000] rusage[mem=128000] span[hosts=1]" \
     -G team222 \
     -J "star_assembly[1-${sample_count}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/star_assembly.sh" \
     "mosquito_virome_yucatan_LEVE"

echo "STAR assembly job submitted!"