#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/star_assembly_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

STAR_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/aligned/STAR_alignment"
cd "$STAR_DIR" || exit 1

# Count all sample types (paired, R1_unpaired, R2_unpaired)
paired_count=$(ls -1 *_paired_Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)
r1_unpaired_count=$(ls -1 *_R1_unpaired_Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)
r2_unpaired_count=$(ls -1 *_R2_unpaired_Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)
total_samples=$((paired_count + r1_unpaired_count + r2_unpaired_count))

if [[ $total_samples -eq 0 ]]; then
    echo "No STAR samples found (paired or unpaired)"
    exit 1
fi

echo "Found ${total_samples} STAR samples:"
echo "  - Paired: ${paired_count}"
echo "  - R1 unpaired: ${r1_unpaired_count}"
echo "  - R2 unpaired: ${r2_unpaired_count}"

bsub -o "${DATE_DIR}/star_assembly_%J_%I.out" \
     -e "${DATE_DIR}/star_assembly_%J_%I.err" \
     -q long \
     -n 32 \
     -M 128000 \
     -R "select[mem>128000] rusage[mem=128000] span[hosts=1]" \
     -G team222 \
     -J "star_assembly[1-${total_samples}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/star_assembly.sh" \
     "mosquito_virome_yucatan_LEVE"

echo "STAR assembly job submitted!"
