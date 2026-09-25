#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/mapping_submit.sh
# Submit mapping as a job array

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

INPUT_FASTQC="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/trimmed"
cd "$INPUT_FASTQC" || exit 1

sample_count=$(find . -maxdepth 1 -type d -name "PM*" | wc -l)

if [[ $sample_count -eq 0 ]]; then
    echo "ERROR: No sample directories found in ${INPUT_FASTQC}"
    exit 1
fi

echo "========================================="
echo "  Submitting Mapping  Array Job"
echo "  Date: ${DATE}"
echo "  Files: ${sample_count}"
echo "  Logs: ${DATE_DIR}"
echo "========================================="

# Submit array job
bsub -o "${DATE_DIR}/mapping_%J_%I.out" \
     -e "${DATE_DIR}/mapping_%J_%I.err" \
     -q long \
     -n 12 \
     -M 80000 \
     -R "select[mem>80000] rusage[mem=80000]" \
     -G team222 \
     -J "mapping[1-${sample_count}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/mapping.sh"
