#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/quality_array_post-trimming_submit.sh
# Submit FastQC as a job array

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

INPUT_FASTQC="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/trimmed"
cd "$INPUT_FASTQC" || exit 1

# Count only paired files (or all .fastq files)
file_count=$(find . -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) 2>/dev/null | wc -l)

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in ${INPUT_FASTQC}"
    exit 1
fi

echo "========================================="
echo "  Submitting FastQC post-trimming Array Job"
echo "  Date: ${DATE}"
echo "  Files: ${file_count}"
echo "  Logs: ${DATE_DIR}"
echo "========================================="

# Submit array job
bsub -o "${DATE_DIR}/quality_array_%J_%I.out" \
     -e "${DATE_DIR}/quality_array_%J_%I.err" \
     -q normal \
     -n 4 \
     -M 4000 \
     -R "select[mem>4000] rusage[mem=4000]" \
     -G team222 \
     -J "quality_array[1-${file_count}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/quality_array_post-trimming.sh"
