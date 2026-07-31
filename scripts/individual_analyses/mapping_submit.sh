
#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/mapping_submit.sh
# Submit mapping as a job array

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
echo "  Submitting Mapping  Array Job"
echo "  Date: ${DATE}"
echo "  Files: ${file_count}"
echo "  Logs: ${DATE_DIR}"
echo "========================================="

# Submit array job
bsub -o "${DATE_DIR}/mapping_%J_%I.out" \
     -e "${DATE_DIR}/mapping_%J_%I.err" \
     -q long \
     -n 16 \
     -M 128000 \
     -R "select[mem>128000] rusage[mem=128000]" \
     -G team222 \
     -J "mapping[1-${file_count}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/mapping.sh"
