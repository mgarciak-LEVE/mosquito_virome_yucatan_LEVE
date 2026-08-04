#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/trimming_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

# Count R1 files
INPUT_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/data/raw/total_RNA/cat_files"
cd "$INPUT_DIR" || exit 1
file_count=$(ls -1 *_R1.fastq *_R1.fastq.gz 2>/dev/null | wc -l)

echo "Found ${file_count} R1 files"

bsub -o "${DATE_DIR}/trimming_%J_%I.out" \
     -e "${DATE_DIR}/trimming_%J_%I.err" \
     -q normal \
     -n 4 \
     -M 8000 \
     -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
     -G team222 \
     -J "trimming_[1-${file_count}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/trimming.sh"

echo "Trimming array job submitted!"
