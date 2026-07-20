#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/sequence_quality_array_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

INPUT_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/data/raw/total_RNA/cat_files"
cd "$INPUT_DIR" || exit 1

file_count=$(ls -1 *.fastq.gz *.fastq 2>/dev/null | wc -l)

echo "Found ${file_count} files"

bsub -o "${DATE_DIR}/sequence_quality_%J_%I.out" \
     -e "${DATE_DIR}/sequence_quality_%J_%I.err" \
     -q normal \
     -n 1 \
     -M 4000 \
     -R "select[mem>4000] rusage[mem=4000]" \
     -G team222 \
     -J "fastqc[1-${file_count}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/sequence_quality.sh"
