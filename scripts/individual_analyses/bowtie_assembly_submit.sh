#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/bowtie_assembly_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

BOWTIE_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/aligned/Bowtie2_alignment"
cd "$BOWTIE_DIR" || exit 1

paired_count=$(ls -1 *_unmapped_mixed.1.fastq.gz 2>/dev/null | wc -l)
r1_count=$(ls -1 *_R1_unpaired_unmapped.fastq.gz 2>/dev/null | wc -l)
r2_count=$(ls -1 *_R2_unpaired_unmapped.fastq.gz 2>/dev/null | wc -l)
total_count=$((paired_count + r1_count + r2_count))

if [[ $total_count -eq 0 ]]; then
    echo "No Bowtie unmapped samples found"
    exit 1
fi

echo "Found ${total_count} Bowtie samples"
echo "Paired mixed: ${paired_count}, R1 unpaired: ${r1_count}, R2 unpaired: ${r2_count}"

bsub -o "${DATE_DIR}/bowtie_assembly_%J_%I.out" \
     -e "${DATE_DIR}/bowtie_assembly_%J_%I.err" \
     -q long \
     -n 32 \
     -M 128000 \
     -R "select[mem>128000] rusage[mem=128000] span[hosts=1]" \
     -G team222 \
     -J "bowtie_assembly[1-${total_count}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/bowtie_assembly.sh"

echo "Bowtie assembly job submitted!"
