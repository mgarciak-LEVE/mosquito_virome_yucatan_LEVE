#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/assembly_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

# We'll count combined samples from both STAR and Bowtie
# The script itself will handle the counting, so we'll just submit
# and let the script figure out how many tasks to run

echo "========================================="
echo "  Submitting Unmapped Assembly Array Job"
echo "  Date: ${DATE}"
echo "  Logs: ${DATE_DIR}"
echo "========================================="

# Count samples in STAR directory
STAR_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/aligned/STAR_alignment"
cd "$STAR_DIR" 2>/dev/null || exit 1
star_count=$(ls -1 *_paired_Unmapped.out.mate1 2>/dev/null | wc -l)
star_r1_count=$(ls -1 *_R1_unpaired_Unmapped.out.mate1 2>/dev/null | wc -l)
star_r2_count=$(ls -1 *_R2_unpaired_Unmapped.out.mate1 2>/dev/null | wc -l)
star_total=$((star_count + star_r1_count + star_r2_count))

# Count samples in Bowtie directory
BOWTIE_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/aligned/Bowtie2_alignment"
cd "$BOWTIE_DIR" 2>/dev/null || exit 1
bowtie_count=$(ls -1 *_both_unmapped.fastq 2>/dev/null | wc -l)
if [[ $bowtie_count -eq 0 ]]; then
    bowtie_count=$(ls -1 *_unmapped_mixed.fastq 2>/dev/null | wc -l)
fi
bowtie_r1_count=$(ls -1 *_R1_unpaired_unmapped.fastq 2>/dev/null | wc -l)
bowtie_r2_count=$(ls -1 *_R2_unpaired_unmapped.fastq 2>/dev/null | wc -l)
bowtie_total=$((bowtie_count + bowtie_r1_count + bowtie_r2_count))

total_samples=$((star_total + bowtie_total))

echo "Found ${total_samples} total unmapped samples"
echo "  - STAR: ${star_total}"
echo "  - Bowtie: ${bowtie_total}"

if [[ $total_samples -eq 0 ]]; then
    echo "No unmapped samples found!"
    exit 1
fi

bsub -o "${DATE_DIR}/assembly_%J_%I.out" \
     -e "${DATE_DIR}/assembly_%J_%I.err" \
     -q week \
     -n 32 \
     -M 128000 \
     -R "select[mem>128000] rusage[mem=128000]" \
     -G team222 \
     -J "assembly[1-${total_samples}]%4" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/assembly.sh"

echo ""
echo "Unmapped assembly array job submitted!"
echo "   - Total samples: ${total_samples}"
