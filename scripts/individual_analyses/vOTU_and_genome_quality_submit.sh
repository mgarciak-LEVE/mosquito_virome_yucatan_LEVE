#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/vOTU_and_genome_quality_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

PROJECT_NAME="mosquito_virome_yucatan_LEVE"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/contig_concatenation"

file_count=0

count_samples() {
    local mapping_tool=$1
    local base_dir="${INPUT_DIR_BASE}/${mapping_tool}"

    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    # Two levels: read_type -> sample. No assembler level.
    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type
        read_type=$(basename "$read_type_dir")

        for sample_dir in "$read_type_dir"/*/; do
            [[ -d "$sample_dir" ]] || continue
            local sample
            sample=$(basename "$sample_dir")

            if [[ -f "${sample_dir}/all.fa" && -s "${sample_dir}/all.fa" ]]; then
                file_count=$((file_count + 1))
                echo "      Found: ${mapping_tool}/${read_type}/${sample}"
            fi
        done
    done
}

echo ""
echo "Counting merged contig files..."
echo "------------------------"

count_samples "STAR"
count_samples "Bowtie"

echo "------------------------"
echo ""
echo "Total merged contig files found: ${file_count}"

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No merged contig files found!"
    echo "Check directory: ${INPUT_DIR_BASE}"
    echo "Expected structure: ${INPUT_DIR_BASE}/{STAR,Bowtie}/{read_type}/{sample}/all.fa"
    exit 1
fi

bsub -o "${DATE_DIR}/vOTU_and_genome_quality_%J_%I.out" \
     -e "${DATE_DIR}/vOTU_and_genome_quality_%J_%I.err" \
     -q normal \
     -n 12 \
     -M 24000 \
     -R "select[mem>=24000] rusage[mem=24000] span[hosts=1]" \
     -G team222 \
     -J "vOTU_and_genome_quality[1-${file_count}]%10" \
     "${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses/vOTU_and_genome_quality.sh"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   vOTU and genome quality array job submitted successfully!"
    echo "   Submitted ${file_count} jobs (one per mapper/read_type/sample)"
else
    echo "Job submission failed!"
    exit 1
fi
