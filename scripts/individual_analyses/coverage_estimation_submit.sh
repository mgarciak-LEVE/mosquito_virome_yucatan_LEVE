#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/coverage_estimation_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

PROJECT_NAME="mosquito_virome_yucatan_LEVE"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/reads_vs_vOTU"

# Count BAM files
file_count=0

# Function to count BAMs for a given mapping tool
# Layout: <mapping_tool>/<assembler>/<read_type>/<SAMPLE>_vOTU_sorted.bam
count_bams() {
    local mapping_tool=$1
    local base_dir="${INPUT_DIR_BASE}/${mapping_tool}"

    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    # Iterate assembler first
    for assembler_dir in "$base_dir"/*/; do
        [[ -d "$assembler_dir" ]] || continue
        assembler=$(basename "$assembler_dir")
        echo "    Processing assembler: $assembler"

        # Then read_type
        for read_type_dir in "$assembler_dir"/*/; do
            [[ -d "$read_type_dir" ]] || continue
            read_type=$(basename "$read_type_dir")
            echo "      Processing read type: $read_type"

            # BAM files are flat inside read_type/
            shopt -s nullglob
            for bam in "$read_type_dir"/*_vOTU_sorted.bam; do
                [[ -s "$bam" ]] || continue
                sample=$(basename "$bam" _vOTU_sorted.bam)
                file_count=$((file_count + 1))
                echo "      Found: ${mapping_tool}/${assembler}/${read_type}/${sample}"
            done
            shopt -u nullglob
        done
    done
}

# Count BAMs for each mapping tool
echo ""
echo "Counting vOTU BAM files..."
echo "------------------------"

# Count STAR mappings
count_bams "STAR"

# Count Bowtie mappings
count_bams "Bowtie"

echo "------------------------"
echo ""
echo "Total vOTU BAM files found: ${file_count}"

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No vOTU BAM files found!"
    echo "Check directory: ${INPUT_DIR_BASE}"
    echo "Expected structure: ${INPUT_DIR_BASE}/{STAR,Bowtie}/{assembler}/{read_type}/<SAMPLE>_vOTU_sorted.bam"
    exit 1
fi

# Submit array job
bsub -o "${DATE_DIR}/coverage_estimation_%J_%I.out" \
     -e "${DATE_DIR}/coverage_estimation_%J_%I.err" \
     -q normal \
     -n 24 \
     -M 2000 \
     -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" \
     -G team222 \
     -J "coverage_estimation[1-${file_count}]%15" \
     "${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses/coverage_estimation.sh"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   vOTU coverage estimation array job submitted successfully!"
    echo "   Submitted ${file_count} jobs"
else
    echo "Job submission failed!"
    exit 1
fi
