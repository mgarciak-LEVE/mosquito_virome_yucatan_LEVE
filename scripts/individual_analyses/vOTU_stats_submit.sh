#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/vOTU_stats_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

PROJECT_NAME="mosquito_virome_yucatan_LEVE"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/vOTU"

# Count assembly files
file_count=0

# Function to count assemblies for a given mapping tool and assembler
count_samples() {
    local mapping_tool=$1
    local base_dir="${INPUT_DIR_BASE}/${mapping_tool}"

    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    # Two levels only: read_type -> sample
    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type
        read_type=$(basename "$read_type_dir")
        echo "    Processing read type: $read_type"

        for sample_dir in "$read_type_dir"/*/; do           # ← FIXED
            [[ -d "$sample_dir" ]] || continue
            local sample
            sample=$(basename "$sample_dir")

            if [[ -f "${sample_dir}/vOTUs.fa" && -s "${sample_dir}/vOTUs.fa" ]]; then
                file_count=$((file_count + 1))
                echo "      Found: ${mapping_tool}/${read_type}/${sample}"
            fi
        done
    done
}

# Count assemblies for each mapping tool
echo ""
echo "Counting OTU files..."
echo "------------------------"

# Count STAR assemblies
count_samples "STAR"

# Count Bowtie assemblies
count_samples "Bowtie"

echo "------------------------"
echo ""
echo "Total OTU files found: ${file_count}"

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No vOTU files found!"
    echo "Check directory: ${INPUT_DIR_BASE}"
    echo "Expected structure: ${INPUT_DIR_BASE}/{STAR,Bowtie}/{read_type}/{sample}/vOTUs.fa"
    exit 1
fi

# Submit array job
bsub -o "${DATE_DIR}/vOTU_stats_%J_%I.out" \
     -e "${DATE_DIR}/vOTU_stats_%J_%I.err" \
     -q normal \
     -n 1 \
     -M 4000 \
     -R "select[mem>=4000] rusage[mem=4000] span[hosts=1]" \
     -G team222 \
     -J "vOTU_stats[1-${file_count}]%15" \
     "${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses/vOTU_stats.sh"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   vOTU stats array job submitted successfully!"
    echo "   Submitted ${file_count} jobs"
else
    echo "Job submission failed!"
    exit 1
fi

