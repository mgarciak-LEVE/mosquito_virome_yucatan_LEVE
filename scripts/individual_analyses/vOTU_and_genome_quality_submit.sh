#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/vOTU_and_genome_quality_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

PROJECT_NAME="mosquito_virome_yucatan_LEVE"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results/contig_concatenation"

# Count assembly files
file_count=0

# Function to count assemblies for a given mapping tool and assembler
count_assemblies() {
    local mapping_tool=$1
    local base_dir="${INPUT_DIR_BASE}/${mapping_tool}"

    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    # Iterate read_type first
    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        read_type=$(basename "$read_type_dir")
        echo "    Processing read type: $read_type"

        for assembler_dir in "$read_type_dir"/*/; do
            [[ -d "$assembler_dir" ]] || continue
            assembler=$(basename "$assembler_dir")
            echo "      Processing assembler: $assembler"

            for sample_dir in "$assembler_dir"/*/; do
                [[ -d "$sample_dir" ]] || continue
                sample=$(basename "$sample_dir")
                if [[ -f "${sample_dir}/all.fa" && -s "${sample_dir}/all.fa" ]]; then
                    file_count=$((file_count + 1))
                    echo "      Found: ${mapping_tool}/${read_type}/${assembler}/${sample}"

                fi
            done
        done
   done
}

# Count assemblies for each mapping tool
echo ""
echo "Counting assembly files..."
echo "------------------------"

# Count STAR assemblies
count_assemblies "STAR"

# Count Bowtie assemblies
count_assemblies "Bowtie"

echo "------------------------"
echo ""
echo "Total assembly files found: ${file_count}"

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No assembly files found!"
    echo "Check directory: ${INPUT_DIR_BASE}"
    echo "Expected structure: ${INPUT_DIR_BASE}/{star_assembly,bowtie_assembly}/{read_type}/{assembler}/{sample}/all.fa"
    exit 1
fi

# Submit array job
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
    echo "   vOTU and genome quality  array job submitted successfully!"
    echo "   Submitted ${file_count} jobs"
else
    echo "Job submission failed!"
    exit 1
fi
