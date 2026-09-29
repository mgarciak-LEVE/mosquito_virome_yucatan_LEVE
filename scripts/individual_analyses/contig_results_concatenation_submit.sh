#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/contig_results_concatenation_submit.sh

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

PROJECT_NAME="mosquito_virome_yucatan_LEVE"
ASSEMBLY_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/${PROJECT_NAME}/results/assembly"

# Count assembly files
file_count=0

# Function to count assemblies for a given mapping tool and assembler
count_assemblies() {
    local mapping_tool=$1
    local base_dir="${ASSEMBLY_DIR}/${mapping_tool}"

    echo "  Checking: ${mapping_tool}"

    if [[ ! -d "$base_dir" ]]; then
        echo "    WARNING: Directory not found: ${base_dir}"
        return
    fi

    # For each read type (e.g. ribo, polyA)
    for read_type_dir in "$base_dir"/*/; do
        if [[ -d "$read_type_dir" ]]; then
            read_type=$(basename "$read_type_dir")

            # For each assembler directory
            for assembler_dir in "$read_type_dir"/*/; do
                if [[ -d "$assembler_dir" ]]; then
                    assembler=$(basename "$assembler_dir")

                    # For each sample directory
                    for sample_dir in "$assembler_dir"/*/; do
                        if [[ -d "$sample_dir" ]]; then
                            sample=$(basename "$sample_dir")

                            # Check based on assembler type
                            case $assembler in
                                "MEGAhit")
                                    if [[ -f "${sample_dir}/final.contigs.fa" && -s "${sample_dir}/final.contigs.fa" ]]; then
                                        file_count=$((file_count + 1))
                                        echo "      Found: ${mapping_tool}/${read_type}/${assembler}/${sample}"
                                    fi
                                    ;;
                                "metaSPAdes"|"metaviralSPAdes")
                                    if [[ -f "${sample_dir}/contigs.fasta" && -s "${sample_dir}/contigs.fasta" ]]; then
                                        file_count=$((file_count + 1))
                                        echo "      Found: ${mapping_tool}/${read_type}/${assembler}/${sample}"
                                    fi
                                    ;;
                                "rnaSPAdes")
                                    if [[ -f "${sample_dir}/hard_filtered_transcripts.fasta" && -s "${sample_dir}/hard_filtered_transcripts.fasta" ]] \
                                       || [[ -f "${sample_dir}/transcripts.fasta" && -s "${sample_dir}/transcripts.fasta" ]]; then
                                        file_count=$((file_count + 1))
                                        echo "      Found: ${mapping_tool}/${read_type}/${assembler}/${sample}"
                                    fi
                                    ;;
                            esac
                        fi
                    done
                fi
            done
        fi
    done
}

# Count assemblies for each mapping tool
echo ""
echo "Counting assembly files..."
echo "------------------------"

count_assemblies "star_assembly"
count_assemblies "bowtie_assembly"

echo "------------------------"
echo ""
echo "Total assembly files found: ${file_count}"

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No assembly files found!"
    echo "Check directory: ${ASSEMBLY_DIR}"
    exit 1
fi

# Submit array job
bsub -o "${DATE_DIR}/contig_results_concatenation_%J_%I.out" \
     -e "${DATE_DIR}/contig_results_concatenation_%J_%I.err" \
     -q normal \
     -n 1 \
     -M 2000 \
     -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" \
     -G team222 \
     -J "contig_concat[1-${file_count}]%20" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/contig_results_concatenation.sh"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   Contig concatenation array job submitted successfully!"
    echo "   Submitted ${file_count} jobs"
else
    echo "Job submission failed!"
    exit 1
fi
