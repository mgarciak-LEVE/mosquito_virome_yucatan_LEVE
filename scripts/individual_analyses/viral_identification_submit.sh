#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/viral_identification_submit.sh

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
    
    # For each assembler directory (MEGAhit, metaSPAdes, metaviralSPAdes, rnaSPAdes)
    for assembler_dir in "$base_dir"/*/; do
        if [[ -d "$assembler_dir" ]]; then
            assembler=$(basename "$assembler_dir")
            echo "    Assembler: ${assembler}"
            
            # For each sample directory
            for sample_dir in "$assembler_dir"/*/; do
                if [[ -d "$sample_dir%/" ]]; then
                    sample=$(basename "$sample_dir")
                    
                    # Check based on assembler type
                    case $assembler in
                        "MEGAhit")
                            if [[ -f "${sample_dir}/final.contigs.fa" ]]; then
                                file_count=$((file_count + 1))
                                echo "      Found: ${mapping_tool}/${assembler}/${sample}"
                            else
                                echo "      WARNING: No final.contigs.fa in ${sample_dir}"
                            fi
                            ;;
                        "metaSPAdes")
                            if [[ -f "${sample_dir}/contigs.fasta" ]]; then
                                file_count=$((file_count + 1))
                                echo "      Found: ${mapping_tool}/${assembler}/${sample}"
                            else
                                echo "      WARNING: No contigs.fasta in ${sample_dir}"
                            fi
                            ;;
                        "metaviralSPAdes")
                            if [[ -f "${sample_dir}/contigs.fasta" ]]; then
                                file_count=$((file_count + 1))
                                echo "      Found: ${mapping_tool}/${assembler}/${sample}"
                            else
                                echo "      WARNING: No contigs.fasta in ${sample_dir}"
                            fi
                            ;;
                        "rnaSPAdes")
                            if [[ -f "${sample_dir}/hard_filtered_transcripts.fasta" ]]; then
                                file_count=$((file_count + 1))
                                echo "      Found: ${mapping_tool}/${assembler}/${sample} (hard_filtered)"
                            elif [[ -f "${sample_dir}/transcripts.fasta" ]]; then
                                file_count=$((file_count + 1))
                                echo "      Found: ${mapping_tool}/${assembler}/${sample} (transcripts)"
                            else
                                echo "      WARNING: No transcripts in ${sample_dir}"
                            fi
                            ;;
                        *)
                            echo "      WARNING: Unknown assembler: ${assembler}"
                            ;;
                    esac
                fi
            done
        fi
    done
}

# Count assemblies for each mapping tool
echo ""
echo "Counting assembly files..."
echo "------------------------"

# Count STAR assemblies
count_assemblies "star_assembly"

# Count Bowtie assemblies
count_assemblies "bowtie_assembly"

echo "------------------------"
echo ""
echo "Total assembly files found: ${file_count}"

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No assembly files found!"
    echo "Check directory: ${ASSEMBLY_DIR}"
    echo "Expected structure: ${ASSEMBLY_DIR}/{star_assembly,bowtie_assembly}/{assembler}/{sample}/*.fasta"
    exit 1
fi


bsub -o "${DATE_DIR}/viral_identification_%J_%I.out" \
     -e "${DATE_DIR}/viral_identification_%J_%I.err" \
     -q normal \
     -n 12 \
     -M 32768 \
     -R "select[mem>=32768] rusage[mem=32768] span[hosts=1]" \
     -G team222 \
     -J "viral_identification[1-${file_count}]%10" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/viral_identification.sh"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   Viral identification array job submitted successfully!"
else
    echo "Job submission failed!"
    exit 1
fi
