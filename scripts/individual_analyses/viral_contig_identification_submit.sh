#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/viral_contig_identification_submit.sh

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

    # Level 1: sample_type (paired, unpaired_R1, unpaired_R2)
    for sample_type_dir in "$base_dir"/*/; do
        [[ -d "$sample_type_dir" ]] || continue
        local sample_type
        sample_type=$(basename "$sample_type_dir")
        echo "    Sample type: ${sample_type}"

        # Level 2: assembler (MEGAhit, metaSPAdes, metaviralSPAdes, rnaSPAdes)
        for assembler_dir in "$sample_type_dir"*/; do
            [[ -d "$assembler_dir" ]] || continue
            local assembler
            assembler=$(basename "$assembler_dir")
            echo "      Assembler: ${assembler}"

            # Level 3: sample
            for sample_dir in "$assembler_dir"*/; do
                [[ -d "$sample_dir" ]] || continue
                local sample
                sample=$(basename "$sample_dir")

                # Skip k-mer dirs (K21, K33, ...) just in case
                if [[ "$sample" =~ ^K[0-9]+$ ]]; then
                    continue
                fi

                local fasta=""
                case $assembler in
                    "MEGAhit")
                        fasta="${sample_dir}final.contigs.fa"
                        ;;
                    "metaSPAdes"|"metaviralSPAdes")
                        fasta="${sample_dir}contigs.fasta"
                        ;;
                    "rnaSPAdes")
                        if [[ -s "${sample_dir}hard_filtered_transcripts.fasta" ]]; then
                            fasta="${sample_dir}hard_filtered_transcripts.fasta"
                        elif [[ -s "${sample_dir}transcripts.fasta" ]]; then
                            fasta="${sample_dir}transcripts.fasta"
                        fi
                        ;;
                    *)
                        echo "        WARNING: Unknown assembler: ${assembler}"
                        continue
                        ;;
                esac

                if [[ -n "$fasta" && -s "$fasta" ]]; then
                    file_count=$((file_count + 1))
                    echo "        Found: ${mapping_tool}/${sample_type}/${assembler}/${sample}"
                else
                    echo "        WARNING: No valid assembly in ${sample_dir}"
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

# Submit array job
bsub -o "${DATE_DIR}/viral_contig_identification_%J_%I.out" \
     -e "${DATE_DIR}/viral_contig_identification_%J_%I.err" \
     -q long \
     -n 12 \
     -M 48000 \
     -R "select[mem>=48000] rusage[mem=48000] span[hosts=1]" \
     -G team222 \
     -J "viral_identification[1-${file_count}]%10" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/viral_contig_identification.sh"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   Viral identification array job submitted successfully!"
    echo "   Submitted ${file_count} jobs"
else
    echo "Job submission failed!"
    exit 1
fi
