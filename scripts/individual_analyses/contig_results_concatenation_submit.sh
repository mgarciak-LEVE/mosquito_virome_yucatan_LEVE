#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/contig_results_concatenation_submit.sh
#
# Submits one job per (mapper, read_type, sample) — the pipeline script
# pools all four assemblers per sample, then runs CAP3 + CD-HIT-EST.

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

PROJECT_NAME="mosquito_virome_yucatan_LEVE"
IDENT_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/${PROJECT_NAME}/results/identification"

# Assemblers we expect to find under each sample
ASSEMBLERS=(MEGAhit metaSPAdes metaviralSPAdes rnaSPAdes)

# Count (mapper, read_type, sample) tuples — one job each
tuple_count=0

count_tuples() {
    local mapping_tool=$1
    local base_dir="${IDENT_DIR}/${mapping_tool}"

    echo "  Checking: ${mapping_tool}"

    if [[ ! -d "$base_dir" ]]; then
        echo "    WARNING: Directory not found: ${base_dir}"
        return
    fi

    # Walk read_type → sample directly (no assembler level in the loop)
    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type
        read_type=$(basename "$read_type_dir")

        # Enumerate samples using any assembler that has an identification dir
        local sample_root=""
        for asm in "${ASSEMBLERS[@]}"; do
            if [[ -d "${read_type_dir}${asm}" ]]; then
                sample_root="${read_type_dir}${asm}"
                break
            fi
        done
        [[ -n "$sample_root" ]] || continue

        for sample_dir in "$sample_root"/*/; do
            [[ -d "$sample_root" ]] || continue
            local sample
            sample=$(basename "$sample_dir")

            # At least one assembler must have produced something for this sample
            local any_found=0
            for asm in "${ASSEMBLERS[@]}"; do
                local asm_dir="${read_type_dir}${asm}/${sample}"
                if [[ -d "${asm_dir}/virsorter2" ]] \
                || [[ -d "${asm_dir}/deepvirfinder" ]] \
                || [[ -d "${asm_dir}/vitap" ]]; then
                    any_found=1
                    break
                fi
            done

            if [[ $any_found -eq 1 ]]; then
                tuple_count=$((tuple_count + 1))
                echo "      Found: ${mapping_tool}/${read_type}/${sample}"
            fi
        done
    done
}

echo ""
echo "Counting (mapper, read_type, sample) tuples..."
echo "------------------------"

count_tuples "STAR"
count_tuples "Bowtie"

echo "------------------------"
echo ""
echo "Total jobs to submit: ${tuple_count}"

if [[ $tuple_count -eq 0 ]]; then
    echo "ERROR: No samples found!"
    echo "Check directory: ${IDENT_DIR}"
    exit 1
fi

# Submit array job — one task per tuple
bsub -o "${DATE_DIR}/contig_results_concatenation_%J_%I.out" \
     -e "${DATE_DIR}/contig_results_concatenation_%J_%I.err" \
     -q normal \
     -n 1 \
     -M 2000 \
     -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" \
     -G team222 \
     -J "contig_concat[1-${tuple_count}]%20" \
     "${HOME}/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/contig_results_concatenation.sh"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   Contig concatenation + merge array job submitted successfully!"
    echo "   Submitted ${tuple_count} jobs (one per mapper/read_type/sample)"
else
    echo "Job submission failed!"
    exit 1
fi
