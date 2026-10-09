#!/bin/bash
# ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/mapping_vOTUs_submit.sh
#
# Submission wrapper for mapping_vOTUs.sh
# Counts (index + unmapped reads) pairs across Bowtie/STAR branches,
# then submits one LSF array job with that exact size.

set -uo pipefail

####==================================####
####           CONFIGURATION          ####
####==================================####

DATE=$(date +%d_%m_%Y)
DATE_DIR="${HOME}/lsf_logs/${DATE}"
mkdir -p "${DATE_DIR}"

PROJECT_NAME="mosquito_virome_yucatan_LEVE"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

INPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/indexed_vOTU/Bowtie"
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/indexed_vOTU/STAR"

READS_DIR="${INPUT_DIR_BASE}/aligned"
BOWTIE_INPUT_DIR="${READS_DIR}/Bowtie2_alignment"
STAR_INPUT_DIR="${READS_DIR}/STAR_alignment"

MAPPING_SCRIPT="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses/mapping_vOTUs.sh"

# LSF resources for the mapping array
LSF_QUEUE="normal"
LSF_THREADS=8                    # matches THREADS=8 inside the mapping script
LSF_MEM_MB=8000
LSF_GROUP="team222"
LSF_ARRAY_CONCURRENCY="20"       # %20 → at most 20 jobs running simultaneously

####==================================####
####          COUNT JOBS              ####
####==================================####

file_count=0

count_indexes() {
    local tool=$1
    local index_root=$2
    local read_root=$3

    if [[ ! -d "$index_root" ]]; then
        echo "  WARNING: index dir not found: $index_root"
        return
    fi
    if [[ ! -d "$read_root" ]]; then
        echo "  WARNING: read dir not found: $read_root"
        return
    fi

    echo "  Counting [${tool}] under ${index_root}"

    for read_type_dir in "$index_root"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        read_type=$(basename "$read_type_dir")

        for asm_dir in "$read_type_dir"/*/; do
            [[ -d "$asm_dir" ]] || continue
            assembler=$(basename "$asm_dir")

            for sample_dir in "$asm_dir"/*/; do
                [[ -d "$sample_dir" ]] || continue
                sample=$(basename "$sample_dir")

                idx="${sample_dir}/vOTUs.1.bt2"
                [[ -s "$idx" ]] || continue

                case "${tool}_${read_type}" in
                    "Bowtie_paired")
                        [[ -s "${read_root}/${sample}_unmapped_mixed.1.fastq.gz" \
                        && -s "${read_root}/${sample}_unmapped_mixed.2.fastq.gz" ]] \
                          && file_count=$((file_count + 1))
                        ;;
                    "Bowtie_unpaired_R1")
                        [[ -s "${read_root}/${sample}_R1_unpaired_unmapped.fastq.gz" ]] \
                          && file_count=$((file_count + 1))
                        ;;
                    "Bowtie_unpaired_R2")
                        [[ -s "${read_root}/${sample}_R2_unpaired_unmapped.fastq.gz" ]] \
                          && file_count=$((file_count + 1))
                        ;;
                    "STAR_paired")
                        [[ -s "${read_root}/${sample}_paired_Unmapped.out.mate1.gz" \
                        && -s "${read_root}/${sample}_paired_Unmapped.out.mate2.gz" ]] \
                          && file_count=$((file_count + 1))
                        ;;
                    "STAR_unpaired_R1")
                        [[ -s "${read_root}/${sample}_R1_unpaired_Unmapped.out.mate1.gz" ]] \
                          && file_count=$((file_count + 1))
                        ;;
                    "STAR_unpaired_R2")
                        [[ -s "${read_root}/${sample}_R2_unpaired_Unmapped.out.mate1.gz" ]] \
                          && file_count=$((file_count + 1))
                        ;;
                    *)
                        # unknown read_type name → same warn the mapping script emits
                        :
                        ;;
                esac
            done
        done
    done
}

echo ""
echo "Counting mappable (index, reads) pairs..."
echo "-----------------------------------------"
count_indexes "Bowtie" "$BOWTIE_ASSEMBLY_DIR" "$BOWTIE_INPUT_DIR"
count_indexes "STAR"   "$STAR_ASSEMBLY_DIR"   "$STAR_INPUT_DIR"
echo "-----------------------------------------"
echo "Total mappable jobs: ${file_count}"
echo ""

if [[ $file_count -eq 0 ]]; then
    echo "ERROR: No (index, unmapped-reads) pairs found."
    echo "Expected index layout:"
    echo "  ${BOWTIE_ASSEMBLY_DIR}/<read_type>/<assembler>/<sample>/vOTUs.1.bt2"
    echo "  ${STAR_ASSEMBLY_DIR}/<read_type>/<assembler>/<sample>/vOTUs.1.bt2"
    echo "Expected read layout:"
    echo "  ${BOWTIE_INPUT_DIR}/<sample>_unmapped_mixed.1.fastq.gz  (etc.)"
    echo "  ${STAR_INPUT_DIR}/<sample>_paired_Unmapped.out.mate1.gz (etc.)"
    exit 1
fi

####==================================####
####          SANITY CHECKS           ####
####==================================####

if [[ ! -x "$MAPPING_SCRIPT" ]]; then
    echo "WARNING: mapping script is not executable: $MAPPING_SCRIPT"
    echo "         Run: chmod +x $MAPPING_SCRIPT"
fi
if [[ ! -f "$MAPPING_SCRIPT" ]]; then
    echo "ERROR: mapping script not found: $MAPPING_SCRIPT"
    exit 1
fi

####==================================####
####          SUBMIT ARRAY            ####
####==================================####

JOB_NAME="mapping_vOTUs"
bsub -o "${DATE_DIR}/${JOB_NAME}_%J_%I.out" \
     -e "${DATE_DIR}/${JOB_NAME}_%J_%I.err" \
     -q "${LSF_QUEUE}" \
     -n "${LSF_THREADS}" \
     -M "${LSF_MEM_MB}" \
     -R "select[mem>=${LSF_MEM_MB}] rusage[mem=${LSF_MEM_MB}] span[hosts=1]" \
     -G "${LSF_GROUP}" \
     -J "${JOB_NAME}[1-${file_count}]%${LSF_ARRAY_CONCURRENCY}" \
     "$MAPPING_SCRIPT" "${PROJECT_NAME}"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "   vOTU mapping array job submitted successfully!"
    echo "   Jobs submitted : ${file_count}"
    echo "   Concurrency    : ${LSF_ARRAY_CONCURRENCY}"
    echo "   Logs           : ${DATE_DIR}/${JOB_NAME}_<jobid>_<arrayidx>.out/.err"
else
    echo "Job submission failed!"
    exit 1
fi
