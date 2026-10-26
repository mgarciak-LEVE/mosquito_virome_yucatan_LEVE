#!/bin/bash
# Author: Jorge Alberto Castro Rodríguez
# Map unmapped reads against the merged vOTU catalogue per sample.
# Ver. 2.0.0  

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

INPUT_DIR_BASE="${PROJECT_SCRATCH}/results"

# Index tree: <mapper>/<read_type>/<sample>/vOTUs.1.bt2
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/indexed_vOTU/STAR"
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/indexed_vOTU/Bowtie"

# Unmapped read sets (flat per sample, no assembler)
READS_DIR="${INPUT_DIR_BASE}/aligned"
BOWTIE_INPUT_DIR="${READS_DIR}/Bowtie2_alignment"
STAR_INPUT_DIR="${READS_DIR}/STAR_alignment"

# Output
OUTPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/reads_vs_vOTU"

CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
BOWTIE_CONTAINER="${CONTAINERS}/bowtie_samtools.sif"

SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=8

echo "========================================="
echo "  Read Mapping to vOTU Files (per sample, no assembler)"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  Bowtie2 container: ${BOWTIE_CONTAINER}"
echo "  Index dirs:"
echo "    - STAR:   ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output:     ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting vOTU mapping pipeline" 2>/dev/null || true

####==================================================####
####   COLLECT MAPPABLE (index, reads) PAIRS         ####
####==================================================####

all_sources=()
all_samples=()
all_types=()
all_read_types=()
all_input_dirs=()

register() {
    local src=$1 samp=$2 typ=$3 rt=$4 dir=$5
    all_sources+=("$src")
    all_samples+=("$samp")
    all_types+=("$typ")
    all_read_types+=("$rt")
    all_input_dirs+=("$dir")
}

INDEX_TOOL_MAP=(
    "Bowtie:${BOWTIE_ASSEMBLY_DIR}:${BOWTIE_INPUT_DIR}"
    "STAR:${STAR_ASSEMBLY_DIR}:${STAR_INPUT_DIR}"
)

for entry in "${INDEX_TOOL_MAP[@]}"; do
    IFS=":" read -r tool index_root read_root <<< "$entry"

    [[ -d "$index_root" ]] || { echo "  skip ${tool}: no index dir $index_root"; continue; }

    # Two levels only: read_type -> sample
    for read_type_dir in "$index_root"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        read_type=$(basename "$read_type_dir")

        for sample_dir in "$read_type_dir"/*/; do
            [[ -d "$sample_dir" ]] || continue
            sample=$(basename "$sample_dir")

            idx="${sample_dir}/vOTUs.1.bt2"
            [[ -s "$idx" ]] || continue

            case "${tool}_${read_type}" in
                "Bowtie_paired")
                    [[ -s "${read_root}/${sample}_unmapped_mixed.1.fastq.gz" \
                    && -s "${read_root}/${sample}_unmapped_mixed.2.fastq.gz" ]] \
                      && register "Bowtie" "$sample" "paired_interleaved" "paired" "$read_root"
                    ;;
                "Bowtie_unpaired_R1")
                    [[ -s "${read_root}/${sample}_R1_unpaired_unmapped.fastq.gz" ]] \
                      && register "Bowtie" "$sample" "R1_unpaired" "unpaired_R1" "$read_root"
                    ;;
                "Bowtie_unpaired_R2")
                    [[ -s "${read_root}/${sample}_R2_unpaired_unmapped.fastq.gz" ]] \
                      && register "Bowtie" "$sample" "R2_unpaired" "unpaired_R2" "$read_root"
                    ;;
                "STAR_paired")
                    [[ -s "${read_root}/${sample}_paired_Unmapped.out.mate1.gz" \
                    && -s "${read_root}/${sample}_paired_Unmapped.out.mate2.gz" ]] \
                      && register "STAR" "$sample" "paired_interleaved" "paired" "$read_root"
                    ;;
                "STAR_unpaired_R1")
                    [[ -s "${read_root}/${sample}_R1_unpaired_Unmapped.out.mate1.gz" ]] \
                      && register "STAR" "$sample" "R1_unpaired" "unpaired_R1" "$read_root"
                    ;;
                "STAR_unpaired_R2")
                    [[ -s "${read_root}/${sample}_R2_unpaired_Unmapped.out.mate1.gz" ]] \
                      && register "STAR" "$sample" "R2_unpaired" "unpaired_R2" "$read_root"
                    ;;
                *)
                    echo "  WARN: ${tool}/${read_type} has no read pattern; skipping ${sample}"
                    ;;
            esac
        done
    done
done

TOTAL=${#all_samples[@]}
if [[ $TOTAL -eq 0 ]]; then
    echo "ERROR: no mappable (index, reads) pairs found"
    tg_send "ERROR: no mappable pairs found" 2>/dev/null || true
    exit 1
fi

# Array index
if [[ -n "${LSB_JOBINDEX:-}" ]]; then
    IDX=$((LSB_JOBINDEX - 1))
else
    IDX=0
fi

if [[ $IDX -ge $TOTAL ]]; then
    echo "ERROR: job index ${LSB_JOBINDEX} out of range (max ${TOTAL})"
    exit 1
fi

source_tool="${all_sources[$IDX]}"
sample="${all_samples[$IDX]}"
sample_type="${all_types[$IDX]}"
read_type="${all_read_types[$IDX]}"
input_dir="${all_input_dirs[$IDX]}"

echo "========================================="
echo "  Source:    ${source_tool}"
echo "  Sample:    ${sample}"
echo "  Type:      ${sample_type}"
echo "  Read type: ${read_type}"
echo "  Task:      $((IDX+1))/${TOTAL}"
echo "========================================="

tg_send "Mapping ${source_tool}/${sample} (${read_type}) [$((IDX+1))/${TOTAL}]" 2>/dev/null || true

####==================================================####
####   RESOLVE READ PATHS                            ####
####==================================================####

case "${source_tool}_${sample_type}" in
    "Bowtie_paired_interleaved")
        R1="${input_dir}/${sample}_unmapped_mixed.1.fastq.gz"
        R2="${input_dir}/${sample}_unmapped_mixed.2.fastq.gz"
        MODE="pe"
        ;;
    "Bowtie_R1_unpaired")
        R1="${input_dir}/${sample}_R1_unpaired_unmapped.fastq.gz"
        R2=""
        MODE="se"
        ;;
    "Bowtie_R2_unpaired")
        R1="${input_dir}/${sample}_R2_unpaired_unmapped.fastq.gz"
        R2=""
        MODE="se"
        ;;
    "STAR_paired_interleaved")
        R1="${input_dir}/${sample}_paired_Unmapped.out.mate1.gz"
        R2="${input_dir}/${sample}_paired_Unmapped.out.mate2.gz"
        MODE="pe"
        ;;
    "STAR_R1_unpaired")
        R1="${input_dir}/${sample}_R1_unpaired_Unmapped.out.mate1.gz"
        R2=""
        MODE="se"
        ;;
    "STAR_R2_unpaired")
        R1="${input_dir}/${sample}_R2_unpaired_Unmapped.out.mate1.gz"
        R2=""
        MODE="se"
        ;;
    *)
        echo "ERROR: unknown combination ${source_tool}/${sample_type}"
        exit 1
        ;;
esac

[[ -s "$R1" ]] || { echo "ERROR: missing/empty $R1"; exit 1; }
[[ "$MODE" == "pe" && ! -s "$R2" ]] && { echo "ERROR: missing/empty $R2"; exit 1; }

####==================================================####
####   RESOLVE vOTU INDEX (3-level)                  ####
####==================================================####

if [[ "$source_tool" == "STAR" ]]; then
    INDEX_BASE="${STAR_ASSEMBLY_DIR}/${read_type}/${sample}"
else
    INDEX_BASE="${BOWTIE_ASSEMBLY_DIR}/${read_type}/${sample}"
fi
INDEX="${INDEX_BASE}/vOTUs"

if [[ ! -f "${INDEX}.1.bt2" ]]; then
    echo "ERROR: Bowtie2 index not found: ${INDEX}.1.bt2"
    exit 1
fi
echo "Index prefix: ${INDEX}"

####==================================================####
####   OUTPUT PATHS (3-level)                        ####
####==================================================####

OUT_SUBDIR="${OUTPUT_DIR}/${source_tool}/${read_type}"
mkdir -p "$OUT_SUBDIR"

BAM_OUT="${OUT_SUBDIR}/${sample}_vOTU_sorted.bam"
LOG_OUT="${OUT_SUBDIR}/${sample}_bowtie2.log"

####==================================================####
####   RUN BOWTIE2                                   ####
####==================================================####

APPTAINER_BINDS="\
${PROJECT_SCRATCH}:${PROJECT_SCRATCH},\
${OUTPUT_DIR}:${OUTPUT_DIR}"

if [[ "$MODE" == "pe" ]]; then
    echo "Running Bowtie2 PE: $R1 + $R2"
    apptainer exec \
        --bind "${APPTAINER_BINDS}" \
        "${BOWTIE_CONTAINER}" \
        bash -c "
            set -euo pipefail
            bowtie2 \
                -x '${INDEX}' \
                -1 '${R1}' \
                -2 '${R2}' \
                --end-to-end \
                --very-sensitive \
                -p ${THREADS} \
                2> '${LOG_OUT}' \
              | samtools sort -@ 4 -o '${BAM_OUT}' -
            samtools index '${BAM_OUT}'
        "
else
    echo "Running Bowtie2 SE: $R1"
    apptainer exec \
        --bind "${APPTAINER_BINDS}" \
        "${BOWTIE_CONTAINER}" \
        bash -c "
            set -euo pipefail
            bowtie2 \
                -x '${INDEX}' \
                -U '${R1}' \
                --end-to-end \
                --very-sensitive \
                -p ${THREADS} \
                2> '${LOG_OUT}' \
              | samtools sort -@ 4 -o '${BAM_OUT}' -
            samtools index '${BAM_OUT}'
        "
fi

echo "Done: ${BAM_OUT}"
tg_send "Finished mapping ${source_tool}/${sample} (${read_type})" 2>/dev/null || true
