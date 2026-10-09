#!/bin/bash

# Author: Jorge Alberto Castro Rodríguez
# Script to map unamapped reads to vOTU sequences.
# 07/10/2026
# Ver. 1.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

set -uo pipefail

# Directory where scripts are.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Project configuration
PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

# --- STORAGE LOCATIONS ---
PERMANENT_BASE="/nfs/users/nfs_j/jr46"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"

# --- PROJECT DIRECTORIES ---
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

# Input directories ---
INPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
# STAR assembly directory (parent indicates mapping tool)
STAR_ASSEMBLY_DIR="${INPUT_DIR_BASE}/indexed_vOTU/STAR"
# Bowtie assembly directory (parent indicates mapping tool)
BOWTIE_ASSEMBLY_DIR="${INPUT_DIR_BASE}/indexed_vOTU/Bowtie"
# Unmapped reads directory 
READS_DIR="${INPUT_DIR_BASE}/aligned"
# Unmapped reads Bowtie
BOWTIE_INPUT_DIR="${READS_DIR}/Bowtie2_alignment"
# Unmapped reads STAR
STAR_INPUT_DIR="${READS_DIR}/STAR_alignment"

# Output directories ---
OUTPUT_DIR_BASE="${PROJECT_SCRATCH}/results"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/reads_vs_vOTU"

# Container configuration
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
BOWTIE_CONTAINER="${CONTAINERS}/bowtie_samtools.sif"

# Scripts directory
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"

# Telegram bot
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

#==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

THREADS=8 

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Read Mapping to vOTU Files Script"
echo "  Project: ${PROJECT_NAME}"
echo "  Threads: ${THREADS}"
echo "  Bowtie2 container: ${BOWTIE_CONTAINER}"
echo "  Input directories (mapping_tool/assembler/read_type/sample):"
echo "    - STAR: ${STAR_ASSEMBLY_DIR}"
echo "    - Bowtie: ${BOWTIE_ASSEMBLY_DIR}"
echo "  Output for Bowtie indexing: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "========================================="

tg_send "Starting vOTU indexing pipeline" 2>/dev/null || true

####==================================================####
####   COLLECT UNMAPPED READ SETS FROM BOTH BRANCHES  ####
####==================================================####

all_samples=()
all_types=()
all_read_types=()
all_sources=()
all_assemblers=()
all_input_dirs=()

register() {
    local src=$1 samp=$2 typ=$3 rt=$4 asm=$5 dir=$6
    all_sources+=("$src")
    all_samples+=("$samp")
    all_types+=("$typ")
    all_read_types+=("$rt")
    all_assemblers+=("$asm")
    all_input_dirs+=("$dir")
}

####==================================================####
####   ENUMERATE INDEXES -> BUILD JOB LIST            ####
####==================================================####

INDEX_TOOL_MAP=(
    "Bowtie:${BOWTIE_ASSEMBLY_DIR}:${BOWTIE_INPUT_DIR}"
    "STAR:${STAR_ASSEMBLY_DIR}:${STAR_INPUT_DIR}"
)

# For each (tool, index_dir, read_dir):
#   walk index_dir/<read_type>/<assembler>/<sample>/vOTUs.1.bt2
#   register the job

for entry in "${INDEX_TOOL_MAP[@]}"; do
    IFS=":" read -r tool index_root read_root <<< "$entry"

    [[ -d "$index_root" ]] || { echo "  skip ${tool}: no index dir $index_root"; continue; }

    for read_type_dir in "$index_root"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        read_type=$(basename "$read_type_dir")     # e.g. paired, unpaired_R1, unpaired_R2

        for asm_dir in "$read_type_dir"/*/; do
            [[ -d "$asm_dir" ]] || continue
            assembler=$(basename "$asm_dir")       # e.g. megahit, metaviral, metaspades, rnaspades

            for sample_dir in "$asm_dir"/*/; do
                [[ -d "$sample_dir" ]] || continue
                sample=$(basename "$sample_dir")

                idx="${sample_dir}/vOTUs.1.bt2"
                [[ -s "$idx" ]] || continue

                # Map (tool, read_type) -> expected read filename pattern
                case "${tool}_${read_type}" in
                    "Bowtie_paired")
                        r1="${read_root}/${sample}_unmapped_mixed.1.fastq.gz"
                        r2="${read_root}/${sample}_unmapped_mixed.2.fastq.gz"
                        [[ -s "$r1" && -s "$r2" ]] || continue
                        register "Bowtie" "$sample" "paired_interleaved" "paired" "$assembler" "$read_root"
                        ;;
                    "Bowtie_unpaired_R1")
                        r1="${read_root}/${sample}_R1_unpaired_unmapped.fastq.gz"
                        [[ -s "$r1" ]] || continue
                        register "Bowtie" "$sample" "R1_unpaired" "unpaired_R1" "$assembler" "$read_root"
                        ;;
                    "Bowtie_unpaired_R2")
                        r1="${read_root}/${sample}_R2_unpaired_unmapped.fastq.gz"
                        [[ -s "$r1" ]] || continue
                        register "Bowtie" "$sample" "R2_unpaired" "unpaired_R2" "$assembler" "$read_root"
                        ;;
                    "STAR_paired")
                        r1="${read_root}/${sample}_paired_Unmapped.out.mate1.gz"
                        r2="${read_root}/${sample}_paired_Unmapped.out.mate2.gz"
                        [[ -s "$r1" && -s "$r2" ]] || continue
                        register "STAR" "$sample" "paired_interleaved" "paired" "$assembler" "$read_root"
                        ;;
                    "STAR_unpaired_R1")
                        r1="${read_root}/${sample}_R1_unpaired_Unmapped.out.mate1.gz"
                        [[ -s "$r1" ]] || continue
                        register "STAR" "$sample" "R1_unpaired" "unpaired_R1" "$assembler" "$read_root"
                        ;;
                    "STAR_unpaired_R2")
                        r1="${read_root}/${sample}_R2_unpaired_Unmapped.out.mate1.gz"
                        [[ -s "$r1" ]] || continue
                        register "STAR" "$sample" "R2_unpaired" "unpaired_R2" "$assembler" "$read_root"
                        ;;
                    *)
                        # Index read_type doesn't match a known read filename pattern
                        echo "  WARN: ${tool}/${read_type} has no read-file mapping; skipping ${assembler}/${sample}"
                        ;;
                esac
            done
        done
    done
done

TOTAL=${#all_samples[@]}
if [[ $TOTAL -eq 0 ]]; then
    echo "ERROR: no unmapped read sets found in either branch"
    tg_send "ERROR: no unmapped reads found" 2>/dev/null || true
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
assembler="${all_assemblers[$IDX]}"
input_dir="${all_input_dirs[$IDX]}"

echo "========================================="
echo "  Source:    ${source_tool}"
echo "  Sample:    ${sample}"
echo "  Type:      ${sample_type}"
echo "  Assembler: ${assembler}"
echo "  Read type: ${read_type}"
echo "  Task:      $((IDX+1))/${TOTAL}"
echo "========================================="

tg_send "Mapping ${source_tool}/${sample} (${read_type}) [$((IDX+1))/${TOTAL}]" 2>/dev/null || true

case "${source_tool}_${sample_type}" in
    # ---------------- Bowtie ----------------
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

    # ---------------- STAR ----------------
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

####=====================================####
####          DETERMINE vOTU INDEX       ####
####=====================================####
# The vOTU index depends on mapping tool + read type + assembler + sample.
# Match the structure:
#   ${INPUT_DIR_BASE}/indexed_vOTU/<MappingTool>/<read_type>/<assembler>/<sample>/vOTUs
# So we need to know which assembler was used. 

if [[ "$source_tool" == "STAR" ]]; then
    INDEX_BASE="${STAR_ASSEMBLY_DIR}/${read_type}/${assembler}/${sample}"
else
    INDEX_BASE="${BOWTIE_ASSEMBLY_DIR}/${read_type}/${assembler}/${sample}"
fi
INDEX="${INDEX_BASE}/vOTUs"             # Bowtie2 prefix (vOTUs.1.bt2 etc.)

if [[ ! -f "${INDEX}.1.bt2" ]]; then
    echo "ERROR: Bowtie2 index not found: ${INDEX}.1.bt2"
    exit 1
fi
echo "Index prefix: ${INDEX}"

####=====================================####
####          OUTPUT PATHS               ####
####=====================================####
# Layout: reads_vs_vOTU/<Source>/<read_type>/<sample>_sorted.bam
OUT_SUBDIR="${OUTPUT_DIR}/${source_tool}/${assembler}/${read_type}"
mkdir -p "$OUT_SUBDIR"

BAM_OUT="${OUT_SUBDIR}/${sample}_vOTU_sorted.bam"
LOG_OUT="${OUT_SUBDIR}/${sample}_bowtie2.log"

####=====================================####
####          RUN BOWTIE2                ####
####=====================================####
# Bind the relevant host-side dirs into the container read-only,
# and the output dir read-write.
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
                --no-unal \
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
                --no-unal \
                2> '${LOG_OUT}' \
              | samtools sort -@ 4 -o '${BAM_OUT}' -
            samtools index '${BAM_OUT}'
        "
fi

echo "Done: ${BAM_OUT}"
tg_send "Finished mapping ${source_tool}/${sample} (${read_type})" 2>/dev/null || true

####==================================####
####          COMPLETION              ####
####==================================####

echo ""
echo "========================================="
echo "  Completed job ${LSB_JOBINDEX}/${TOTAL}"
echo "  Mapping Tool: ${source_tool}"
echo "  Assembler: ${assembler}"
echo "  Sample: ${INDEX}"
echo "  Output: ${BAM_OUT}"
echo "  Log: ${LOG_OUT}"
echo "============="
