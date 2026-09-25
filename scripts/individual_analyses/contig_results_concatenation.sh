#!/bin/bash
# contig_results_concatenation
# Collect and unify viral contigs per sample for downstream taxonomy.
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

IDENT_DIR="${PROJECT_SCRATCH}/results/identification"
ASM_BASE="${PROJECT_SCRATCH}/results/assembly"
OUTPUT_DIR="${PROJECT_SCRATCH}/results/contig_concatenation"

STAR_ASM_DIR="${ASM_BASE}/star_assembly"
BOWTIE_ASM_DIR="${ASM_BASE}/bowtie_assembly"

# Containers
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
SEQKIT_CONTAINER="${CONTAINERS}/seqtk_1.5.sif"  

# Thresholds
DVF_PVALUE=0.05
VITAP_REALMS='Riboviria|Duplodnaviria|Monodnaviria'

# Telegram
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

####==================================####
####          LOAD MODULES            ####
####==================================####

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Contig Concatenation"
echo "  Project:   ${PROJECT_NAME}"
echo "  Ident dir: ${IDENT_DIR}"
echo "  Asm dirs:  ${STAR_ASM_DIR}"
echo "             ${BOWTIE_ASM_DIR}"
echo "  Output:    ${OUTPUT_DIR}"
echo "  Date:      $(date)"
echo "========================================="

tg_send "Starting viral contig concatenation" 2>/dev/null || true

####==================================####
####          SANITY CHECKS           ####
####==================================####

if [[ ! -f "${SEQKIT_CONTAINER}" ]]; then
    echo "WARNING: seqkit container not found at ${SEQKIT_CONTAINER}"
    exit 1
else
    run_seqkit_grep() {
        # $1 = ID file, $2 = FASTA, $3 = output
        local id_file=$1
        local fasta=$2
        local out=$3

        local id_dir fasta_dir fasta_base
        id_dir=$(dirname "${id_file}")
        fasta_dir=$(dirname "${fasta}")
        fasta_base=$(basename "${fasta}")
        id_file_base=$(basename "${id_file}")

        apptainer exec \
            --bind "${id_dir}":/id:ro \
            --bind "${fasta_dir}":/asm:ro \
            "${SEQKIT_CONTAINER}" \
            seqkit grep -f "/id/${id_file_base}" "/asm/${fasta_base}" \
            > "${out}"
    }
fi

####==================================####
####     BUILD SAMPLE LIST            ####
####==================================####

declare -a MAPPING_ARRAY
declare -a ASSEMBLER_ARRAY
declare -a READ_TYPE_ARRAY
declare -a SAMPLE_ARRAY
declare -a ASSEMBLY_FILE_ARRAY

add_sample() {
    local mapping_tool=$1
    local read_type=$2
    local assembler=$3
    local sample=$4
    local assembly_file=$5

    if [[ -f "$assembly_file" ]] && [[ -s "$assembly_file" ]]; then
        MAPPING_ARRAY+=("$mapping_tool")
        ASSEMBLER_ARRAY+=("$assembler")
        READ_TYPE_ARRAY+=("$read_type")
        SAMPLE_ARRAY+=("$sample")
        ASSEMBLY_FILE_ARRAY+=("$assembly_file")
    else
        echo "  WARNING: missing assembly: $assembly_file"
    fi
}

process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2

    if [[ ! -d "$base_dir" ]]; then
        echo "  WARNING: Directory not found: $base_dir"
        return
    fi

    echo "  Processing mapping tool: $mapping_tool"

    # Read type iteration
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

                case $assembler in
                    "MEGAhit")
                        assembly_file="${sample_dir}/final.contigs.fa"
                        ;;
                    "metaSPAdes"|"metaviralSPAdes")
                        assembly_file="${sample_dir}/contigs.fasta"
                        ;;
                    "rnaSPAdes")
                        if [[ -f "${sample_dir}/hard_filtered_transcripts.fasta" ]]; then
                            assembly_file="${sample_dir}/hard_filtered_transcripts.fasta"
                        else
                            assembly_file="${sample_dir}/transcripts.fasta"
                        fi
                        ;;
                    *)
                        echo "        WARNING: Unknown assembler: $assembler"
                        continue
                        ;;
                esac

                add_sample "$mapping_tool" "$read_type" "$assembler" "$sample" "$assembly_file"
            done
        done
    done
}


echo "Building sample list..."
process_mapping_tool "STAR"   "$STAR_ASM_DIR"
process_mapping_tool "Bowtie" "$BOWTIE_ASM_DIR"

TOTAL_JOBS=${#SAMPLE_ARRAY[@]}
if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: no samples found"
    tg_send "ERROR: no samples found for concatenation" 2>/dev/null || true
    exit 1
fi
echo "  Total jobs: ${TOTAL_JOBS}"

####==================================####
####          GET ARRAY TASK          ####
####==================================####

JOB_INDEX=$((LSB_JOBINDEX - 1))
if [[ $JOB_INDEX -ge ${#SAMPLE_ARRAY[@]} ]]; then
    echo "ERROR: invalid job index ${LSB_JOBINDEX}"
    exit 1
fi

MAPPING_TOOL="${MAPPING_ARRAY[$JOB_INDEX]}"
ASSEMBLER="${ASSEMBLER_ARRAY[$JOB_INDEX]}"
READ_TYPE="${READ_TYPE_ARRAY[$JOB_INDEX]}"
SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"
ASSEMBLY_FILE="${ASSEMBLY_FILE_ARRAY[$JOB_INDEX]}"

SAMPLE_INPUT="${IDENT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
SAMPLE_OUTPUT="${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
mkdir -p "${SAMPLE_OUTPUT}"

echo "========================================="
echo "  Job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
echo "  Input:    ${SAMPLE_INPUT}"
echo "  Output:   ${SAMPLE_OUTPUT}"
echo "  Assembly: ${ASSEMBLY_FILE}"
echo "========================================="

tg_send "Concatenating: ${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true

####==================================####
####       VIRSORTER2                 ####
####==================================####

vs2_in="${SAMPLE_INPUT}/virsorter2/final-viral-combined.fa"
vs2_out="${SAMPLE_OUTPUT}/virsorter2.fa"
: > "${vs2_out}"

if [[ -s "${vs2_in}" ]]; then
    # Headers look like: >k141_123||score=0.95
    # Keep the full header for now; the merge step will normalise IDs.
    cp "${vs2_in}" "${vs2_out}"
    n=$(grep -c '^>' "${vs2_out}")
    echo "  VirSorter2: ${n} contigs -> $(basename "${vs2_out}")"
else
    echo "  VirSorter2: no input (empty output)"
fi

####==================================####
####       DEEPVIRFINDER              ####
####==================================####

dvf_in="${SAMPLE_INPUT}/deepvirfinder/${SAMPLE}_prediction.txt"
dvf_ids="${SAMPLE_OUTPUT}/.dvf_ids.txt"
dvf_out="${SAMPLE_OUTPUT}/deepvirfinder.fa"
: > "${dvf_out}"

if [[ -s "${dvf_in}" ]]; then
    # Columns: name  length  score  pvalue  (tab-separated)
    # Keep IDs with pvalue <= threshold
    awk -v thr="${DVF_PVALUE}" -F'\t' 'NR>1 && $4+0 <= thr {print $1}' \
        "${dvf_in}" > "${dvf_ids}"

    n_ids=$(wc -l < "${dvf_ids}")
    echo "  DeepVirFinder: ${n_ids} viral IDs (p <= ${DVF_PVALUE})"

    if [[ "${n_ids}" -gt 0 ]]; then
        run_seqkit_grep "${dvf_ids}" "${ASSEMBLY_FILE}" "${dvf_out}"
        n=$(grep -c '^>' "${dvf_out}" || true)
        echo "  DeepVirFinder: ${n} sequences extracted"
    fi
else
    echo "  DeepVirFinder: no input (empty output)"
fi

####==================================####
####       VITAP                      ####
####==================================####

vitap_in="${SAMPLE_INPUT}/vitap/best_determined_lineages.tsv"
vitap_ids="${SAMPLE_OUTPUT}/.vitap_ids.txt"
vitap_out="${SAMPLE_OUTPUT}/vitap.fa"
: > "${vitap_out}"

if [[ -s "${vitap_in}" ]]; then
    # Columns: Genome_ID  lineage  lineage_score  Confidence_level
    # Filter viral realms, then drop reference accessions.
    awk -F'\t' -v realms="${VITAP_REALMS}" \
        'NR>1 && $2 ~ realms {print $1}' \
        "${vitap_in}" \
        | grep -v -E '^[A-Z]{1,3}[0-9]{4,}\.[0-9]+$' \
        > "${vitap_ids}"

    n_ids=$(wc -l < "${vitap_ids}")
    echo "  VITAP: ${n_ids} viral contig IDs (realm filter)"

    if [[ "${n_ids}" -gt 0 ]]; then
        run_seqkit_grep "${vitap_ids}" "${ASSEMBLY_FILE}" "${vitap_out}"
        n=$(grep -c '^>' "${vitap_out}" || true)
        echo "  VITAP: ${n} sequences extracted"
    fi
else
    echo "  VITAP: no input (empty output)"
fi

####==================================####
####       MERGE -> all.fa            ####
####==================================####

all_out="${SAMPLE_OUTPUT}/all.fa"

# Python helper for the merge. Uses only stdlib.
python3 - "${vs2_out}" "${vitap_out}" "${dvf_out}" "${all_out}" <<'PY'
import sys
from pathlib import Path

vs2_path, vitap_path, dvf_path, out_path = sys.argv[1:5]

def base_id(header):
    # header includes the leading '>'
    return header[1:].split("||")[0].split()[0]

def read_fasta(path):
    """Yield (header, sequence) pairs. header includes '>'."""
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line
                parts = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)

# Order matters: later writes overwrite header choice only if we want that.
# We want VirSorter2's header to win on collision (it carries score info),
# so process DVF and VITAP first, then VirSorter2 last to override.
order = [
    ("DeepVirFinder", dvf_path),
    ("VITAP",         vitap_path),
    ("VirSorter2",    vs2_path),
]

merged = {}  # base_id -> {"header": str, "seq": str, "tools": set()}

for tool, path in order:
    if not Path(path).exists():
        continue
    for header, seq in read_fasta(path):
        bid = base_id(header)
        if not bid:
            continue
        if bid not in merged:
            merged[bid] = {"header": header, "seq": seq, "tools": set()}
        merged[bid]["tools"].add(tool)
        # If this tool is VirSorter2 and it collides, prefer its header
        if tool == "VirSorter2":
            merged[bid]["header"] = header
            merged[bid]["seq"] = seq

def wrap(seq, width=80):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

with open(out_path, "w") as fh:
    for bid, rec in merged.items():
        tools = ",".join(sorted(rec["tools"]))
        # Strip any existing ||... annotation to avoid stacking
        header = rec["header"].split("||")[0]
        fh.write(f"{header}||tools={tools}\n")
        fh.write(wrap(rec["seq"]) + "\n")

# Report
from collections import Counter
tool_counts = Counter()
for rec in merged.values():
    for t in rec["tools"]:
        tool_counts[t] += 1
print(f"  merge: {len(merged)} unique contigs")
for t, n in sorted(tool_counts.items()):
    print(f"    {t}: {n}")
PY

n_all=$(grep -c '^>' "${all_out}" || true)
echo "  all.fa: ${n_all} unique contigs"

# Cleanup temporary ID lists
rm -f "${dvf_ids}" "${vitap_ids}"

####==================================####
####          COMPLETION              ####
####==================================####

echo "========================================="
echo "  Completed ${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}"
echo "  Outputs in ${SAMPLE_OUTPUT}"
echo "========================================="

tg_send "Completed concatenation: ${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${ASSEMBLER}/${SAMPLE}" 2>/dev/null || true
