#!/bin/bash
# contig_results_concatenation
# Collect and unify viral contigs per sample for downstream taxonomy. This version pools across assemblers with CAP3 and cdhit
# Author: Jorge Alberto Castro Rodríguez
# Ver. 2.0.0

####==================================####
####           CONFIGURATION          ####
####==================================####

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"

SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PROJECT_SCRATCH="${SCRATCH_BASE}/${PROJECT_NAME}"

IDENT_DIR="${PROJECT_SCRATCH}/results/identification"
OUTPUT_DIR="${PROJECT_SCRATCH}/results/contig_concatenation"

STAR_ASM_DIR="${IDENT_DIR}/STAR"
BOWTIE_ASM_DIR="${IDENT_DIR}/Bowtie"

# Assemblers to pool
ASSEMBLERS_TO_POOL=(MEGAhit metaSPAdes metaviralSPAdes rnaSPAdes)

# Containers
CONTAINERS="/lustre/scratch126/tol/teams/lawniczak/users/jr46/containers"
SEQKIT_CONTAINER="${CONTAINERS}/seqkit_2.8.2.sif"
CAP3_CONTAINER="${CONTAINERS}/cap3_10.2011--0.sif"                 
CDHIT_CONTAINER="${CONTAINERS}/cdhit_4.8.1--c697693.sif" 

# Thresholds
DVF_PVALUE=0.05
VITAP_REALMS='Riboviria|Duplodnaviria|Monodnaviria'
CAP3_P=90        # minimum overlap percent identity
CAP3_O=40        # minimum overlap length (bp)
CDHIT_C=0.95     # cluster identity
CDHIT_AS=0.85    # coverage of shorter sequence

# Telegram
SCRIPTS_NFS="${HOME}/git_repos/${PROJECT_NAME}/scripts/individual_analyses"
source "${SCRIPTS_NFS}/bot_telegram.sh" 2>/dev/null || echo "Telegram bot not available"

module load ISG/apptainer/1.4.0 2>/dev/null || echo "Apptainer module not available"

####==================================####
####          PRINT CONFIG            ####
####==================================####

echo "========================================="
echo "  Viral Contig Concatenation + Merge (per sample)"
echo "  Project:   ${PROJECT_NAME}"
echo "  Ident dir: ${IDENT_DIR}"
echo "  Output:    ${OUTPUT_DIR}"
echo "  Date:      $(date)"
echo "========================================="

tg_send "Starting viral contig concatenation (per-sample)" 2>/dev/null || true

####==================================####
####          SANITY CHECKS           ####
####==================================####

if [[ ! -f "${SEQKIT_CONTAINER}" ]]; then
    echo "ERROR: seqkit container not found at ${SEQKIT_CONTAINER}"
    exit 1
fi

if [[ ! -f "${CAP3_CONTAINER}" ]]; then
    echo "ERROR: CAP3 container not found at ${CAP3_CONTAINER}"
    echo "  Build with: apptainer build cap3.sif docker://quay.io/biocontainers/cap3:10.2015.05.15--h5ef6573_4"
    exit 1
fi

if [[ ! -f "${CDHIT_CONTAINER}" ]]; then
    echo "ERROR: CD-HIT container not found at ${CDHIT_CONTAINER}"
    exit 1
fi

run_seqkit_grep() {
    # $1 = ID file, $2 = FASTA, $3 = output
    local id_file=$1 fasta=$2 out=$3
    local id_dir fasta_dir fasta_base id_file_base
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

####==================================####
####     BUILD SAMPLE LIST            ####
####==================================####

declare -a MAPPING_ARRAY
declare -a READ_TYPE_ARRAY
declare -a SAMPLE_ARRAY

add_sample() {
    local mapping_tool=$1
    local read_type=$2
    local sample=$3

    # One entry per (mapper, read_type, sample). Existence is verified later
    # when the actual per-assembler identification dirs are probed.
    MAPPING_ARRAY+=("$mapping_tool")
    READ_TYPE_ARRAY+=("$read_type")
    SAMPLE_ARRAY+=("$sample")
}

process_mapping_tool() {
    local mapping_tool=$1
    local base_dir=$2
    [[ -d "$base_dir" ]] || { echo "  WARNING: missing $base_dir"; return; }

    echo "  Processing mapping tool: $mapping_tool"

    for read_type_dir in "$base_dir"/*/; do
        [[ -d "$read_type_dir" ]] || continue
        local read_type; read_type=$(basename "$read_type_dir")

        # Use the first assembler present to enumerate samples; the rest
        # are checked per-sample when pooling.
        local probe_asm=""
        for asm in "${ASSEMBLERS_TO_POOL[@]}"; do
            if [[ -d "${read_type_dir}${asm}" ]]; then probe_asm="$asm"; break; fi
        done
        [[ -n "$probe_asm" ]] || continue

        local sample_root="${read_type_dir}${probe_asm}"
        for sample_dir in "$sample_root"/*/; do
            [[ -d "$sample_root" ]] || continue
            add_sample "$mapping_tool" "$read_type" "$(basename "$sample_dir")"
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
READ_TYPE="${READ_TYPE_ARRAY[$JOB_INDEX]}"
SAMPLE="${SAMPLE_ARRAY[$JOB_INDEX]}"

SAMPLE_OUTPUT="${OUTPUT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
mkdir -p "${SAMPLE_OUTPUT}"

echo "========================================="
echo "  Job ${LSB_JOBINDEX}/${TOTAL_JOBS}"
echo "  ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
echo "  Output: ${SAMPLE_OUTPUT}"
echo "========================================="

tg_send "Concatenating (merged): ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true

####========================================================####
####   STEP 1 — gather per-assembler all.fa and pool        ####
####========================================================####

POOLED="${SAMPLE_OUTPUT}/_pooled.fa"
: > "${POOLED}"
pooled_count=0

for asm in "${ASSEMBLERS_TO_POOL[@]}"; do
    asm_ident="${IDENT_DIR}/${MAPPING_TOOL}/${READ_TYPE}/${asm}/${SAMPLE}"
    [[ -d "${asm_ident}" ]] || { echo "  skip ${asm}: no ident dir"; continue; }

    asm_assemblies="${PROJECT_SCRATCH}/results/assembly"
    if [[ "${MAPPING_TOOL}" == "STAR" ]]; then
        asm_base="${asm_assemblies}/star_assembly/${READ_TYPE}/${asm}/${SAMPLE}"
    else
        asm_base="${asm_assemblies}/bowtie_assembly/${READ_TYPE}/${asm}/${SAMPLE}"
    fi

    case $asm in
        "MEGAhit")             assembly_file="${asm_base}/final.contigs.fa" ;;
        "metaSPAdes"|"metaviralSPAdes") assembly_file="${asm_base}/contigs.fasta" ;;
        "rnaSPAdes")
            if [[ -f "${asm_base}/hard_filtered_transcripts.fasta" ]]; then
                assembly_file="${asm_base}/hard_filtered_transcripts.fasta"
            else
                assembly_file="${asm_base}/transcripts.fasta"
            fi
            ;;
    esac

    if [[ ! -s "${assembly_file}" ]]; then
        echo "  skip ${asm}: no assembly file (${assembly_file})"
        continue
    fi

    tmp_asm="${SAMPLE_OUTPUT}/_tmp_${asm}"
    mkdir -p "${tmp_asm}"

    vs2_in="${asm_ident}/virsorter2/final-viral-combined.fa"
    vs2_out="${tmp_asm}/virsorter2.fa"
    : > "${vs2_out}"
    [[ -s "${vs2_in}" ]] && cp "${vs2_in}" "${vs2_out}"

    dvf_in=$(find "${asm_ident}/deepvirfinder" -maxdepth 1 -name '*_dvfpred.txt' 2>/dev/null | head -n1)
    dvf_ids="${tmp_asm}/.dvf_ids.txt"; dvf_out="${tmp_asm}/deepvirfinder.fa"; : > "${dvf_out}"
    if [[ -s "${dvf_in:-}" ]]; then
        awk -v thr="${DVF_PVALUE}" -F'\t' 'NR>1 && $4+0 <= thr {sub(/ .*/,"",$1); print $1}' \
            "${dvf_in}" > "${dvf_ids}"
        [[ -s "${dvf_ids}" ]] && run_seqkit_grep "${dvf_ids}" "${assembly_file}" "${dvf_out}"
    fi

    vitap_in="${asm_ident}/vitap/best_determined_lineages.tsv"
    vitap_ids="${tmp_asm}/.vitap_ids.txt"; vitap_out="${tmp_asm}/vitap.fa"; : > "${vitap_out}"
    if [[ -s "${vitap_in}" ]]; then
        awk -F'\t' -v realms="${VITAP_REALMS}" 'NR>1 && $2 ~ realms {print $1}' \
            "${vitap_in}" | tr -d '\r' \
            | grep -v -E '^[A-Z]{1,3}[0-9]{4,}\.[0-9]+$' > "${vitap_ids}"
        [[ -s "${vitap_ids}" ]] && run_seqkit_grep "${vitap_ids}" "${assembly_file}" "${vitap_out}"
    fi

    asm_all="${tmp_asm}/all.fa"
    python3 - "${vs2_out}" "${vitap_out}" "${dvf_out}" "${asm_all}" <<'PY'
import sys
from pathlib import Path

vs2_path, vitap_path, dvf_path, out_path = sys.argv[1:5]

def base_id(header):
    return header[1:].split("||")[0].split()[0]

def read_fasta(path):
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line, []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)

order = [("DeepVirFinder", dvf_path), ("VITAP", vitap_path), ("VirSorter2", vs2_path)]
merged = {}
for tool, path in order:
    if not Path(path).exists(): continue
    for header, seq in read_fasta(path):
        bid = base_id(header)
        if not bid: continue
        if bid not in merged:
            merged[bid] = {"header": header, "seq": seq, "tools": set()}
        merged[bid]["tools"].add(tool)
        if tool == "VirSorter2":
            merged[bid]["header"] = header
            merged[bid]["seq"] = seq

def wrap(seq, width=80):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

with open(out_path, "w") as fh:
    for bid, rec in merged.items():
        tools = ",".join(sorted(rec["tools"]))
        header = rec["header"].split("||")[0]
        fh.write(f"{header}||tools={tools}\n")
        fh.write(wrap(rec["seq"]) + "\n")
PY

    # Tag headers with assembler and append to pooled
    if [[ -s "${asm_all}" ]]; then
        sed "s/^>/>${asm}|/" "${asm_all}" >> "${POOLED}"
        n=$(grep -c '^>' "${asm_all}")
        echo "  pooled ${asm}: ${n} contigs"
        pooled_count=$((pooled_count + n))
    fi
done

if [[ ! -s "${POOLED}" ]]; then
    echo "  No contigs pooled for ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
    rm -rf "${SAMPLE_OUTPUT}"/_tmp_*
    exit 0
fi

echo "  Total pooled: ${pooled_count} contigs across assemblers"

####==============================================####
####   STEP 2 — CAP3 merges overlapping contigs   ####
####==============================================####

# CAP3 is very sensitive to header characters; sanitize first.
# Keep a map so we can restore original headers later if needed.
SANITIZED="${SAMPLE_OUTPUT}/_pooled_clean.fa"
MAP_FILE="${SAMPLE_OUTPUT}/_header_map.tsv"

awk '
  BEGIN { OFS="\t" }
  /^>/ {
    orig = substr($0, 2)
    clean = orig
    gsub(/\|/, "_", clean)
    gsub(/[^A-Za-z0-9_\.\-]/, "_", clean)
    if (clean == "") clean = "seq" NR
    print clean, orig >> "'"${MAP_FILE}"'"
    print ">" clean
    next
  }
  { print }
' "${POOLED}" > "${SANITIZED}"

n_sanitized=$(grep -c '^>' "${SANITIZED}")
echo "  Sanitized contigs for CAP3: ${n_sanitized}"

# Fresh log each run
: > "${SAMPLE_OUTPUT}/cap3.log"
apptainer exec \
    --bind "${SAMPLE_OUTPUT}":/work \
    "${CAP3_CONTAINER}" \
    bash -c "cd /work && cap3 _pooled_clean.fa -p ${CAP3_P} -o ${CAP3_O} > cap3.log 2>&1"

CAP3_OUT="${SAMPLE_OUTPUT}/_pooled_clean.fa.cap.contigs"
if [[ ! -s "${CAP3_OUT}" ]]; then
    echo "  WARNING: CAP3 produced no output; falling back to sanitized pooled sequences"
    echo "  --- cap3.log (first 40 lines) ---"
    head -40 "${SAMPLE_OUTPUT}/cap3.log" 2>/dev/null || echo "(no cap3.log)"
    cp "${SANITIZED}" "${CAP3_OUT}"
fi

n_cap3=$(awk '/^>/{n++} END{print n+0}' "${CAP3_OUT}" 2>/dev/null)
echo "  CAP3: ${n_cap3} contigs"
if [[ "${n_cap3}" -eq 0 ]]; then
    echo "ERROR: CAP3 output is empty even after fallback"
    ls -la "${SAMPLE_OUTPUT}"
    exit 1
fi

####================================================####
####   STEP 3 — CD-HIT-EST dereplicates residuals   ####
####================================================####

CDHIT_IN="${SAMPLE_OUTPUT}/_cap3_for_cdhit.fa"
awk '
  /^>/ {
    line = substr($0, 2)
    gsub(/\|/, "_", line)
    gsub(/[^A-Za-z0-9_\.\-]/, "_", line)
    print ">" line
    next
  }
  { print }
' "${CAP3_OUT}" > "${CDHIT_IN}"

echo "  Sanitized contigs for CD-HIT: $(grep -c '^>' "${CDHIT_IN}")"

# Remove stale outputs so CD-HIT cannot silently succeed against old files
rm -f "${SAMPLE_OUTPUT}/all.fa" "${SAMPLE_OUTPUT}/all.fa.clstr"
: > "${SAMPLE_OUTPUT}/cdhit.log"

apptainer exec \
    --bind "${SAMPLE_OUTPUT}":/work \
    "${CDHIT_CONTAINER}" \
    bash -c "cd /work && \
        cd-hit-est \
            -i _cap3_for_cdhit.fa \
            -o all.fa \
            -c ${CDHIT_C} -aS ${CDHIT_AS} \
            -T ${CORES:-6} \
            > cdhit.log 2>&1"

if [[ ! -s "${SAMPLE_OUTPUT}/all.fa" ]]; then
    echo "ERROR: CD-HIT-EST failed to produce all.fa"
    echo "  --- cdhit.log (last 60 lines) ---"
    tail -60 "${SAMPLE_OUTPUT}/cdhit.log" 2>/dev/null || echo "(no cdhit.log)"
    echo "  --- CDHIT_IN head ---"
    head -4 "${CDHIT_IN}" 2>/dev/null || echo "(missing)"
    exit 1
fi

n_all=$(awk '/^>/{n++} END{print n+0}' "${SAMPLE_OUTPUT}/all.fa")
echo "  all.fa: ${n_all} non-redundant contigs"

# Keep the cluster mapping for provenance, drop the heavy intermediates
mv "${SAMPLE_OUTPUT}/all.fa.clstr" "${SAMPLE_OUTPUT}/all.clstr" 2>/dev/null || true
rm -f "${SANITIZED}" \
      "${SAMPLE_OUTPUT}/_pooled_clean.fa.cap.singlets" \
      "${SAMPLE_OUTPUT}/_pooled_clean.fa.cap.contigs.links" \
      "${SAMPLE_OUTPUT}/_cap3_for_cdhit.fa" \
      "${POOLED}" "${CAP3_OUT}"
rm -rf "${SAMPLE_OUTPUT}"/_tmp_*

####==================================####
####          COMPLETION              ####
####==================================####

echo "========================================="
echo "  Completed ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}"
echo "  Merged all.fa: ${SAMPLE_OUTPUT}/all.fa"
echo "========================================="

tg_send "Completed merged concatenation: ${MAPPING_TOOL}/${READ_TYPE}/${SAMPLE}" 2>/dev/null || true
