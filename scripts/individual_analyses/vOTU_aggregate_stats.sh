#!/usr/bin/env bash
# Aggregate per-sample vOTU stats into project-level tables.
# Run once, after the vOTU_stats array has finished.

set -o pipefail

PROJECT_NAME="${1:-mosquito_virome_yucatan_LEVE}"
SCRATCH_BASE="/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects"
PERMANENT_BASE="/nfs/users/nfs_j/jr46"

STATS_DIR="${SCRATCH_BASE}/${PROJECT_NAME}/docs/vOTU_stats"
NFS_STATS_DIR="${PERMANENT_BASE}/git_repos/${PROJECT_NAME}/docs/vOTU_stats"

OUT_MASTER="${STATS_DIR}/vOTU_master_stats_ALL_SAMPLES.tsv"
OUT_SUMMARY="${STATS_DIR}/vOTU_per_sample_summary.tsv"

echo "Aggregating from: ${STATS_DIR}"

# ---- Pick a reference master file, capture its header ----
first_master=$(find "$STATS_DIR" -name "vOTU_master_stats.tsv" | sort | head -1)
if [[ -z "$first_master" ]]; then
    echo "ERROR: no per-sample master files found"
    exit 1
fi
master_header=$(head -1 "$first_master")

first_summary=$(find "$STATS_DIR" -name "vOTUs_stats_overall.tsv" | sort | head -1)
if [[ -z "$first_summary" ]]; then
    echo "ERROR: no per-sample summary files found"
    exit 1
fi
summary_header=$(head -1 "$first_summary")

# ---- Master table: per-OTU, all samples ----
printf 'mapping_tool\tread_type\tsample\t%s\n' "$master_header" > "$OUT_MASTER"

n_master=0
skipped_master=0
while IFS= read -r f; do
    if [[ "$(head -1 "$f")" != "$master_header" ]]; then
        echo "  SKIP (header mismatch): $f"
        skipped_master=$((skipped_master + 1))
        continue
    fi
    rel="${f#${STATS_DIR}/}"
    map="${rel%%/*}"; rest="${rel#*/}"
    rtype="${rest%%/*}"; rest2="${rest#*/}"
    samp="${rest2%%/*}"

    tail -n +2 "$f" | awk -F'\t' -v OFS='\t' -v m="$map" -v rt="$rtype" -v s="$samp" \
        '{print m, rt, s, $0}'
    n_master=$((n_master + 1))
done < <(find "$STATS_DIR" -name "vOTU_master_stats.tsv" | sort) >> "$OUT_MASTER"

echo "Master table: $OUT_MASTER"
echo "  Files merged: $n_master  (skipped: $skipped_master)"
echo "  OTU rows: $(( $(wc -l < "$OUT_MASTER") - 1 ))"

# ---- Per-sample summary ----
printf 'mapping_tool\tread_type\tsample\t%s\n' "$summary_header" > "$OUT_SUMMARY"

n_sum=0
while IFS= read -r f; do
    rel="${f#${STATS_DIR}/}"
    map="${rel%%/*}"; rest="${rel#*/}"
    rtype="${rest%%/*}"; rest2="${rest#*/}"
    samp="${rest2%%/*}"

    tail -n +2 "$f" | awk -F'\t' -v OFS='\t' -v m="$map" -v rt="$rtype" -v s="$samp" \
        '{print m, rt, s, $0}'
    n_sum=$((n_sum + 1))
done < <(find "$STATS_DIR" -name "vOTUs_stats_overall.tsv" | sort) >> "$OUT_SUMMARY"

echo "Summary table: $OUT_SUMMARY"
echo "  Samples: $n_sum"

# ---- Copy to NFS ----
mkdir -p "$NFS_STATS_DIR"
rsync -a "$OUT_MASTER" "$OUT_SUMMARY" "$NFS_STATS_DIR/"
echo "Copied to NFS: $NFS_STATS_DIR"
