#!/bin/bash
# Compress heavy result directories one at a time, verify integrity, then delete.
set -uo pipefail

RESULTS_DIR=/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/results/
cd "$RESULTS_DIR" || { echo "ERROR: cannot cd to $RESULTS_DIR"; exit 1; }

# Whole directories to compress, biggest first
TARGETS=(
    "aligned"
    "assembly"
    "trimmed"
    "trimmed_qc"
    "untrimmed_qc"
)

for base in "${TARGETS[@]}"; do
    src="$base"
    archive="${base}.tar.zst"
    tmp="${archive}.tmp"

    echo "=========================================="
    echo "  $(date)  Compressing $src"
    echo "=========================================="

    if [[ ! -d "$src" ]]; then
        echo "  SKIP: $src not found (maybe already done)"
        continue
    fi
    if [[ -f "$archive" ]]; then
        echo "  SKIP: $archive already exists"
        continue
    fi

    # Headroom check on the filesystem holding the archive
    avail_kb=$(df -Pk . | awk 'NR==2 {print $4}')
    needed_kb=$(du -sk "$src" | awk '{print $1}')
    echo "  Original size: $(du -sh "$src" | cut -f1)"
    echo "  Free on filesystem: $((avail_kb/1024/1024))G"
    if [[ "$avail_kb" -lt "$needed_kb" ]]; then
        echo "  ABORT: not enough space for temp archive (need ~$((needed_kb/1024/1024))G)"
        continue
    fi

    # Compress
    if ! tar -I 'zstd -T8 -3' -cf "$tmp" "$base"; then
        echo "  ERROR: tar failed for $src"
        rm -f "$tmp"
        continue
    fi

    # Verify integrity
    if ! zstd -t "$tmp" >/dev/null 2>&1; then
        echo "  ERROR: integrity test failed for $tmp"
        rm -f "$tmp"
        continue
    fi

    # Verify file count matches
    orig_count=$(find "$src" -type f | wc -l)
    arch_count=$(tar -tf "$tmp" | grep -v '/$' | wc -l)
    if [[ "$orig_count" -ne "$arch_count" ]]; then
        echo "  ERROR: file count mismatch ($orig_count in dir vs $arch_count in archive)"
        rm -f "$tmp"
        continue
    fi

    # Promote and delete
    mv "$tmp" "$archive"
    echo "  Archive: $(du -sh "$archive" | cut -f1)  (was $(du -sh "$src" | cut -f1))"
    rm -rf "$src"
    echo "  Deleted $src/"
    echo "  Group quota now:"
    lfs quota -h -g team222 /lustre/scratch126 2>/dev/null | tail -1
done

echo "=========================================="
echo "  All done at $(date)"
echo "=========================================="
du -sh "$RESULTS_DIR"
lfs quota -h -g team222 /lustre/scratch126
