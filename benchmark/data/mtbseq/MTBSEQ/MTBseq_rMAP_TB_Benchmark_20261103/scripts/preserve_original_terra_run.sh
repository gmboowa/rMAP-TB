#!/usr/bin/env bash
# Hash already-downloaded original 18-library archive without copying it.
set -euo pipefail
DIR="${1:-$HOME/Downloads/MTBSEQ/mtbseq_out}"
DIR="$(cd "$DIR" && pwd)"
[[ -f "$DIR/mtbseq_outputs.tar.gz" ]] || { echo 'Missing original archive' >&2; exit 1; }
OUT="$DIR/MTBSEQ_18_LIBRARY_BASELINE_SHA256.txt"
[[ ! -e "$OUT" ]] || { echo "Refusing to overwrite existing checksum: $OUT" >&2; exit 1; }
(
 cd "$DIR"
 shasum -a 256 mtbseq_outputs.tar.gz mtbseq_run.log sample_map.tsv samples.tsv script status.txt
) > "$OUT"
echo "Preserved checksums in $OUT"
echo 'Retain original Terra submission and local archive; this script does not make a second copy.'
