#!/usr/bin/env bash
# Recompute every benchmarking value reported in the rMAP-TB manuscript from
# the primary output files of the four workflows. Run from this directory.
#
#   bash run_all.sh [data_dir] [out_dir]
#
# Defaults: data/ and results/. Python 3.9+; matplotlib only for the figure.
set -euo pipefail

DATA="${1:-data}"
OUT="${2:-results}"

if [ ! -d "$DATA" ]; then
    echo "ERROR: no such directory: $DATA" >&2
    exit 1
fi
mkdir -p "$OUT"

echo "==> per-library comparison across the four workflows"
python3 -I scripts/compare.py "$DATA" "$OUT"

echo "==> pair-level analysis and the three unmixed control pairs"
python3 -I scripts/pairs.py "$DATA" "$OUT"

echo "==> per-library table"
python3 -I scripts/tables.py "$DATA" "$OUT"

if python3 -c "import matplotlib" 2>/dev/null; then
    echo "==> Figure 6"
    python3 -I scripts/fig_bench.py "$OUT"
else
    echo "==> Figure 6 skipped (matplotlib not installed)"
fi

echo
echo "Done. Results in $OUT/:"
ls -1 "$OUT"
