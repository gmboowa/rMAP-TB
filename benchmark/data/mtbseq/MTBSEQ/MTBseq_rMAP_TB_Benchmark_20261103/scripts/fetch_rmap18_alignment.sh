#!/usr/bin/env bash
# Retrieve the 18-library rMAP-TB core alignment for coordinate-matched comparison.
set -euo pipefail
ROOT='gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/3d10dc72-6796-4f83-85ac-6ebbbf7a816/rMAP_TB/06425efe-8104-42ef-b8f4-7f4b69f4e3b5'
OBJ="$ROOT/call-SNIPPY_CORE_MTBC/cacheCopy/snippy_core/core.full.aln"
OUT="${1:-$HOME/Downloads/MTBseq_rMAP_matched_sites}"
mkdir -p "$OUT"
echo 'First verify the object in this 18-library rMAP-TB submission:'
gsutil ls -lh "$OBJ"
echo 'Downloading the core.full.aln alignment:'
gsutil cp "$OBJ" "$OUT/core.full.aln"
ls -lh "$OUT/core.full.aln"
