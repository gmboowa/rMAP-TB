#!/usr/bin/env bash
# Run MTBseq 1.1.0 TBjoin->TBamend->TBgroups on two already-processed pure controls.
# Does not realign reads, rerun GATK, or change baseline 18-sample files.
set -euo pipefail
PURE_DIR="${1:-$HOME/Downloads/MTBseq_pure_only}"
PURE_DIR="$(cd "$PURE_DIR" && pwd)"
IMAGE='quay.io/biocontainers/mtbseq:1.1.0--hdfd78af_1'
PROJECT='rMAPTB_pure2'
if [[ ! -f "$PURE_DIR/samples.tsv" ]]; then echo 'Missing samples.tsv; first run prepare_pure_only_from_tar.py' >&2; exit 1; fi
for s in PURE-ERR038255-100x PURE-ERR13260175-100x; do
  compgen -G "$PURE_DIR/Called/$s"'*gatk_position_variants*.tab' >/dev/null || { echo "Missing Called table for $s" >&2; exit 1; }
  compgen -G "$PURE_DIR/Position_Tables/$s"'*gatk_position_table.tab' >/dev/null || { echo "Missing Position_Tables table for $s" >&2; exit 1; }
done
command -v docker >/dev/null || { echo 'Docker CLI unavailable. Start Colima (if installed) and check docker version.' >&2; exit 1; }
docker info >/dev/null || { echo 'Docker daemon unavailable. Start Colima (colima start), then retry.' >&2; exit 1; }
mkdir -p "$PURE_DIR/logs"
common=( --project "$PROJECT" --samples samples.tsv --ref M._tuberculosis_H37Rv_2015-11-13 --mincovf 4 --mincovr 4 --minphred 4 --minfreq 75 --unambig 95 --window 12 --distance 12 )
run_step() {
  local STEP="$1"
  echo "Starting $STEP at $(date -u +%FT%TZ)"
  docker run --rm --platform linux/amd64 -v "$PURE_DIR:/work" -w /work "$IMAGE" \
    MTBseq --step "$STEP" "${common[@]}" 2>&1 | tee "$PURE_DIR/logs/${STEP}.log"
}
run_step TBjoin
run_step TBamend
run_step TBgroups
find "$PURE_DIR/Amend" "$PURE_DIR/Groups" -maxdepth 1 -type f | sort
printf '\nDONE: pure-control-only native MTBseq joint analysis generated; the 18-sample baseline was not changed.\n'
