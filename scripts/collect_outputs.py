#!/usr/bin/env python3
"""Gather the principal rMAP_TB outputs of a Cromwell run into one directory.

Cromwell leaves outputs scattered under cromwell-executions/<workflow>/<id>/call-*/.
This script reads the outputs JSON that `cromwell run` prints (save it with
`cromwell run ... | tee run.log`, or use `--metadata` from the Cromwell server API)
and copies or symlinks every file output into a flat, named results directory:

    results/<run-id>/
        report/     integrated HTML report and run metadata
        tables/     TSV summaries (TB-Profiler, routing, SNP distances, metadata, ...)
        figures/    tree SVG/PNG, heatmap and lineage SVGs
        alignments/ core alignments, VCF, Newick tree
        qc/         MultiQC, FastQC summaries, trimming summary
        logs/       status and log files
        manifest.tsv  output identifier -> original path -> copied path

Example:
    java -jar cromwell.jar run rMAP_TB.wdl --inputs inputs.json | tee run.log
    python3 scripts/collect_outputs.py --cromwell-log run.log --dest results/
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

CATEGORY_RULES = [
    ("report", re.compile(r"final_merged_html_report|run_metadata")),
    ("figures", re.compile(r"tree_svg|tree_image|_svg$")),
    ("alignments", re.compile(r"alignment|vcf|newick|core_|gubbins_filtered")),
    ("qc", re.compile(r"multiqc|fastqc|trimming|mean_depth|coverage_qc")),
    ("logs", re.compile(r"status|_log$|render_log")),
    ("tables", re.compile(r".*")),
]


def categorize(name: str) -> str:
    short = name.split(".")[-1]
    for cat, pat in CATEGORY_RULES:
        if pat.search(short):
            return cat
    return "tables"


def load_outputs(args) -> tuple[dict, str]:
    if args.outputs_json:
        data = json.loads(Path(args.outputs_json).read_text())
        outputs = data.get("outputs", data)
        run_id = data.get("id", "run")
        return outputs, run_id
    if args.metadata:
        data = json.loads(Path(args.metadata).read_text())
        return data.get("outputs", {}), data.get("id", "run")
    text = Path(args.cromwell_log).read_text(errors="replace")
    # `cromwell run` prints a JSON object with "outputs" and "id"; it may be the
    # last JSON object in the log. Find the last balanced object containing "outputs".
    candidates = [m.start() for m in re.finditer(r"\{\s*\n?\s*\"outputs\"", text)]
    for start in reversed(candidates):
        depth = 0
        for i in range(start, len(text)):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
                if depth == 0:
                    try:
                        data = json.loads(text[start:i + 1])
                        return data.get("outputs", {}), data.get("id", "run")
                    except json.JSONDecodeError:
                        break
    sys.exit("error: could not find a Cromwell outputs JSON object in the log")


def iter_files(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for v in value:
            yield from iter_files(v)
    elif isinstance(value, dict):
        for v in value.values():
            yield from iter_files(v)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--cromwell-log", help="Log captured from `cromwell run` (contains the outputs JSON)")
    src.add_argument("--outputs-json", help="JSON with an 'outputs' object (e.g. from the /outputs API)")
    src.add_argument("--metadata", help="Cromwell metadata JSON (from the /metadata API)")
    ap.add_argument("--dest", type=Path, default=Path("results"))
    ap.add_argument("--symlink", action="store_true", help="Symlink instead of copying")
    ap.add_argument("--include-arrays", action="store_true",
                    help="Also collect per-sample array outputs (FastQC reports, Kraken reports, VCFs, ...)")
    a = ap.parse_args(argv)

    outputs, run_id = load_outputs(a)
    if not outputs:
        sys.exit("error: outputs object is empty")

    dest = a.dest / run_id
    dest.mkdir(parents=True, exist_ok=True)
    manifest = []
    missing = 0

    for name, value in outputs.items():
        if value is None:
            continue
        if isinstance(value, list) and not a.include_arrays:
            continue
        cat = categorize(name)
        for path in iter_files(value):
            if re.match(r"^(gs|s3|https?)://", path):
                manifest.append((name, path, "remote; not copied"))
                continue
            p = Path(path)
            if not p.exists():
                missing += 1
                manifest.append((name, path, "MISSING"))
                continue
            outdir = dest / cat
            outdir.mkdir(exist_ok=True)
            target = outdir / p.name
            if target.exists() and target.resolve() != p.resolve():
                target = outdir / f"{name.split('.')[-1]}__{p.name}"
            if a.symlink:
                if target.exists() or target.is_symlink():
                    target.unlink()
                target.symlink_to(p.resolve())
            else:
                shutil.copy2(p, target)
            manifest.append((name, str(p), str(target)))

    with (dest / "manifest.tsv").open("w") as fh:
        fh.write("output\tsource\tdestination\n")
        for row in manifest:
            fh.write("\t".join(row) + "\n")

    print(f"collected {len(manifest)} file(s) into {dest} ({missing} missing)")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
