#!/usr/bin/env python3
"""Build a Cromwell inputs JSON for rMAP_TB from a directory of paired FASTQ files.

Pairs are detected from common Illumina naming conventions
(_R1/_R2, _1/_2, .R1/.R2, with optional _001 suffixes) and sorted by sample
name so that read1[i] and read2[i] always refer to the same sample. Absolute
paths are written, as Cromwell requires.

Example:
    python3 scripts/make_inputs_json.py \
        --fastq-dir /data/run42/fastq \
        --adapters resources/adapters/adapters.fa \
        --reference resources/reference/H37Rv_NC_000962.3.gb \
        --mask resources/mask/H37Rv_NC_000962.3_PE_PPE_repetitive.bed \
        --metadata metadata.tsv \
        --output inputs.json
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

FASTQ_EXT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)
R1_PATTERNS = [
    re.compile(r"^(?P<sample>.+?)[._]R1(?:_001)?$"),
    re.compile(r"^(?P<sample>.+?)[._]1$"),
    re.compile(r"^(?P<sample>.+?)_R1_paired$"),
]
R2_PATTERNS = [
    re.compile(r"^(?P<sample>.+?)[._]R2(?:_001)?$"),
    re.compile(r"^(?P<sample>.+?)[._]2$"),
    re.compile(r"^(?P<sample>.+?)_R2_paired$"),
]

WORKFLOW = "rMAP_TB"


def classify(path: Path):
    stem = FASTQ_EXT.sub("", path.name)
    for pat in R1_PATTERNS:
        m = pat.match(stem)
        if m:
            return m.group("sample"), 1
    for pat in R2_PATTERNS:
        m = pat.match(stem)
        if m:
            return m.group("sample"), 2
    return None, None


def find_pairs(fastq_dir: Path, recursive: bool):
    files = sorted(fastq_dir.rglob("*") if recursive else fastq_dir.iterdir())
    pairs: dict[str, dict[int, Path]] = {}
    unmatched = []
    for f in files:
        if not f.is_file() or not FASTQ_EXT.search(f.name):
            continue
        sample, mate = classify(f)
        if sample is None:
            unmatched.append(f)
            continue
        pairs.setdefault(sample, {})[mate] = f.resolve()
    complete, incomplete = {}, {}
    for sample, mates in pairs.items():
        (complete if 1 in mates and 2 in mates else incomplete)[sample] = mates
    return complete, incomplete, unmatched


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastq-dir", required=True, type=Path, help="Directory containing paired FASTQ files")
    p.add_argument("--recursive", action="store_true", help="Search sub-directories as well")
    p.add_argument("--adapters", required=True, type=Path, help="Adapter FASTA for Trimmomatic")
    p.add_argument("--reference", required=True, type=Path, help="MTBC reference in GenBank format")
    p.add_argument("--mask", type=Path, help="Optional SNP mask BED (strongly recommended)")
    p.add_argument("--metadata", type=Path, help="Optional sample metadata TSV")
    p.add_argument("--set", action="append", default=[], metavar="PARAM=VALUE",
                   help="Additional workflow parameter, e.g. --set use_gubbins=true --set max_cpus=8")
    p.add_argument("--output", type=Path, default=Path("inputs.json"))
    p.add_argument("--allow-incomplete", action="store_true",
                   help="Write the file even if some samples lack a mate (they are skipped)")
    a = p.parse_args(argv)

    if not a.fastq_dir.is_dir():
        sys.exit(f"error: {a.fastq_dir} is not a directory")

    complete, incomplete, unmatched = find_pairs(a.fastq_dir, a.recursive)

    for f in unmatched:
        print(f"warning: could not determine mate for {f.name}; skipped", file=sys.stderr)
    for sample, mates in incomplete.items():
        print(f"warning: sample {sample} has only R{list(mates)[0]}; skipped", file=sys.stderr)
    if incomplete and not a.allow_incomplete:
        sys.exit("error: incomplete pairs found; fix the files or pass --allow-incomplete")
    if not complete:
        sys.exit("error: no complete read pairs found")

    for label, path in (("adapters", a.adapters), ("reference", a.reference), ("mask", a.mask), ("metadata", a.metadata)):
        if path is not None and not path.is_file():
            sys.exit(f"error: {label} file not found: {path}")

    samples = sorted(complete)
    inputs = {
        f"{WORKFLOW}.read1": [str(complete[s][1]) for s in samples],
        f"{WORKFLOW}.read2": [str(complete[s][2]) for s in samples],
        f"{WORKFLOW}.adapters": str(a.adapters.resolve()),
        f"{WORKFLOW}.mtbc_reference_genbank": str(a.reference.resolve()),
    }
    if a.mask:
        inputs[f"{WORKFLOW}.snp_mask_bed"] = str(a.mask.resolve())
    if a.metadata:
        inputs[f"{WORKFLOW}.sample_metadata_tsv"] = str(a.metadata.resolve())

    for item in a.set:
        if "=" not in item:
            sys.exit(f"error: --set expects PARAM=VALUE, got {item!r}")
        key, raw = item.split("=", 1)
        key = key if key.startswith(f"{WORKFLOW}.") else f"{WORKFLOW}.{key}"
        low = raw.strip().lower()
        if low in {"true", "false"}:
            value = low == "true"
        else:
            try:
                value = int(raw)
            except ValueError:
                try:
                    value = float(raw)
                except ValueError:
                    value = raw
        inputs[key] = value

    a.output.write_text(json.dumps(inputs, indent=2) + "\n")
    print(f"wrote {a.output} with {len(samples)} sample pair(s): {', '.join(samples)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
