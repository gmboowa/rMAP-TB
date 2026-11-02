#!/usr/bin/env python3
"""Validate an rMAP_TB inputs JSON before submitting it to Cromwell or Terra.

Checks performed:
  * every key belongs to the rMAP_TB namespace and is a declared workflow input
    (declared inputs are parsed from the WDL, so the check follows the workflow version)
  * required inputs are present
  * read1 and read2 have the same length, no duplicates, and each pair shares a sample stem
  * local files exist and are non-empty (gs:// and other remote URIs are reported, not opened)
  * gzip FASTQ files start with the gzip magic bytes
  * booleans/numbers have the expected JSON type
  * phylogeny is only possible when at least min_mtbc_samples_for_tree pairs are supplied

Exit status is 0 when no errors are found (warnings do not fail the check).

Example:
    python3 scripts/validate_inputs.py inputs.json --wdl rMAP_TB.wdl
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

WORKFLOW = "rMAP_TB"
REQUIRED = ["read1", "read2", "adapters", "mtbc_reference_genbank"]
MATE_RE = re.compile(r"([._]R?[12](?:_001)?|_R[12]_paired)(\.(fastq|fq)(\.gz)?)$", re.IGNORECASE)
REMOTE_RE = re.compile(r"^(gs|s3|https?|drs)://")


def parse_wdl_inputs(wdl_path: Path) -> dict[str, tuple[str, bool]]:
    """Return {name: (wdl_type, is_optional)} for the workflow input block."""
    text = wdl_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(rf"workflow\s+{WORKFLOW}\s*\{{\s*input\s*\{{(.*?)\n\s*\}}", text, re.S)
    if not m:
        sys.exit(f"error: could not find the {WORKFLOW} input block in {wdl_path}")
    declared = {}
    for line in m.group(1).splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        dm = re.match(r"([A-Za-z\[\]\?\+ ]+?)\s+([A-Za-z_][A-Za-z0-9_]*)\s*(=.*)?$", line)
        if not dm:
            continue
        wdl_type, name, default = dm.group(1).strip(), dm.group(2), dm.group(3)
        optional = wdl_type.endswith("?") or default is not None
        declared[name] = (wdl_type, optional)
    return declared


def stem(path: str) -> str:
    name = path.rsplit("/", 1)[-1]
    return MATE_RE.sub("", name)


def check_file(path: str, label: str, errors: list, warnings: list, fastq: bool = False):
    if REMOTE_RE.match(path):
        warnings.append(f"{label}: remote path not checked ({path})")
        return
    p = Path(path)
    if not p.is_absolute():
        warnings.append(f"{label}: relative path; Cromwell resolves relative to its working directory ({path})")
    if not p.is_file():
        errors.append(f"{label}: file not found ({path})")
        return
    if p.stat().st_size == 0:
        errors.append(f"{label}: file is empty ({path})")
        return
    if fastq and path.lower().endswith(".gz"):
        with p.open("rb") as fh:
            if fh.read(2) != b"\x1f\x8b":
                errors.append(f"{label}: .gz file does not start with gzip magic bytes ({path})")


def expected_python_type(wdl_type: str):
    base = wdl_type.rstrip("?+").strip()
    if base == "Boolean":
        return (bool,)
    if base == "Int":
        return (int,)
    if base == "Float":
        return (int, float)
    if base == "String" or base == "File":
        return (str,)
    if base.startswith("Array"):
        return (list,)
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", type=Path)
    ap.add_argument("--wdl", type=Path, default=Path(__file__).resolve().parents[1] / "rMAP_TB.wdl")
    a = ap.parse_args(argv)

    errors: list[str] = []
    warnings: list[str] = []

    try:
        data = json.loads(a.inputs.read_text())
    except Exception as exc:
        sys.exit(f"error: cannot parse {a.inputs}: {exc}")
    if not isinstance(data, dict):
        sys.exit("error: inputs JSON must be an object")

    declared = parse_wdl_inputs(a.wdl)

    # Namespace and unknown keys
    values = {}
    for key, value in data.items():
        if not key.startswith(f"{WORKFLOW}."):
            errors.append(f"key {key!r} is not in the {WORKFLOW}. namespace")
            continue
        name = key[len(WORKFLOW) + 1:]
        if "." in name:
            warnings.append(f"key {key!r} targets a task-level input; only workflow-level inputs are validated")
            continue
        if name not in declared:
            errors.append(f"unknown workflow input {name!r}")
            continue
        wdl_type, _ = declared[name]
        exp = expected_python_type(wdl_type)
        if exp and not isinstance(value, exp) or (exp == (int,) and isinstance(value, bool)):
            errors.append(f"{name}: expected {wdl_type}, got {type(value).__name__}")
        values[name] = value

    for name in REQUIRED:
        if name not in values:
            errors.append(f"required input {name!r} is missing")

    r1 = values.get("read1", []) or []
    r2 = values.get("read2", []) or []
    if isinstance(r1, list) and isinstance(r2, list):
        if len(r1) != len(r2):
            errors.append(f"read1 has {len(r1)} files but read2 has {len(r2)}")
        if len(set(r1)) != len(r1) or len(set(r2)) != len(r2):
            errors.append("duplicate FASTQ paths in read1/read2")
        for i, (x, y) in enumerate(zip(r1, r2)):
            if not isinstance(x, str) or not isinstance(y, str):
                errors.append(f"read pair {i}: paths must be strings")
                continue
            if stem(x) != stem(y):
                warnings.append(f"read pair {i}: stems differ ({x.rsplit('/',1)[-1]} vs {y.rsplit('/',1)[-1]})")
            check_file(x, f"read1[{i}]", errors, warnings, fastq=True)
            check_file(y, f"read2[{i}]", errors, warnings, fastq=True)
        min_pairs = values.get("min_mtbc_samples_for_tree", 3)
        if 0 < len(r1) < min_pairs:
            warnings.append(f"only {len(r1)} pair(s) supplied; the phylogenomics branch needs at least {min_pairs} MTBC pairs")

    for name in ("adapters", "mtbc_reference_genbank", "snp_mask_bed", "sample_metadata_tsv"):
        if isinstance(values.get(name), str):
            check_file(values[name], name, errors, warnings)
    if "snp_mask_bed" not in values:
        warnings.append("snp_mask_bed not set; SNP distances will include PE/PPE regions and are not comparable with published thresholds")

    for w in warnings:
        print(f"WARNING  {w}")
    for e in errors:
        print(f"ERROR    {e}")
    print(f"{len(errors)} error(s), {len(warnings)} warning(s); {len(r1) if isinstance(r1, list) else 0} read pair(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
