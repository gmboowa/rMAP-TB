#!/usr/bin/env python3
"""Copy an integrated report into docs/reports/ under the repository naming convention.

Name pattern:
    rMAP-TB_<version>_run<UTC stamp>_<n>-sample_<panel>_integrated_report.html

The run stamp and sample count are read from the report itself, so the file
name always agrees with its content. The panel label is supplied by the user
(for example `in-silico-mixture-panel` or `NTM-contamination-panel`). An entry
is appended to docs/reports/README.md with the headline figures.

Example:
    python3 scripts/archive_report.py \
        results/<id>/report/integrated_tb_amr_mtbc_phylogenomics_report.html \
        --version v22 --panel NTM-contamination-panel \
        --description "12 MTBC isolates, one isolate spiked with M. avium reads at 10 % and 50 %, six NTM species"
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize_report import summarize  # noqa: E402

README_HEADER = """# Example integrated reports

Each file is the complete, self-contained HTML report produced by `MERGE_TB_REPORTS`
for one validation panel. Open it in any browser; no server or network access is
required. File names follow
`rMAP-TB_<version>_run<UTC stamp>_<n>-sample_<panel>_integrated_report.html`,
where the run stamp is the one printed in the report header.

| File | Workflow | Commit | Samples (total / MTBC / non-MTBC) | Drug-resistant | In tree | Panel |
|---|---|---|---|---|---|---|
"""


def git_commit(repo: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "n/a"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", type=Path)
    ap.add_argument("--version", required=True, help="Workflow version label, e.g. v22")
    ap.add_argument("--panel", required=True, help="Short panel label, letters/digits/hyphens only")
    ap.add_argument("--description", default="", help="One-line description for docs/reports/README.md")
    ap.add_argument("--dest", type=Path, default=Path(__file__).resolve().parents[1] / "docs" / "reports")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    if not a.report.is_file():
        sys.exit(f"error: {a.report} not found")
    if not re.fullmatch(r"[A-Za-z0-9-]+", a.panel):
        sys.exit("error: --panel may contain only letters, digits and hyphens")

    s = summarize(a.report)
    stamp = s["run_stamp"].replace("_UTC", "").replace("_", "-")
    if not re.fullmatch(r"\d{8}-\d{6}", stamp):
        sys.exit(f"error: could not read a run stamp from the report (got {s['run_stamp']!r})")
    n = s["counts"].get("total_paired_samples")
    if n is None:
        sys.exit("error: could not read the sample count from the report")

    name = f"rMAP-TB_{a.version}_run{stamp}_{n}-sample_{a.panel}_integrated_report.html"
    target = a.dest / name
    repo = Path(__file__).resolve().parents[1]
    commit = git_commit(repo)

    c = s["counts"]
    row = (
        f"| `{name}` | {a.version} | {commit} | "
        f"{c.get('total_paired_samples', '?')} / {c.get('mtbc_isolates_retained', '?')} / {c.get('non-mtbc_excluded', c.get('non_mtbc_excluded', '?'))} | "
        f"{c.get('drug-resistant_isolates', c.get('drug_resistant_isolates', '?'))} | {s['samples_in_tree']} | "
        f"{a.description or a.panel} |\n"
    )

    if a.dry_run:
        print(f"would copy {a.report} -> {target}")
        print(row, end="")
        return 0

    a.dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(a.report, target)
    readme = a.dest / "README.md"
    if not readme.exists():
        readme.write_text(README_HEADER)
    content = readme.read_text()
    if name not in content:
        with readme.open("a") as fh:
            fh.write(row)
    print(f"archived {target}")
    if s["tree_tips_unannotated"]:
        print(f"WARNING: {s['tree_tips_unannotated']} tree tip(s) have no lineage/resistance annotation in this report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
