#!/usr/bin/env python3
"""Extract the headline results from an rMAP-TB integrated HTML report.

The integrated report is self-contained HTML. This script reads it without a
browser and prints (or writes as JSON/TSV) the figures a reader would otherwise
have to look up by hand: run stamp, sample counts, routing, resistance-category
distribution, lineage counts, number of genomically close pairs, IQ-TREE and
tree-rendering status, and the per-sample tree annotation. It is used by
`archive_report.py` to name files and by the test-suite to verify reports.

Example:
    python3 scripts/summarize_report.py integrated_tb_amr_mtbc_phylogenomics_report.html
    python3 scripts/summarize_report.py report.html --json summary.json
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from pathlib import Path


def text(s: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", " ", s)).replace("\xa0", " ").strip()


def first(pattern: str, s: str, flags=0, default=""):
    m = re.search(pattern, s, flags)
    return m.group(1).strip() if m else default


def summarize(path: Path) -> dict:
    s = path.read_text(encoding="utf-8", errors="replace")
    out: dict = {"file": str(path), "bytes": len(s.encode("utf-8", errors="replace"))}

    out["run_stamp"] = first(r"Run stamp:</strong>\s*<span class=\"mono\">([^<]+)", s) or first(r"Run stamp:</strong>\s*([^<&]+)", s)
    out["run_generated"] = first(r"Run generated:</strong>\s*([^<]+)", s)

    cards = re.findall(r'<div class="card [a-z]+"><h3>([^<]+)</h3><div class="num">(\d+)</div>', s)
    out["counts"] = {k.strip().lower().replace(" ", "_"): int(v) for k, v in cards}

    out["status"] = {
        k.strip(): v.strip()
        for k, v in re.findall(r'<div class="pill [^"]*"><span class="pill-dot"></span><div><strong>([^<]+)</strong><small>([^<]*)</small>', s)
    }

    out["resistance_distribution"] = {
        text(k): int(v)
        for k, v in re.findall(r'<span class="dist-chip"><i[^>]*></i>([^<]+) <b>(\d+)</b>', s)
    }

    out["lineage_distribution"] = {
        text(k): int(v)
        for k, v in re.findall(r'<div class="lineage-label">([^<]+)</div>.*?<div class="lineage-value">(\d+)', s, re.S)
    }

    out["iqtree_status"] = first(r"IQ-TREE status:</strong>\s*([^<&]+)", s)
    out["samples_in_tree"] = int(first(r"Samples retained for tree:</strong>\s*(\d+)", s, default="0") or 0)

    chip = re.search(r'class="state-chip [^"]*">([^<]+)</span><span class="muted">([^<]*)', s)
    out["tree_rendering"] = {"state": chip.group(1).strip(), "detail": html.unescape(chip.group(2)).strip()} if chip else {}

    tree_block = ""
    if 'id="treeFigure"' in s:
        tree_block = s.split('id="treeFigure"', 1)[1]
        tree_block = tree_block.split("</svg>", 1)[0]
    tips = []
    if tree_block:
        for m in re.finditer(r'<circle[^>]*><title>([^<]+)</title></circle>', tree_block):
            parts = [p.strip() for p in html.unescape(m.group(1)).split(" — ")]
            if len(parts) == 3:
                tips.append({"sample": parts[0], "lineage": parts[1], "resistance": parts[2]})
    out["tree_tips"] = tips
    out["tree_tips_unannotated"] = sum(1 for t in tips if t["lineage"] == "Lineage NA" and t["resistance"] == "Unknown")

    out["close_pairs"] = len(re.findall(r"<tr><td>Cluster_\d+</td>", s)) or (s.count("Genomically close; review epidemiological linkage") // 2)
    out["mixed_infection_flags"] = len(re.findall(r">Possible mixed infection<", s))
    out["unparsed_inputs"] = "Report inputs that could not be parsed" in s

    excluded = re.findall(r"<li><strong>([^<]+)</strong>\s*(Excluded before IQ-TREE[^<]{0,200}|was excluded from IQ-TREE[^<]{0,200})", s)
    out["excluded_samples"] = [{"sample": text(a), "reason": text(b)} for a, b in excluded]

    species = re.findall(r"<tr>\s*<td class='mono'>([^<]+)</td>\s*<td><span class=\"badge\"[^>]*>(?:<em>)?([^<]+)", s)
    out["species_calls"] = {text(k): text(v) for k, v in species}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", type=Path)
    ap.add_argument("--json", type=Path, help="Write the full summary as JSON")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    if not a.report.is_file():
        sys.exit(f"error: {a.report} not found")
    summary = summarize(a.report)
    if a.json:
        a.json.write_text(json.dumps(summary, indent=2) + "\n")
    if not a.quiet:
        c = summary["counts"]
        print(f"Report:            {a.report.name}")
        print(f"Run stamp:         {summary['run_stamp']}")
        print(f"Samples:           {c.get('total_paired_samples', '?')} total, {c.get('mtbc_isolates_retained', '?')} MTBC, {c.get('non-mtbc_excluded', c.get('non_mtbc_excluded', '?'))} non-MTBC, {c.get('drug-resistant_isolates', c.get('drug_resistant_isolates', '?'))} drug-resistant")
        print(f"IQ-TREE:           {summary['iqtree_status']} ({summary['samples_in_tree']} in tree)")
        tr = summary["tree_rendering"]
        print(f"Tree rendering:    {tr.get('state', '?')} - {tr.get('detail', '')}")
        print(f"Tree tips:         {len(summary['tree_tips'])} drawn, {summary['tree_tips_unannotated']} without annotation")
        print(f"Resistance:        {summary['resistance_distribution']}")
        print(f"Lineages:          {summary['lineage_distribution']}")
        print(f"Close pairs:       {summary['close_pairs']}   Mixed-infection flags: {summary['mixed_infection_flags']}")
        if summary["excluded_samples"]:
            print("Excluded samples:  " + "; ".join(f"{e['sample']}: {e['reason'][:80]}" for e in summary["excluded_samples"]))
        if summary["unparsed_inputs"]:
            print("WARNING: the report lists inputs that could not be parsed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
