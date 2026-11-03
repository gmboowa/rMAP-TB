"""Emit the per-library benchmarking table in markdown, from primary outputs.

Run:  python3 -I tables.py <raw_dir> <out_dir>
"""
import csv
import os
import re
import sys

RAW, OUT = sys.argv[1], sys.argv[2]


def key(name):
    s = str(name).strip().strip("'").strip('"')
    s = re.sub(r"^rMAPTB_benchmark\.", "", s)
    s = re.sub(r"^\d{4}_", "", s)
    s = re.sub(r"_lib\d+$", "", s)
    return s.replace("-", "_")


def tsv(p):
    with open(p, encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


ORDER = [
    ("MIX_ERR038255_ERR13260175_05pctA_100x", "5 % L4 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_10pctA_100x", "10 % L4 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_20pctA_100x", "20 % L4 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_30pctA_100x", "30 % L4 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_05pctB_100x", "5 % L2 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_10pctB_40x", "10 % L2 minor", "40x"),
    ("MIX_ERR038255_ERR13260175_10pctB_70x", "10 % L2 minor", "70x"),
    ("MIX_ERR038255_ERR13260175_10pctB_100x", "10 % L2 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_10pctB_140x", "10 % L2 minor", "140x"),
    ("MIX_ERR038255_ERR13260175_20pctB_100x", "20 % L2 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_30pctB_100x", "30 % L2 minor", "100x"),
    ("MIX_ERR038255_ERR13260175_50pctB_100x", "50 % L2 minor", "100x"),
    ("MIX_ERR038266_ERR13260106_10pctB_100x", "10 % L4 minor (pair 2)", "100x"),
    ("MIX_ERR038266_ERR13260106_20pctB_100x", "20 % L4 minor (pair 2)", "100x"),
    ("MIX_ERR038266_ERR13260106_50pctB_100x", "50 % L4 minor (pair 2)", "100x"),
    ("PURE_ERR038255_100x", "control: pure L4", "100x"),
    ("PURE_ERR13260175_100x", "control: pure L2", "100x"),
    ("CTRL_SAMELIN_ERR038254_ERR038255_50pct_100x",
     "control: same strain 50:50", "100x"),
]

# rMAP-TB
rm = {}
for r in tsv(os.path.join(RAW, "rmap", "rMAP-TB-V22", "cacheCopy-1",
                          "snp_distance", "clustering_eligibility.tsv")):
    rm[key(r["sample"])] = r["tbprofiler_mixed_infection_status"]

# MAGMA
mg = {}
for r in tsv(os.path.join(RAW, "magma", "collected", "magma_cohort_stats.tsv")):
    lins = [x for x in (r.get("LINEAGES") or "").split(";") if x]
    frs = [x for x in (r.get("FREQUENCIES") or "").split(";") if x]
    mg[key(r["SAMPLE"])] = {
        "mapped": float(r["MAPPED_PERCENTAGE"]),
        "lin": lins, "fr": frs,
        "pass": r["ALL_THRESHOLDS_MET"] == "1",
    }

# MTBseq
mt = {}
for r in tsv(os.path.join(RAW, "mtbseq", "MTBSEQ", "mtbseq_classification.tsv")):
    mt[key(r["SampleID"])] = {
        "lin": (r.get("Coll lineage (branch)") or "").strip("'"),
        "q": (r.get("Coll quality (branch)") or "").strip("'"),
    }
mstat = {}
for r in tsv(os.path.join(RAW, "mtbseq", "MTBSEQ", "mtbseq_statistics.tsv")):
    mstat[key(r["SampleID"])] = float((r["% Mapped Reads"] or "0").strip("'"))

# UShER
uc = {}
for r in tsv(os.path.join(RAW, "usher", "USHER", "call-USHER_DISTANCES",
                          "usher_dist", "usher_clusters.tsv")):
    uc[key(r["sample"])] = r["cluster"]

lines = []
lines.append("| Library (design) | Depth | rMAP-TB v22 | MAGMA | MTBseq | UShER |")
lines.append("|---|---|---|---|---|---|")
for k, design, depth in ORDER:
    # rMAP-TB
    st = rm.get(k, "")
    a = "**flagged**" if st.lower().startswith("possible mixed") else "not flagged"
    # MAGMA
    m = mg.get(k, {})
    if m.get("mapped", 0) == 0.0:
        b = "no reads mapped"
    else:
        nl = len(m.get("lin", []))
        fr = ";".join(m.get("fr", []))
        tag = "**>1 lineage**" if nl > 1 else ("%d lineage" % nl if nl else "none")
        b = "%s (%s)%s" % (tag, fr or "-", "" if m.get("pass") else " *")
    # MTBseq
    t = mt.get(k, {})
    c = "1 lineage (%s)%s" % (t.get("lin", "?"),
                              ", quality bad" if t.get("q") == "bad" else "")
    # UShER
    d = "cluster 1" if uc.get(k) else "-"
    lines.append("| %s<br>%s | %s | %s | %s | %s | %s |"
                 % (k.replace("_", "\\_")[:46], design, depth, a, b, c, d))

lines.append("")
lines.append("\\* did not pass MAGMA's own quality thresholds "
             "(`ALL_THRESHOLDS_MET = 0`).")

out = os.path.join(OUT, "table_per_library.md")
open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("\nwrote", out)
