"""Pair-level analysis: what each tool reported for the 153 pairs, and for the
three pairs whose truth is unambiguous (the unmixed controls).

Run:  python3 -I pairs.py <raw_dir> <out_dir>
"""
import csv
import itertools
import json
import os
import re
import sys

RAW, OUT = sys.argv[1], sys.argv[2]

PA = "PURE_ERR038255_100x"          # lineage 4.8.2
PB = "PURE_ERR13260175_100x"        # lineage 2.2.1
CT = "CTRL_SAMELIN_ERR038254_ERR038255_50pct_100x"   # same strain as PA

# Truth for the three pairs that contain no mixture of different lineages.
TRUTH = {
    frozenset((CT, PA)): "close",      # two libraries of the same strain
    frozenset((PA, PB)): "distant",    # different lineages, >1,000 SNPs
    frozenset((CT, PB)): "distant",
}


def key(name):
    s = str(name).strip().strip("'").strip('"')
    s = re.sub(r"^rMAPTB_benchmark\.", "", s)
    s = re.sub(r"^\d{4}_", "", s)
    s = re.sub(r"_lib\d+$", "", s)
    return s.replace("-", "_")


res = {}

# ---------------------------------------------------------------- rMAP-TB
p = os.path.join(RAW, "rmap", "rMAP-TB-V22", "cacheCopy-1", "snp_distance",
                 "pairwise_snp_distance_pairs.tsv")
rm = {}
with open(p, encoding="utf-8", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rm[frozenset((key(r["sample1"]), key(r["sample2"])))] = (
            int(r["snp_distance"]), r["cluster_class"])

res["rMAP-TB"] = {
    "pairs_scored": sum(1 for v in rm.values() if v[1] != "uninterpretable"),
    "pairs_refused": sum(1 for v in rm.values() if v[1] == "uninterpretable"),
    "called_close": sum(1 for v in rm.values() if v[1] == "close"),
    "controls": {" vs ".join(sorted(k)): {"snps": rm[k][0], "called": rm[k][1],
                                          "truth": t}
                 for k, t in TRUTH.items()},
}

# ----------------------------------------------------------------- MTBseq
mp = os.path.join(RAW, "mtbseq", "MTBSEQ", "mtbseq_distance_matrix.tsv")
names, tri = [], []
for line in open(mp, encoding="utf-8"):
    q = line.rstrip("\n").split("\t")
    names.append(key(q[0]))
    tri.append([int(x) for x in q[1:] if x.strip() != ""])


def md(a, b):
    ia, ib = names.index(a), names.index(b)
    i, j = max(ia, ib), min(ia, ib)
    return tri[i][j]


allm = [md(a, b) for a, b in itertools.combinations(names, 2)]
res["MTBseq"] = {
    "pairs_scored": len(allm),
    "pairs_refused": 0,
    "max_distance": max(allm),
    "called_close_at_5": sum(1 for d in allm if d <= 5),
    "called_close_at_12": sum(1 for d in allm if d <= 12),
    "controls": {" vs ".join(sorted(k)): {"snps": md(*tuple(k)),
                                          "called": "close" if md(*tuple(k)) <= 5
                                          else "distant", "truth": t}
                 for k, t in TRUTH.items()},
}

# ------------------------------------------------------------------ UShER
up = os.path.join(RAW, "usher", "USHER", "call-USHER_DISTANCES", "usher_dist",
                  "usher_clusters.tsv")
uc = {}
with open(up, encoding="utf-8", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        uc[key(r["sample"])] = r["cluster"]
samples = sorted(uc)
linked = sum(1 for a, b in itertools.combinations(samples, 2)
             if uc[a] == uc[b])
res["UShER"] = {
    "pairs_scored": len(samples) * (len(samples) - 1) // 2,
    "pairs_refused": 0,
    "pairs_in_same_cluster": linked,
    "controls": {" vs ".join(sorted(k)):
                 {"called": "close" if uc[list(k)[0]] == uc[list(k)[1]]
                  else "distant", "truth": t}
                 for k, t in TRUTH.items()},
}

# ------------------------------------------------------------------ MAGMA
cl = open(os.path.join(RAW, "magma", "collected", "magma_clusters.tsv"),
          encoding="utf-8", errors="replace").read()
mem = {}
for m in re.finditer(r"^(\d+)\t(\d+)\t\d+\t\[(.*?)\]\t", cl, re.M | re.S):
    for x in m.group(3).split(","):
        mem[key(x.strip())] = int(m.group(1))
res["MAGMA"] = {
    "note": "PURE_ERR038255_100x and the same-strain control produced no "
            "mapped reads, so none of the three control pairs can be scored",
    "clustered_samples": len(mem),
    "pure_parent_in_cluster": {k: v for k, v in mem.items()
                               if k.startswith("PURE_")},
    "controls_assessable": sum(1 for k in TRUTH
                               if all(x in mem for x in k)),
}


def errors(tool):
    c = res[tool].get("controls", {})
    return sum(1 for v in c.values() if v.get("called") != v["truth"])


print("Truth for the three unmixed control pairs: 1 close, 2 distant\n")
for t in ("rMAP-TB", "MTBseq", "UShER"):
    print("--- %s ---" % t)
    for name, v in res[t]["controls"].items():
        short = name.replace("_100x", "").replace(
            "CTRL_SAMELIN_ERR038254_ERR038255_50pct", "CTRL")
        mark = "OK " if v.get("called") == v["truth"] else "ERR"
        print("  %s %-44s truth=%-7s called=%-8s %s"
              % (mark, short, v["truth"], v.get("called"),
                 ("%d SNPs" % v["snps"]) if "snps" in v else ""))
    print("  errors: %d of 3\n" % errors(t))

print("--- MAGMA ---")
print("  %s" % res["MAGMA"]["note"])
print("  control pairs assessable: %d of 3\n"
      % res["MAGMA"]["controls_assessable"])

print("Pair-level totals over the 153 pairs")
print("  rMAP-TB : %3d scored, %3d refused, %3d called close"
      % (res["rMAP-TB"]["pairs_scored"], res["rMAP-TB"]["pairs_refused"],
         res["rMAP-TB"]["called_close"]))
print("  MTBseq  : %3d scored,   0 refused, %3d within 5 SNPs "
      "(largest pairwise distance in the cohort: %d SNPs)"
      % (res["MTBseq"]["pairs_scored"], res["MTBseq"]["called_close_at_5"],
         res["MTBseq"]["max_distance"]))
print("  UShER   : %3d scored,   0 refused, %3d in the same cluster"
      % (res["UShER"]["pairs_scored"], res["UShER"]["pairs_in_same_cluster"]))
print("  MAGMA   : %d of 18 libraries reached the clustering step"
      % res["MAGMA"]["clustered_samples"])

json.dump(res, open(os.path.join(OUT, "pairs.json"), "w"), indent=1)
print("\nwrote %s" % os.path.join(OUT, "pairs.json"))
