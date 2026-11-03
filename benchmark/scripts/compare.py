"""Head-to-head comparison of rMAP-TB, MAGMA, MTBseq and UShER on the shared
18-sample in-silico mixture benchmark.

Every value is read from a primary output file of the tool concerned. Nothing is
taken from the summary documents shipped alongside those outputs.

Run:  python3 -I compare.py <raw_dir> <out_dir>
"""
import csv
import json
import os
import re
import sys

RAW, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- the design
# Ground truth is the construction, not any tool's output.
# key -> (is_mixture, minor_pct, direction, nominal_depth, parent_pair)
DESIGN = {
    "CTRL_SAMELIN_ERR038254_ERR038255_50pct_100x":
        (False, 50, "same-lineage", 100, "ERR038254/ERR038255"),
    "MIX_ERR038255_ERR13260175_05pctA_100x": (True, 5, "A", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_05pctB_100x": (True, 5, "B", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_10pctA_100x": (True, 10, "A", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_10pctB_100x": (True, 10, "B", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_10pctB_140x": (True, 10, "B", 140, "pair1"),
    "MIX_ERR038255_ERR13260175_10pctB_40x": (True, 10, "B", 40, "pair1"),
    "MIX_ERR038255_ERR13260175_10pctB_70x": (True, 10, "B", 70, "pair1"),
    "MIX_ERR038255_ERR13260175_20pctA_100x": (True, 20, "A", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_20pctB_100x": (True, 20, "B", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_30pctA_100x": (True, 30, "A", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_30pctB_100x": (True, 30, "B", 100, "pair1"),
    "MIX_ERR038255_ERR13260175_50pctB_100x": (True, 50, "B", 100, "pair1"),
    "MIX_ERR038266_ERR13260106_10pctB_100x": (True, 10, "B", 100, "pair2"),
    "MIX_ERR038266_ERR13260106_20pctB_100x": (True, 20, "B", 100, "pair2"),
    "MIX_ERR038266_ERR13260106_50pctB_100x": (True, 50, "B", 100, "pair2"),
    "PURE_ERR038255_100x": (False, 0, "pure", 100, "ERR038255"),
    "PURE_ERR13260175_100x": (False, 0, "pure", 100, "ERR13260175"),
}
MIXTURES = [k for k, v in DESIGN.items() if v[0]]
CONTROLS = [k for k, v in DESIGN.items() if not v[0]]


def key(name):
    """Normalise a sample identifier from any of the four tools."""
    s = str(name).strip().strip("'").strip('"')
    s = re.sub(r"^rMAPTB_benchmark\.", "", s)
    s = re.sub(r"^\d{4}_", "", s)
    s = re.sub(r"_lib\d+$", "", s)
    s = s.replace("-", "_")
    return s


def tsv(path, delim="\t"):
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


R = {}            # results, tool -> {...}

# ================================================================== rMAP-TB
rm = os.path.join(RAW, "rmap", "rMAP-TB-V22")
elig = tsv(os.path.join(rm, "cacheCopy-1", "snp_distance",
                        "clustering_eligibility.tsv"))
rmap_mixed, rmap_elig = {}, {}
for r in elig:
    k = key(r["sample"])
    rmap_mixed[k] = r["tbprofiler_mixed_infection_status"]
    rmap_elig[k] = r["eligible_for_clustering"]

pairs = tsv(os.path.join(rm, "cacheCopy-1", "snp_distance",
                         "pairwise_snp_distance_pairs.tsv"))
rmap_pairs = {}
for r in pairs:
    a, b = key(r["sample1"]), key(r["sample2"])
    rmap_pairs[frozenset((a, b))] = {
        "snps": int(r["snp_distance"]),
        "sites": int(r["comparable_sites"]),
        "class": r["cluster_class"],
    }

R["rMAP-TB"] = {
    "version": "v22",
    "samples_with_output": len(rmap_mixed),
    "detected_mixture": [k for k in MIXTURES
                         if rmap_mixed.get(k, "").lower().startswith("possible mixed")],
    "controls_flagged": [k for k in CONTROLS
                         if rmap_mixed.get(k, "").lower().startswith("possible mixed")],
    "pairs_total": len(rmap_pairs),
    "pairs_refused": sum(1 for v in rmap_pairs.values()
                         if v["class"] == "uninterpretable"),
    "pairs_close": sum(1 for v in rmap_pairs.values() if v["class"] == "close"),
    "pairs_distant": sum(1 for v in rmap_pairs.values() if v["class"] == "distant"),
    "parent_parent_snps": rmap_pairs[frozenset(
        ("PURE_ERR038255_100x", "PURE_ERR13260175_100x"))]["snps"],
    "parent_parent_sites": rmap_pairs[frozenset(
        ("PURE_ERR038255_100x", "PURE_ERR13260175_100x"))]["sites"],
    "parent_parent_class": rmap_pairs[frozenset(
        ("PURE_ERR038255_100x", "PURE_ERR13260175_100x"))]["class"],
}

# ==================================================================== MAGMA
mg = os.path.join(RAW, "magma", "collected")
stats = tsv(os.path.join(mg, "magma_cohort_stats.tsv"))
magma = {}
for r in stats:
    k = key(r["SAMPLE"])
    lins = [x for x in (r.get("LINEAGES") or "").split(";") if x]
    freqs = [x for x in (r.get("FREQUENCIES") or "").split(";") if x]
    fnum = []
    for f in freqs:
        try:
            fnum.append(float(f.rstrip("%")))
        except ValueError:
            pass
    magma[k] = {
        "mapped_pct": float(r["MAPPED_PERCENTAGE"]),
        "mean_cov": float(r["MEAN_COVERAGE"]),
        "lineages": lins,
        "freqs": fnum,
        "all_thresholds_met": r["ALL_THRESHOLDS_MET"] == "1",
    }

magma_pass = [k for k, v in magma.items() if v["all_thresholds_met"]]
magma_pass_mix = [k for k in magma_pass if k in MIXTURES]


def magma_multi(k, lenient=False):
    v = magma.get(k, {})
    if len(v.get("lineages", [])) > 1:
        return True
    if lenient and v.get("freqs") and max(v["freqs"]) < 90.0:
        return True
    return False


R["MAGMA"] = {
    "samples_with_output": sum(1 for v in magma.values() if v["mapped_pct"] > 0),
    "samples_passing_own_thresholds": len(magma_pass),
    "samples_zero_mapped": [k for k, v in magma.items() if v["mapped_pct"] == 0.0],
    "detected_mixture_strict": [k for k in MIXTURES if magma_multi(k)],
    "detected_mixture_lenient": [k for k in MIXTURES if magma_multi(k, True)],
    "passing_mixtures": magma_pass_mix,
    "detected_among_passing_strict": [k for k in magma_pass_mix if magma_multi(k)],
    "detected_among_passing_lenient": [k for k in magma_pass_mix
                                       if magma_multi(k, True)],
    "controls_flagged": [k for k in CONTROLS if magma_multi(k, True)],
    "mapped_pct": {k: v["mapped_pct"] for k, v in magma.items()},
}

# read MAGMA's ClusterPicker output
cl = open(os.path.join(mg, "magma_clusters.tsv"), encoding="utf-8",
          errors="replace").read()
clusters = []
for m in re.finditer(r"^(\d+)\t(\d+)\t\d+\t\[(.*?)\]\t([\d.]+)\t",
                     cl, re.M | re.S):
    members = [key(x.strip()) for x in m.group(3).split(",")]
    clusters.append({"id": int(m.group(1)), "n": int(m.group(2)),
                     "members": members, "bootstrap": float(m.group(4))})
R["MAGMA"]["clusters"] = clusters
R["MAGMA"]["cluster_with_pure_parent"] = [
    c for c in clusters
    if any(x.startswith("PURE_") for x in c["members"]) and c["n"] > 1]

# =================================================================== MTBseq
ms = os.path.join(RAW, "mtbseq", "MTBSEQ")
cls = tsv(os.path.join(ms, "mtbseq_classification.tsv"))
mtbseq_cls = {}
for r in cls:
    k = key(r["SampleID"])
    mtbseq_cls[k] = {
        "lineage": (r.get("Coll lineage (branch)") or "").strip("'"),
        "name": (r.get("Coll lineage_name (branch)") or "").strip("'"),
        "quality": (r.get("Coll quality (branch)") or "").strip("'"),
    }

st = tsv(os.path.join(ms, "mtbseq_statistics.tsv"))
mtbseq_stats = {}
for r in st:
    k = key(r["SampleID"])
    mtbseq_stats[k] = {
        "mapped_pct": float((r["% Mapped Reads"] or "0").strip("'")),
        "cov": float((r["(Unambiguous) Coverage mean"] or "0").strip("'")),
    }

# the 18-sample group assignment
grp = open(os.path.join(ms, "mtbseq_groups.tsv"), encoding="utf-8").read()
glist = {}
for line in grp.split("### Output as lists:")[-1].strip().split("\n"):
    p = line.split("\t")
    if len(p) == 2:
        glist[key(p[0])] = p[1]

# the 18-sample distance matrix (lower triangle, MTBseq format)
mat_path = os.path.join(ms, "mtbseq_distance_matrix.tsv")
names, tri = [], []
for line in open(mat_path, encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    names.append(key(p[0]))
    tri.append([int(x) for x in p[1:] if x.strip() != ""])


def mt_dist(a, b):
    ia, ib = names.index(a), names.index(b)
    i, j = max(ia, ib), min(ia, ib)
    return tri[i][j]


def aln_len(path):
    n, L = 0, None
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            n += 1
        elif L is None:
            L = len(line.strip())
    return n, L


amend18 = os.path.join(
    ms, "MTBseq_Amend_Audit",
    "rMAPTB_benchmark_joint_cf4_cr4_fr75_ph4_samples18_amended_u95_phylo"
    ".plainIDs.fasta")
amend2 = os.path.join(
    ms, "MTBseq_pure_only", "Amend",
    "rMAPTB_pure2_joint_cf4_cr4_fr75_ph4_samples2_amended_u95_phylo"
    ".plainIDs.fasta")
n18, L18 = aln_len(amend18)
n2, L2 = aln_len(amend2)


def tab_sites(path):
    """MTBseq .tab site tables: data rows below the two header lines."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return max(sum(1 for _ in fh) - 2, 0)


joint18 = tab_sites(os.path.join(
    ms, "MTBseq_Extra_Audit",
    "rMAPTB_benchmark_joint_cf4_cr4_fr75_ph4_samples18.tab"))
joint2 = tab_sites(os.path.join(
    ms, "MTBseq_pure_only", "Amend",
    "rMAPTB_pure2_joint_cf4_cr4_fr75_ph4_samples2_amended.tab"))

# pure-only cohort matrix
pm = os.path.join(ms, "MTBseq_pure_only", "Groups")
pure_matrix = [f for f in os.listdir(pm) if f.endswith(".matrix")][0]
pure_tri = []
for line in open(os.path.join(pm, pure_matrix), encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    pure_tri.append([x for x in p[1:] if x.strip() != ""])
pure_parent_dist = int(pure_tri[1][0])
pure_groups = open(
    [os.path.join(pm, f) for f in os.listdir(pm) if f.endswith(".groups")][0],
    encoding="utf-8").read()

R["MTBseq"] = {
    "samples_with_output": len(mtbseq_stats),
    "mapped_pct_min": min(v["mapped_pct"] for v in mtbseq_stats.values()),
    "mapped_pct_max": max(v["mapped_pct"] for v in mtbseq_stats.values()),
    "cov_min": min(v["cov"] for v in mtbseq_stats.values()),
    "cov_max": max(v["cov"] for v in mtbseq_stats.values()),
    "detected_mixture": [],           # MTBseq reports one lineage per sample
    "quality_bad": [k for k, v in mtbseq_cls.items() if v["quality"] == "bad"],
    "quality_bad_mixtures": [k for k in MIXTURES
                             if mtbseq_cls.get(k, {}).get("quality") == "bad"],
    "groups_18": {g: sum(1 for v in glist.values() if v == g)
                  for g in set(glist.values())},
    "parent_parent_snps_18cohort": mt_dist("PURE_ERR038255_100x",
                                           "PURE_ERR13260175_100x"),
    "parent_parent_same_group_18cohort":
        glist["PURE_ERR038255_100x"] == glist["PURE_ERR13260175_100x"],
    "alignment_cols_18cohort": L18,
    "alignment_seqs_18cohort": n18,
    "joint_candidate_sites_18cohort": joint18,
    "parent_parent_snps_2cohort": pure_parent_dist,
    "alignment_cols_2cohort": L2,
    "joint_candidate_sites_2cohort": joint2,
    "pure_cohort_grouping": "ungrouped" in pure_groups,
}

# ==================================================================== UShER
us = os.path.join(RAW, "usher", "USHER")
ucl = tsv(os.path.join(us, "call-USHER_DISTANCES", "usher_dist",
                       "usher_clusters.tsv"))
uclusters = {}
for r in ucl:
    uclusters.setdefault(r["cluster"], []).append(key(r["sample"]))
usamp = tsv(os.path.join(us, "call-USHER_PLACE", "usher_out", "samples.tsv"))
audit = open(os.path.join(us, "call-USHER_DISTANCES", "usher_dist",
                          "audit_summary.txt"), encoding="utf-8").read()

R["UShER"] = {
    "samples_with_output": len([r for r in usamp
                                if key(r["sample"]) in DESIGN]),
    "detected_mixture": [],           # placement only; no mixture model
    "clusters": {k: len(v) for k, v in uclusters.items()},
    "largest_cluster": max(len(v) for v in uclusters.values()),
    "threshold": ucl[0]["threshold_mutation_events"],
    "parents_same_cluster":
        uclusters[[c for c, v in uclusters.items()
                   if "PURE_ERR038255_100x" in v][0]]
        == uclusters[[c for c, v in uclusters.items()
                      if "PURE_ERR13260175_100x" in v][0]],
    "parsimony_zero": [key(r["sample"]) for r in usamp
                       if r["parsimony"] == "0" and key(r["sample"]) in DESIGN],
    "vcf_snp_records": int(re.search(r"VCF SNP records: (\d+)", audit).group(1)),
    "fasta_columns": int(re.search(r"FASTA columns: (\d+)", audit).group(1)),
    "caveat": [ln.strip() for ln in audit.strip().split("\n")
               if ln.strip().startswith(("UShER parsimony", "Shared input"))],
}

R["_design"] = {"n_samples": len(DESIGN), "n_mixtures": len(MIXTURES),
                "n_controls": len(CONTROLS), "mixtures": MIXTURES,
                "controls": CONTROLS}

json.dump(R, open(os.path.join(OUT, "comparison.json"), "w"), indent=1)

# ------------------------------------------------------------------- report
def pr(*a):
    print(*a)


pr("DESIGN: %d samples = %d mixtures + %d controls\n"
   % (len(DESIGN), len(MIXTURES), len(CONTROLS)))

pr("--- mixed-infection detection (of %d mixtures) ---" % len(MIXTURES))
pr("rMAP-TB v22 : %2d/%d" % (len(R["rMAP-TB"]["detected_mixture"]), len(MIXTURES)))
pr("MAGMA strict: %2d/%d   lenient (<90%% freq counts): %d/%d"
   % (len(R["MAGMA"]["detected_mixture_strict"]), len(MIXTURES),
      len(R["MAGMA"]["detected_mixture_lenient"]), len(MIXTURES)))
pr("  among the %d mixtures MAGMA itself passed: strict %d, lenient %d"
   % (len(magma_pass_mix), len(R["MAGMA"]["detected_among_passing_strict"]),
      len(R["MAGMA"]["detected_among_passing_lenient"])))
pr("MTBseq      :  0/%d  (reports one lineage per sample; "
   "'bad' quality on %d mixtures)"
   % (len(MIXTURES), len(R["MTBseq"]["quality_bad_mixtures"])))
pr("UShER       :  0/%d  (no mixture model)" % len(MIXTURES))
pr("false positives on the %d controls: rMAP-TB %d, MAGMA %d"
   % (len(CONTROLS), len(R["rMAP-TB"]["controls_flagged"]),
      len(R["MAGMA"]["controls_flagged"])))

pr("\n--- the ground-truth control: two unrelated parents ---")
pr("rMAP-TB v22 : %d SNPs over %s comparable sites -> %s"
   % (R["rMAP-TB"]["parent_parent_snps"],
      format(R["rMAP-TB"]["parent_parent_sites"], ","),
      R["rMAP-TB"]["parent_parent_class"]))
pr("MTBseq 18   : %d SNPs over %d alignment columns -> same group: %s"
   % (R["MTBseq"]["parent_parent_snps_18cohort"],
      R["MTBseq"]["alignment_cols_18cohort"],
      R["MTBseq"]["parent_parent_same_group_18cohort"]))
pr("MTBseq  2   : %d SNPs over %d alignment columns -> ungrouped: %s"
   % (R["MTBseq"]["parent_parent_snps_2cohort"],
      R["MTBseq"]["alignment_cols_2cohort"],
      R["MTBseq"]["pure_cohort_grouping"]))
pr("UShER       : same cluster: %s (threshold %s mutation events)"
   % (R["UShER"]["parents_same_cluster"], R["UShER"]["threshold"]))
pr("MAGMA       : not measurable - PURE_ERR038255_100x had %.2f%% reads mapped"
   % magma["PURE_ERR038255_100x"]["mapped_pct"])

pr("\n--- clustering output ---")
pr("rMAP-TB v22 : %d pairs = %d refused, %d close, %d distant"
   % (R["rMAP-TB"]["pairs_total"], R["rMAP-TB"]["pairs_refused"],
      R["rMAP-TB"]["pairs_close"], R["rMAP-TB"]["pairs_distant"]))
pr("MTBseq      : groups %s" % R["MTBseq"]["groups_18"])
pr("UShER       : clusters %s" % R["UShER"]["clusters"])
pr("MAGMA       : %d cluster(s); containing a pure parent: %s"
   % (len(clusters),
      [(c["id"], c["n"]) for c in R["MAGMA"]["cluster_with_pure_parent"]]))

pr("\n--- cohort completeness ---")
for t in ("rMAP-TB", "MAGMA", "MTBseq", "UShER"):
    pr("%-12s: %d/18 samples with usable output" % (t, R[t]["samples_with_output"]))
pr("MAGMA zero-mapped: %s" % R["MAGMA"]["samples_zero_mapped"])

pr("\n--- sites available to the clustering decision ---")
pr("rMAP-TB v22 : %s comparable alignment columns (parent pair)"
   % format(R["rMAP-TB"]["parent_parent_sites"], ","))
pr("MTBseq 18   : %d of %d joint candidate sites survived callability (%.1f%% lost)"
   % (R["MTBseq"]["alignment_cols_18cohort"],
      R["MTBseq"]["joint_candidate_sites_18cohort"],
      100 * (1 - R["MTBseq"]["alignment_cols_18cohort"]
             / R["MTBseq"]["joint_candidate_sites_18cohort"])))
pr("MTBseq  2   : %d of %d (%.1f%% lost)"
   % (R["MTBseq"]["alignment_cols_2cohort"],
      R["MTBseq"]["joint_candidate_sites_2cohort"],
      100 * (1 - R["MTBseq"]["alignment_cols_2cohort"]
             / R["MTBseq"]["joint_candidate_sites_2cohort"])))
pr("UShER       : %d VCF SNP records from %s alignment columns"
   % (R["UShER"]["vcf_snp_records"], format(R["UShER"]["fasta_columns"], ",")))

pr("\nwrote %s" % os.path.join(OUT, "comparison.json"))
