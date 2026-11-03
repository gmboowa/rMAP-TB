"""Benchmark figure: what four workflows did with the same 18 libraries.

A  mixed-infection detection per workflow
B  MTBseq's cohort control: sites retained and the parent-parent distance
C  the three unmixed control pairs, by workflow

Run:  python3 -I fig_bench.py <out_dir>
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = sys.argv[1]
C = json.load(open(os.path.join(OUT, "comparison.json")))
P = json.load(open(os.path.join(OUT, "pairs.json")))

INK, MUT = "#16202b", "#6b7785"
OK, BAD, NA = "#1f4e79", "#b91c1c", "#b6bec9"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": "#9aa7b4", "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUT, "ytick.color": MUT,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
})

fig = plt.figure(figsize=(7.1, 5.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92],
                      hspace=0.54, wspace=0.58,
                      left=0.145, right=0.985, top=0.925, bottom=0.095)

# ------------------------------------------------------------------ A
ax = fig.add_subplot(gs[0, :])
tools = ["rMAP-TB v22", "MAGMA", "MTBseq", "UShER"]
det = [len(C["rMAP-TB"]["detected_mixture"]),
       len(C["MAGMA"]["detected_mixture_strict"]),
       0, 0]
notes = ["",
         "+1 more had a single lineage\nbelow 90%; none of the 8 mixtures\nMAGMA itself passed were flagged",
         "one lineage per sample;\n'bad' classification quality on 6",
         "no mixture model"]
y = list(range(len(tools)))[::-1]
bars = ax.barh(y, det, height=0.52, color=[OK, BAD, BAD, BAD], zorder=3)
ax.axvline(15, color=MUT, lw=1.0, ls=(0, (4, 3)), zorder=2)
ax.text(15.15, max(y) + 0.52, "all 15 mixtures", fontsize=7.2, color=MUT,
        va="center")
ax.set_yticks(y)
ax.set_yticklabels(tools, fontsize=8.4)
ax.set_xlim(0, 26)
ax.set_xticks([0, 5, 10, 15])
ax.set_xlabel("Mixtures flagged as containing more than one strain (of 15)")
ax.set_title("A   Mixed-infection detection", fontsize=9, fontweight="bold",
             loc="left", pad=7)
for i, (v, n) in enumerate(zip(det, notes)):
    yy = y[i]
    ax.text(v + 0.35, yy, str(v), va="center", fontsize=8.6,
            fontweight="bold", color=OK if i == 0 else BAD)
    if n:
        ax.text(16.1, yy, n, va="center", fontsize=6.3, color=MUT,
                linespacing=1.25)
ax.set_axisbelow(True)
ax.xaxis.grid(True, color="#e4e9ef", lw=0.7)
ax.tick_params(axis="y", length=0)

# ------------------------------------------------------------------ B
ax = fig.add_subplot(gs[1, 0])
lab = ["parents alone\n(2 libraries)\n\u2192 ungrouped",
       "parents + 16 mixtures\n(18 libraries)\n\u2192 one group of 18"]
sites = [C["MTBseq"]["alignment_cols_2cohort"],
         C["MTBseq"]["alignment_cols_18cohort"]]
dist = [C["MTBseq"]["parent_parent_snps_2cohort"],
        C["MTBseq"]["parent_parent_snps_18cohort"]]
x = [0, 1]
ax.bar(x, sites, width=0.5, color=[OK, BAD], zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(lab, fontsize=7.1, linespacing=1.35)
ax.set_ylabel("Variant sites retained", fontsize=8.2)
ax.set_ylim(0, 1900)
ax.set_yticks([0, 500, 1000, 1500])
ax.set_title("B   MTBseq: the same two isolates,\n     two cohorts",
             fontsize=9, fontweight="bold", loc="left", pad=10)
for xi, sv, d in zip(x, sites, dist):
    inside = sv > 400
    ax.text(xi, sv / 2 if inside else sv + 80, format(sv, ","),
            ha="center", va="center", fontsize=7.6,
            fontweight="bold", color="white" if inside else MUT)
    ax.text(xi, sv + (150 if inside else 300),
            "%s SNPs apart" % format(d, ","), ha="center", fontsize=8.0,
            fontweight="bold", color=OK if d > 100 else BAD)
ax.set_axisbelow(True)
ax.yaxis.grid(True, color="#e4e9ef", lw=0.7)

# ------------------------------------------------------------------ C
ax = fig.add_subplot(gs[1, 1])
rows = ["CTRL vs own parent\ntruth: close",
        "parent vs parent\ntruth: distant",
        "CTRL vs other parent\ntruth: distant"]
cols = ["rMAP-TB", "MTBseq", "UShER", "MAGMA"]
# 1 correct, 0 wrong, -1 not assessable
grid = [[1, 1, 1, -1],
        [1, 0, 0, -1],
        [1, 0, 0, -1]]
for i, r in enumerate(grid):
    for j, v in enumerate(r):
        col = OK if v == 1 else (BAD if v == 0 else NA)
        ax.add_patch(plt.Rectangle((j, -i), 0.92, 0.92, color=col,
                                   zorder=3))
        txt = {1: "✓", 0: "✗", -1: "n/a"}[v]
        ax.text(j + 0.46, -i + 0.46, txt, ha="center", va="center",
                color="white", fontsize=11 if v != -1 else 7.5,
                fontweight="bold", zorder=4)
ax.set_xlim(-0.06, 4.0)
ax.set_ylim(-2.1, 1.0)
ax.set_xticks([j + 0.46 for j in range(4)])
ax.set_xticklabels(cols, fontsize=7.6)
ax.set_yticks([-i + 0.46 for i in range(3)])
ax.set_yticklabels(rows, fontsize=6.9, linespacing=1.3)
ax.set_title("C   Unmixed control pairs", fontsize=9,
             fontweight="bold", loc="left", pad=7)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
ax.legend(handles=[Patch(color=OK, label="correct"),
                   Patch(color=BAD, label="incorrect"),
                   Patch(color=NA, label="not assessable")],
          frameon=False, fontsize=6.9, ncol=3, loc="upper center",
          bbox_to_anchor=(0.5, -0.21), handlelength=1.1, columnspacing=1.0)

fig.savefig(os.path.join(OUT, "Figure6_benchmark.png"), dpi=600,
            facecolor="white")
fig.savefig(os.path.join(OUT, "Figure6_benchmark.pdf"), facecolor="white")
print("written")
