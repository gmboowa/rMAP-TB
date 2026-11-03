"""Figure S2: pairwise core-SNP distance heatmap for Panel 2.

Built from the workflow's own pairwise table (the 78 pairs among the 13
libraries that entered the core alignment). Distances span four orders of
magnitude because M. canettii sits far outside the rest, so the colour scale
is logarithmic and every cell carries its value.

Run:  python3 fig_s2.py <panel2_json> <out_dir>
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

SRC = Path(sys.argv[1])
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "out")
OUT.mkdir(parents=True, exist_ok=True)

rows = json.load(open(SRC))["snpPairsTable"][1:]

LINEAGE = {
    "CONTAM_ERR038255_ERR12861768_10pctNTM_100x": ("10 % NTM contam.", "L4.8.2"),
    "CONTAM_ERR038255_ERR12861768_50pctNTM_100x": ("50 % NTM contam.", "L4.8.2"),
    "ERR1203067": ("ERR1203067", "L4"),
    "ERR14187598": ("ERR14187598", "La1"),
    "ERR14187918": ("ERR14187918", "La1"),
    "ERR266109": ("ERR266109", "M. canettii"),
    "ERR3804031": ("ERR3804031", "L6"),
    "ERR970410": ("ERR970410", "La1"),
    "SRR25473322": ("SRR25473322", "La1"),
    "SRR29440878": ("SRR29440878", "La1 (BCG)"),
    "SRR3500411": ("SRR3500411", "La1"),
    "SRR35607801": ("SRR35607801", "La2"),
    "SRR9157804": ("SRR9157804", "La3"),
}

names = sorted({r[0] for r in rows} | {r[1] for r in rows},
               key=lambda s: (LINEAGE[s][1] == "M. canettii", LINEAGE[s][1], s))
idx = {n: i for i, n in enumerate(names)}
N = len(names)

M = np.full((N, N), np.nan)
for a, b, d, _sites, _interp, _cls in rows:
    i, j = idx[a], idx[b]
    M[i, j] = M[j, i] = int(d)
np.fill_diagonal(M, np.nan)   # self-comparisons are not pairs

n_pairs = int(np.isfinite(M).sum()) // 2
assert n_pairs == len(rows) == 78, n_pairs

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.linewidth": 0.6, "savefig.dpi": 600,
})
INK, MUT = "#16202b", "#6b7785"

fig, ax = plt.subplots(figsize=(8.1, 6.6))

# Discrete bins, because the epidemiologically meaningful contrast is between
# "inside a transmission threshold" and "not", and then between degrees of
# relatedness within the complex. A continuous scale is dominated by the
# M. canettii outgroup and hides all of that.
BOUNDS = [0, 6, 500, 1000, 1500, 2000, 2500, 20000]
COLORS = ["#b91c1c",   # <= 5 SNPs: inside a transmission threshold
          "#e8f1f8", "#bcd8ea", "#7fb3d5", "#4a8bbd", "#2e6ca3", "#17365d"]
cmap = ListedColormap(COLORS)
cmap.set_bad("#f4f6f9")
norm = BoundaryNorm(BOUNDS, cmap.N)

im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, norm=norm)

labels = ["%s  \u00b7  %s" % (LINEAGE[n][0], LINEAGE[n][1]) for n in names]
ax.set_xticks(range(N))
ax.set_yticks(range(N))
ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=7.2, color=INK)
ax.set_yticklabels(labels, fontsize=7.2, color=INK)
ax.set_xticks(np.arange(-0.5, N, 1), minor=True)
ax.set_yticks(np.arange(-0.5, N, 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.1)
ax.tick_params(which="minor", length=0)
ax.tick_params(which="major", length=2, color=MUT)
for sp in ax.spines.values():
    sp.set_color("#c9d2dc")

for i in range(N):
    for j in range(N):
        if i == j:
            ax.text(j, i, "\u00b7", ha="center", va="center",
                    color="#aab4c0", fontsize=9)
            continue
        v = M[i, j]
        if not np.isfinite(v):
            continue
        txt = "%d" % v if v < 10000 else "%.1fk" % (v / 1000)
        dark = v <= 5 or v >= 1500
        ax.text(j, i, txt, ha="center", va="center", fontsize=6.0,
                color="white" if dark else INK)

cb = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.015,
                  boundaries=BOUNDS, spacing="uniform")
cb.set_ticks(BOUNDS)
cb.ax.set_yticklabels(["0", "5", "500", "1,000", "1,500", "2,000", "2,500",
                       "20,000"], fontsize=6.8, color=INK)
cb.set_label("pairwise core-SNP distance", fontsize=7.4, color=INK, labelpad=6)
cb.outline.set_linewidth(0.5)
cb.ax.tick_params(length=2, color=MUT)

ax.set_title("Panel 2: pairwise core-SNP distances among the 13 libraries\n"
             "that entered the core alignment",
             fontsize=9.0, color=INK, pad=10)

close = [(a, b, int(d)) for a, b, d, _s, _i, c in rows if c == "close"]
fig.subplots_adjust(left=0.235, right=0.90, top=0.915, bottom=0.235)

png = OUT / "Figure_S2_panel2_heatmap.png"
fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.06)
fig.savefig(OUT / "Figure_S2_panel2_heatmap.pdf", bbox_inches="tight", pad_inches=0.06)
plt.close(fig)

print("%d libraries, %d pairs" % (N, n_pairs))
print("close pairs (<=5 SNPs): %s" % close)
fin = M[np.isfinite(M)]
print("range %d - %d SNPs" % (fin.min(), fin.max()))
print("wrote %s" % png.name)
