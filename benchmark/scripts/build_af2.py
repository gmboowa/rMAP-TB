"""Assemble Additional file 2: supplementary figures S1-S2, one figure per page
with its legend, A4 portrait.

Run:  python3 build_af2.py <out_dir>
"""
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "out")

A4 = (8.27, 11.69)
INK, MUT = "#16202b", "#4a5462"

PAGES = [
    (OUT / "Figure_S1_report_states.png",
     "Figure S1. Report states for rendering failure, incomplete annotation "
     "and the legacy placeholder image.",
     "Each panel is the phylogeny section of an integrated report built from "
     "synthetic inputs by the repository's test-suite, which asserts the "
     "state shown. (A) The tree renderer failed: the report gives the error, "
     "the renderer status keys and the renderer log in place of an empty "
     "image, and states separately that IQ-TREE itself succeeded, so a "
     "drawing failure is not read as an analysis failure. (B) The tree "
     "rendered but no tip matched the TB-Profiler metadata; the topology is "
     "shown, the annotation failure is named with its underlying error, and "
     "the reader is directed to the per-isolate section for lineage and "
     "resistance. (C) A 68-byte single-pixel placeholder image, produced by "
     "an earlier version of the renderer, is recognised as unusable and is "
     "not embedded as though it were a tree. Produced by "
     "benchmark/scripts/gen_states.py, which calls the same fixtures and "
     "report builder as tests/test_merge_report.py and re-checks each test's "
     "assertions."),

    (OUT / "Figure_S2_panel2_heatmap.png",
     "Figure S2. Pairwise core-SNP distance heatmap for Panel 2.",
     "All 78 pairwise distances among the 13 libraries that entered the core "
     "alignment, from the workflow's own pairwise distance table. Labels give "
     "the accession and the lineage call. Colour is binned rather than "
     "continuous, because a continuous scale is dominated by the "
     "M. canettii outgroup at roughly 16,000 SNPs and hides the structure "
     "within the complex. Red marks the only pair at or below the 5-SNP "
     "transmission threshold: the two contamination mixtures built from the "
     "same host library, 0 SNPs apart. SRR1239339 was removed by the "
     "post-mapping coverage gate before alignment and so does not appear. "
     "Comparable sites ranged from 3,703,794 to 3,962,824 across the pairs."),
]

pdf_path = OUT / "Additional_file_2_Figures_S1-S2.pdf"
with PdfPages(pdf_path) as pdf:
    for img_path, title, legend in PAGES:
        img = mpimg.imread(img_path)
        ar = img.shape[0] / img.shape[1]

        fig = plt.figure(figsize=A4)
        wrapped = textwrap.wrap(legend, 108)
        legend_h = (len(wrapped) + 1) * 0.155 / A4[1]

        title_h = 0.030 + 0.017 * (len(textwrap.wrap(title, 78)) - 1)
        avail_h = 1 - title_h - 0.030 - 0.030 - legend_h - 0.028
        avail_w = 1 - 2 * 0.055
        w = min(avail_w, avail_h / ar * (A4[1] / A4[0]))
        h = w * ar * (A4[0] / A4[1])
        left = (1 - w) / 2
        box_bottom = 0.030 + legend_h + 0.028
        bottom = box_bottom + (avail_h - h) / 2    # centre in the free space

        ax = fig.add_axes([left, bottom, w, h])
        ax.imshow(img)
        ax.axis("off")

        tlines = textwrap.wrap(title, 78)
        ty = 0.972
        for tline in tlines:
            fig.text(0.06, ty, tline, fontsize=9.6, color=INK,
                     fontweight="bold", va="top", family="DejaVu Sans")
            ty -= 0.017
        y = 0.030 + legend_h
        for line in wrapped:
            y -= 0.0133
            fig.text(0.06, y, line, fontsize=7.6, color=MUT,
                     va="top", family="DejaVu Sans")
        pdf.savefig(fig, dpi=600)
        plt.close(fig)
        print("page: %s  (%d x %d px)" % (img_path.name, img.shape[1], img.shape[0]))

    d = pdf.infodict()
    d["Title"] = "Additional file 2. Supplementary figures S1-S2"
    d["Subject"] = ("rMAP-TB: a routed, eligibility-gated WDL/Cromwell "
                    "workflow for MTBC and NTM genomic surveillance")

print("\nwrote %s  (%.2f MB)" % (pdf_path.name, pdf_path.stat().st_size / 1e6))
