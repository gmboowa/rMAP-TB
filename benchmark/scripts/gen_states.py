"""Produce the three report states shown in Additional file 2, Figure S1.

The states are exactly those asserted by the repository's own test-suite:

  test_renderer_failure_is_shown_with_log        -> rendering failure
  test_unannotated_tree_is_flagged               -> incomplete annotation
  test_legacy_one_pixel_png_is_not_shown_as_a_tree -> legacy placeholder image

This driver calls the same fixture and the same build_report() helper those
tests use, so the HTML written here is what the test-suite checks. It then
re-asserts each test's own conditions, so a state that silently stopped
behaving as documented would fail here rather than be photographed.

Run:  python3 gen_states.py <out_dir>
"""
import os
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "_pytest_stub"))      # stub pytest
REPO = HERE.parents[1]                              # benchmark/scripts -> repo root
sys.path.insert(0, str(REPO / "tests"))
sys.path.insert(0, str(REPO / "scripts"))

import conftest                                      # noqa: E402
import test_merge_report as T                        # noqa: E402
from summarize_report import summarize               # noqa: E402

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else HERE / "out")
OUT.mkdir(parents=True, exist_ok=True)


class Factory:
    """Stands in for pytest's tmp_path_factory."""

    def __init__(self, root):
        self.root = Path(root)

    def mktemp(self, name):
        p = Path(tempfile.mkdtemp(prefix=name + "_", dir=self.root))
        return p


work = Path(tempfile.mkdtemp(prefix="rmaptb_states_"))
fac = Factory(work)

print("building the synthetic panel ...")
panel = conftest.panel(fac)
print("  %d MTBC + %d NTM libraries"
      % (len(conftest.SAMPLES), len(conftest.NTM_SAMPLES)))

print("rendering the reference tree ...")
tree = T.rendered_tree(fac, panel)
print("  %s" % (tree / "phylogenetic_tree.svg").name)

states = {}

# ------------------------------------------------- A. renderer failure
d = work / "a"
d.mkdir()
status = ("svg_status=failed\npng_status=failed\npng_renderer=matplotlib-placeholder\n"
          "tips_rendered=0\nmidpoint_rooted=false\nreference_tips_removed=0\n"
          "metadata_status=ok\nmetadata_rows=13\ntips_matched=0\n"
          "error=ValueError('IQ-TREE did not produce a valid phylogenetic tree.')\n")
log = d / "render.log"
log.write_text("TREE_VISUALIZATION started\nERROR: ValueError(...)\n"
               "Traceback (most recent call last):\n")
ph = d / "placeholder.svg"
ph.write_text('<svg xmlns="http://www.w3.org/2000/svg">'
              '<text>Phylogenetic tree was not rendered</text></svg>')
path, html = T.build_report(d, panel, svg=ph, status=status, log=log)
s = summarize(path)
assert s["tree_rendering"]["state"] == "Rendering failed", s["tree_rendering"]
assert "Phylogenetic tree image was not rendered" in html
assert "Tree renderer log (render.log)" in html
assert "IQ-TREE did not produce a valid phylogenetic tree" in html
assert 'id="treeFigure"' not in html
states["A_rendering_failure"] = path
print("A  rendering failure            -> %s" % s["tree_rendering"]["state"])

# --------------------------------------------- B. annotations incomplete
d = work / "b"
d.mkdir()
status = (tree / "render_status.txt").read_text()
status = status.replace(
    "metadata_status=ok",
    "metadata_status=failed: Error('field larger than field limit (131072)')")
status = status.replace("tips_matched=13", "tips_matched=0")
path, html = T.build_report(d, panel, svg=tree / "phylogenetic_tree.svg",
                            status=status, log=tree / "render.log")
s = summarize(path)
assert s["tree_rendering"]["state"] == "Rendered, annotations incomplete", s["tree_rendering"]
assert "Tree annotations incomplete" in html
assert "field larger than field limit" in html
states["B_annotations_incomplete"] = path
print("B  incomplete annotation        -> %s" % s["tree_rendering"]["state"])

# ------------------------------------------- C. legacy placeholder image
d = work / "c"
d.mkdir()
path, html = T.build_report(d, panel, png=T.ONE_PIXEL_PNG)
s = summarize(path)
assert s["tree_rendering"]["state"] == "Rendering failed", s["tree_rendering"]
assert "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAAB" not in html, \
    "the 68-byte placeholder was embedded as if it were a tree"
states["C_legacy_placeholder"] = path
print("C  legacy placeholder image     -> %s (not embedded)"
      % s["tree_rendering"]["state"])

for name, p in states.items():
    dest = OUT / ("state_%s.html" % name)
    shutil.copy(p, dest)
    print("wrote %s  (%.1f KB)" % (dest.name, dest.stat().st_size / 1024))

print("\nall three states reproduced and re-asserted")
