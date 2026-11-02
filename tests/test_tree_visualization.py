"""Runs the TREE_VISUALIZATION Python program on synthetic trees.

ETE3 is not required: when it is absent the task must still produce a complete
SVG tree, a legible placeholder raster image and an honest render_status.txt.
"""

import re
import struct
import zlib

import pytest

from wdl_blocks import python_heredoc, run_python, substitute


def render(tmp_path, panel, newick=None, tb=None, res=None, midpoint="true", fmt="png", label_max=100):
    cwd = tmp_path
    tree = cwd / "input.nwk"
    tree.write_text(newick if newick is not None else panel["tree"].read_text())
    src = substitute(python_heredoc("TREE_VISUALIZATION"), {
        "image_format": fmt,
        "input_tree": str(tree),
        "tbprofiler_summary_tsv": str(tb if tb is not None else panel["tb"]),
        "resistance_profile_summary_tsv": str(res if res is not None else panel["res"]),
        "iqtree_excluded_samples_tsv": "",
        "width": "2400", "height": "1600", "max_bootstrap_to_label": str(label_max),
        "title": "MTBC Core-SNP Phylogenetic Tree", "midpoint_root_tree": midpoint,
    })
    r = run_python(src, cwd)
    assert r.returncode == 0, r.stderr
    out = cwd / "tree_visualization"
    status = dict(line.split("=", 1) for line in (out / "render_status.txt").read_text().splitlines() if "=" in line)
    return out, status


def png_size(path):
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    w, h = struct.unpack(">II", data[16:24])
    return w, h


def test_svg_tree_rendered_with_annotations(tmp_path, panel):
    out, status = render(tmp_path, panel)
    svg = (out / "phylogenetic_tree.svg").read_text()
    assert status["svg_status"] == "success"
    assert status["tips_rendered"] == "13"
    assert status["tips_matched"] == "13"
    assert status["reference_tips_removed"] == "1"
    assert status["midpoint_rooted"] == "true"
    assert svg.count("<circle") >= 13
    for name in ("ERR266109", "SRR9157804", "CONTAM_ERR038255_ERR12861768_10pctNTM_100x"):
        assert name in svg
    assert "Reference" not in re.sub(r"reference tip removed: Reference", "", svg)
    assert "Lineage NA" not in svg
    for cat in ("MDR-TB", "RR-TB", "Pre-XDR-TB", "XDR-TB", "Hr-TB", "Sensitive"):
        assert f">{cat}</text>" in svg
    assert "substitutions / site" in svg
    # Support labels: the SH-aLRT/UFBoot pair must never be concatenated.
    assert "100100" not in svg
    assert ">58<" in svg and ">51<" in svg


def test_large_tbprofiler_field_does_not_strip_annotations(tmp_path, panel):
    # Regression for the 20-sample run: a >131072-byte field raised csv.Error
    # inside the renderer and every tip became "Lineage NA / Unknown".
    out, status = render(tmp_path, panel, tb=panel["tb_big"])
    assert status["metadata_status"] == "ok"
    assert status["tips_matched"] == "13"
    assert "Lineage NA" not in (out / "phylogenetic_tree.svg").read_text()


def test_placeholder_raster_is_legible_not_one_pixel(tmp_path, panel):
    out, status = render(tmp_path, panel)
    png = out / "phylogenetic_tree.png"
    if status["png_status"] == "success":
        pytest.skip("ETE3 available locally; placeholder path not exercised")
    assert status["png_renderer"].endswith("placeholder")
    assert png.stat().st_size > 5000
    w, h = png_size(png)
    assert w >= 1000 and h >= 300


def test_cleaned_newick_is_reparsable_and_reference_free(tmp_path, panel):
    out, _ = render(tmp_path, panel)
    nwk = (out / "phylogenetic_tree.cleaned.nwk").read_text()
    assert nwk.strip().endswith(";")
    assert "Reference" not in nwk
    assert nwk.count("(") == nwk.count(")")
    assert "0000_" not in nwk  # Terra shard prefixes removed


def test_unrooted_when_midpoint_disabled(tmp_path, panel):
    _, status = render(tmp_path, panel, midpoint="false")
    assert status["midpoint_rooted"] == "false"
    assert status["svg_status"] == "success"


def test_svg_output_format_skips_qt(tmp_path, panel):
    out, status = render(tmp_path, panel, fmt="svg")
    assert status["png_status"] == "success"
    assert status["png_renderer"] == "builtin-svg"
    assert (out / "phylogenetic_tree.svg").read_text().startswith("<svg")


def test_bootstrap_threshold_hides_fully_supported_nodes(tmp_path, panel):
    out, _ = render(tmp_path, panel, label_max=99)
    svg = (out / "phylogenetic_tree.svg").read_text()
    assert ">100</text>" not in svg
    assert ">58</text>" in svg


@pytest.mark.parametrize("newick,tips", [
    ("('A sample':0.1,'B''s':0.2,(C:0.05,D:0.07)0.95:0.1);", 4),
    ("(A:0.1,B:0.2);", 2),
    ("((A:0.1,B:0.2):0.05,(C:0.05,D:0.07):0.1);", 4),
    ("(A:0.1,B:0.2,C:0.3,D:0.4,E:0.5);", 5),
    ("((A:0.0,B:0.0)100/100:0.0,(C:0.0,D:0.0)100/100:0.0);", 4),
    ("((A:0.1,B:0.2)100/100:0.05,(Reference:0.05,D:0.07)88/90:0.1);", 3),
])
def test_newick_edge_cases(tmp_path, panel, newick, tips):
    _, status = render(tmp_path, panel, newick=newick)
    assert status["svg_status"] == "success", status
    assert status["tips_rendered"] == str(tips)


def test_iqtree_fallback_tree_is_reported_as_failure(tmp_path, panel):
    out, status = render(tmp_path, panel, newick="(IQTREE_failed:0.0);")
    assert status["svg_status"] == "failed"
    assert status["png_status"] == "failed"
    assert "fallback Newick" in status["error"]
    assert "not rendered" in (out / "phylogenetic_tree.svg").read_text()
    log = (out / "render.log").read_text()
    assert "Traceback" in log


def test_three_hundred_tips(tmp_path, panel):
    rnd = __import__("random").Random(1)
    nodes = [f"S{i}:{rnd.random() * 0.01:.5f}" for i in range(300)]
    while len(nodes) > 1:
        a = nodes.pop(rnd.randrange(len(nodes)))
        b = nodes.pop(rnd.randrange(len(nodes)))
        nodes.append(f"({a},{b}){rnd.randint(40, 100)}/{rnd.randint(40, 100)}:{rnd.random() * 0.01:.5f}")
    _, status = render(tmp_path, panel, newick=nodes[0] + ";")
    assert status["svg_status"] == "success"
    assert status["tips_rendered"] == "300"
