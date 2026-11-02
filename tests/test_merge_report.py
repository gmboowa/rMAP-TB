"""Runs the MERGE_TB_REPORTS Python program and checks what the report says.

The scenarios cover the failure modes the report must make visible: a tree
that rendered, a tree whose annotations failed to join, a renderer failure,
the legacy 68-byte placeholder image, and a sample removed by the coverage gate.
"""

import base64
import shutil
import sys
from pathlib import Path

import pytest

from wdl_blocks import python_heredoc, run_python, substitute

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from summarize_report import summarize  # noqa: E402

ONE_PIXEL_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/p9sAAAAASUVORK5CYII=")


@pytest.fixture(scope="module")
def rendered_tree(tmp_path_factory, panel):
    """A real SVG tree produced by the TREE_VISUALIZATION program."""
    cwd = tmp_path_factory.mktemp("tree")
    src = substitute(python_heredoc("TREE_VISUALIZATION"), {
        "image_format": "png", "input_tree": str(panel["tree"]), "tbprofiler_summary_tsv": str(panel["tb"]),
        "resistance_profile_summary_tsv": str(panel["res"]), "iqtree_excluded_samples_tsv": "",
        "width": "2400", "height": "1600", "max_bootstrap_to_label": "100",
        "title": "MTBC Core-SNP Phylogenetic Tree", "midpoint_root_tree": "true"})
    r = run_python(src, cwd)
    assert r.returncode == 0, r.stderr
    return cwd / "tree_visualization"


def build_report(tmp_path, panel, *, svg=None, png=None, status=None, log=None, coverage=None, excluded=None, newick=None):
    final = tmp_path / "final_report"
    final.mkdir(exist_ok=True)
    if svg:
        shutil.copy(svg, final / "mtbc_tree.svg")
    if png is not None:
        if isinstance(png, bytes):
            (final / "mtbc_tree.png").write_bytes(png)
        else:
            shutil.copy(png, final / "mtbc_tree.png")
    if status:
        (final / "tree_render_status.txt").write_text(status if isinstance(status, str) else status.read_text())
    if log:
        shutil.copy(log, final / "tree_render.log")
    if newick:
        shutil.copy(newick, final / "mtbc_tree.cleaned.nwk")
    src = python_heredoc("MERGE_TB_REPORTS")
    args = [str(panel["tb"]), str(panel["res"]), str(panel["species"]), str(panel["ntm"]), str(panel["routing"]),
            str(panel["routing_status"]), str(panel["nonsyn"]), str(panel["mut"]), str(panel["pairs"]), str(panel["clusters"]),
            "", "", "", "", "", "", str(panel["iqreport"]), "", "",
            str(excluded or panel["excluded_empty"]), str(panel["included"]), str(panel["filtsum"]), str(panel["iqstatus"]),
            "final_report/mtbc_tree.svg", "final_report/tree_render.log", "final_report/tree_render_status.txt",
            str(coverage) if coverage else ""]
    r = run_python(src, tmp_path, args)
    assert r.returncode == 0, r.stderr
    html_path = final / "integrated_tb_amr_mtbc_phylogenomics_report.html"
    return html_path, html_path.read_text(encoding="utf-8")


def test_successful_tree_is_embedded_as_svg(tmp_path, panel, rendered_tree):
    path, html = build_report(tmp_path, panel, svg=rendered_tree / "phylogenetic_tree.svg",
                              png=rendered_tree / "phylogenetic_tree.png",
                              status=rendered_tree / "render_status.txt", log=rendered_tree / "render.log",
                              newick=rendered_tree / "phylogenetic_tree.cleaned.nwk")
    s = summarize(path)
    assert s["tree_rendering"]["state"] == "Rendered"
    assert len(s["tree_tips"]) == 13 and s["tree_tips_unannotated"] == 0
    assert s["counts"]["total_paired_samples"] == 16
    assert s["counts"]["mtbc_isolates_retained"] == 13
    assert s["iqtree_status"] == "success"
    assert 'id="treeFigure"' in html and "<svg" in html.split('id="treeFigure"')[1][:500]
    assert "Download SVG (vector)" in html and "Download Newick" in html
    assert "Pairwise SNP Distance Heatmap" not in html
    assert s["resistance_distribution"]["MDR-TB"] == 1
    assert s["mixed_infection_flags"] == 1
    assert s["close_pairs"] == 1
    # Every section present, numbered and navigable.
    for sec in ("sec-qc", "sec-species", "sec-ntm", "sec-tbprofiler", "sec-mutation-evidence", "sec-nonsyn",
                "sec-lineage", "sec-snp-distance", "sec-heatmap", "sec-tree", "sec-surveillance", "sec-provenance"):
        assert f'id="{sec}"' in html, sec


def test_renderer_failure_is_shown_with_log(tmp_path, panel):
    status = ("svg_status=failed\npng_status=failed\npng_renderer=matplotlib-placeholder\ntips_rendered=0\n"
              "midpoint_rooted=false\nreference_tips_removed=0\nmetadata_status=ok\nmetadata_rows=13\ntips_matched=0\n"
              "error=ValueError('IQ-TREE did not produce a valid phylogenetic tree.')\n")
    log = tmp_path / "render.log"
    log.write_text("TREE_VISUALIZATION started\nERROR: ValueError(...)\nTraceback (most recent call last):\n")
    placeholder_svg = tmp_path / "placeholder.svg"
    placeholder_svg.write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>Phylogenetic tree was not rendered</text></svg>')
    path, html = build_report(tmp_path, panel, svg=placeholder_svg, status=status, log=log)
    s = summarize(path)
    assert s["tree_rendering"]["state"] == "Rendering failed"
    assert "Phylogenetic tree image was not rendered" in html
    assert "Tree renderer log (render.log)" in html
    assert "IQ-TREE did not produce a valid phylogenetic tree" in html
    assert 'id="treeFigure"' not in html


def test_legacy_one_pixel_png_is_not_shown_as_a_tree(tmp_path, panel):
    path, html = build_report(tmp_path, panel, png=ONE_PIXEL_PNG)
    s = summarize(path)
    assert s["tree_rendering"]["state"] == "Rendering failed"
    assert "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAAB" not in html


def test_unannotated_tree_is_flagged(tmp_path, panel, rendered_tree):
    status = (rendered_tree / "render_status.txt").read_text()
    status = status.replace("metadata_status=ok", "metadata_status=failed: Error('field larger than field limit (131072)')")
    status = status.replace("tips_matched=13", "tips_matched=0")
    path, html = build_report(tmp_path, panel, svg=rendered_tree / "phylogenetic_tree.svg", status=status,
                              log=rendered_tree / "render.log")
    s = summarize(path)
    assert s["tree_rendering"]["state"] == "Rendered, annotations incomplete"
    assert "Tree annotations incomplete" in html
    assert "field larger than field limit" in html


def test_coverage_gate_exclusion_is_explained(tmp_path, panel, rendered_tree):
    path, html = build_report(tmp_path, panel, svg=rendered_tree / "phylogenetic_tree.svg",
                              status=rendered_tree / "render_status.txt", coverage=panel["coverage_qc"])
    assert "failed the post-mapping coverage gate" in html
    assert "SRR9157804" in html
    assert "mean depth 16.52" in html
    assert "1 MTBC-routed sample(s) excluded from the core-SNP alignment" in html


def test_iqtree_exclusion_notes(tmp_path, panel, rendered_tree):
    excluded = tmp_path / "excluded.tsv"
    excluded.write_text("sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note\n"
                        "SRR9157804\t120000\t30000\t90000\t0.75\t0.5\tmissing_fraction_ge_0.5\t\n")
    path, html = build_report(tmp_path, panel, svg=rendered_tree / "phylogenetic_tree.svg",
                              status=rendered_tree / "render_status.txt", excluded=excluded)
    assert "75.00% of its core-SNP alignment" in html
    s = summarize(path)
    assert any(e["sample"] == "SRR9157804" for e in s["excluded_samples"])


def test_report_is_self_contained(tmp_path, panel, rendered_tree):
    _, html = build_report(tmp_path, panel, svg=rendered_tree / "phylogenetic_tree.svg", status=rendered_tree / "render_status.txt")
    assert 'src="http' not in html and 'href="http' not in html.split("<body", 1)[0]
    assert "<link rel=\"stylesheet\"" not in html
    assert html.count("<script") == 1
