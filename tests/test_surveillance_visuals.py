"""Runs the TB_SURVEILLANCE_SUMMARY_VISUALS Python program on the synthetic panel."""

import re

from wdl_blocks import python_heredoc, run_python


def run_visuals(tmp_path, panel, matrix=None, metadata=None):
    (tmp_path / "surveillance_summary").mkdir(exist_ok=True)
    src = python_heredoc("TB_SURVEILLANCE_SUMMARY_VISUALS")
    args = [str(panel["tb"]), str(panel["species"]), str(matrix if matrix is not None else panel["matrix"]),
            str(panel["depth"]), str(metadata if metadata is not None else panel["metadata"])]
    r = run_python(src, tmp_path, args)
    assert r.returncode == 0, r.stderr
    return tmp_path / "surveillance_summary"


def test_heatmap_has_no_title_and_is_not_clipped(tmp_path, panel):
    out = run_visuals(tmp_path, panel)
    svg = (out / "snp_distance_heatmap.svg").read_text()
    assert svg.startswith("<svg")
    assert "Pairwise SNP Distance Heatmap" not in svg
    assert "SNP distance (log scale)" in svg
    assert 'id="snpHeatGrad"' in svg
    # Every sample appears as a row label and the long names fit in the left margin.
    longest = max(len(s) for s in re.findall(r'text-anchor="end" font-size="\d+" fill="#0f172a">([^<]+)</text>', svg))
    width = int(re.search(r'<svg[^>]*\swidth="(\d+)"', svg).group(1))
    left = min(float(x) for x in re.findall(r'<rect x="([\d.]+)" y="[\d.]+" width="\d+" height="\d+" fill="#', svg))
    assert left >= longest * 6.5
    assert width > left
    # Diagonal cells are blanked; off-diagonal cells carry tooltips.
    assert svg.count('fill="#f1f5f9"') == 13
    assert "<title>" in svg


def test_heatmap_without_matrix_explains_itself(tmp_path, panel):
    empty = tmp_path / "empty_matrix.tsv"
    empty.write_text("")
    out = run_visuals(tmp_path, panel, matrix=empty)
    svg = (out / "snp_distance_heatmap.svg").read_text()
    assert "No SNP distance matrix available" in svg


def test_metadata_and_rationale_tables(tmp_path, panel):
    out = run_visuals(tmp_path, panel)
    meta = (out / "tb_surveillance_metadata.tsv").read_text().splitlines()
    header = meta[0].split("\t")
    assert "sample" in header and "resistance_profile" in header
    assert any(h.startswith("user_") for h in header), "user metadata columns should be carried with a user_ prefix"
    assert len(meta) > 1
    rationale = (out / "qc_filtering_rationale.tsv").read_text().splitlines()
    assert rationale[0].split("\t")[0] == "sample"
    assert len(rationale) > 1
    lineage = (out / "lineage_distribution.tsv").read_text()
    assert "lineage" in lineage.splitlines()[0]
