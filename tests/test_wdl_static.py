"""Static checks on rMAP_TB.wdl: syntax of embedded programs and call wiring."""

import re
import subprocess
import sys

from wdl_blocks import REPO_ROOT, WDL_PATH, wdl_text, task_body

EXPECTED_TASKS = [
    "PREPARE_READ_PAIRS", "TRIMMING", "FASTQC", "MULTIQC", "SPECIES_TYPING", "MYCOBACTERIA_SAMPLE_ROUTER",
    "TB_PROFILER_AND_MTBC_FILTER", "SNIPPY_CORE_MTBC", "TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS",
    "SNP_DISTANCE_CLUSTERING", "GUBBINS_RECOMBINATION", "IQTREE2_PHYLOGENY", "TREE_VISUALIZATION",
    "TB_SURVEILLANCE_SUMMARY_VISUALS", "MERGE_TB_REPORTS",
]


def test_wdl_exists_and_is_version_1():
    assert WDL_PATH.is_file()
    assert wdl_text().lstrip().startswith("version 1.0")


def test_all_tasks_present():
    text = wdl_text()
    found = re.findall(r"^task\s+(\w+)", text, re.M)
    assert found == EXPECTED_TASKS


def test_check_wdl_script_passes():
    r = subprocess.run([sys.executable, str(REPO_ROOT / "scripts" / "check_wdl.py"), str(WDL_PATH)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "0 problem(s)" in r.stdout


def test_check_wdl_script_detects_broken_wiring(tmp_path):
    broken = wdl_text().replace("tree_render_status = TREE_VISUALIZATION.render_status,",
                                "tree_render_status = TREE_VISUALIZATION.no_such_output,")
    p = tmp_path / "broken.wdl"
    p.write_text(broken)
    r = subprocess.run([sys.executable, str(REPO_ROOT / "scripts" / "check_wdl.py"), str(p)], capture_output=True, text=True)
    assert r.returncode == 1
    assert "no_such_output" in r.stdout


def test_tree_visualization_outputs_declared():
    body = task_body("TREE_VISUALIZATION")
    for name in ("tree_image", "tree_svg", "cleaned_tree", "render_log", "render_status", "excluded_from_iqtree"):
        assert re.search(rf"File\s+{name}\s*=", body), name


def test_tree_visualization_raises_csv_field_limit():
    # Regression: TB-Profiler fields above 131072 bytes made csv.DictReader raise
    # and the tree lost all lineage/resistance annotations.
    assert "csv.field_size_limit" in task_body("TREE_VISUALIZATION")
    assert "csv.field_size_limit" in task_body("MERGE_TB_REPORTS")


def test_merge_receives_tree_diagnostics_and_coverage_gate():
    text = wdl_text()
    wf = text.split("task PREPARE_READ_PAIRS")[0]
    for line in ("tree_svg = TREE_VISUALIZATION.tree_svg", "tree_render_log = TREE_VISUALIZATION.render_log",
                 "tree_render_status = TREE_VISUALIZATION.render_status", "coverage_qc_tsv = SNIPPY_CORE_MTBC.coverage_qc_tsv"):
        assert line in wf, line


def test_no_one_by_one_placeholder_png():
    # The 68-byte transparent PNG must never come back.
    assert "iVBORw0KGgoAAAANSUhEUgAAAAEAAAAB" not in task_body("TREE_VISUALIZATION")


def test_heatmap_has_no_title_text():
    body = task_body("TB_SURVEILLANCE_SUMMARY_VISUALS")
    assert "Pairwise SNP Distance Heatmap</text>" not in body
