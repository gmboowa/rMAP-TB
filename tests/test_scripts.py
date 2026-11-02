"""Tests for the helper scripts in scripts/."""

import gzip
import json
import subprocess
import sys
from pathlib import Path

from wdl_blocks import REPO_ROOT, WDL_PATH

SCRIPTS = REPO_ROOT / "scripts"


def run(script, *args):
    return subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], capture_output=True, text=True)


def make_fastq_dir(tmp_path, names):
    d = tmp_path / "fastq"
    d.mkdir()
    for n in names:
        with gzip.open(d / n, "wt") as fh:
            fh.write("@r1\nACGT\n+\nIIII\n")
    return d


def resources(tmp_path):
    res = tmp_path / "res"
    res.mkdir()
    (res / "adapters.fa").write_text(">a\nACGT\n")
    (res / "H37Rv.gb").write_text("LOCUS NC_000962 4411532 bp\n//\n")
    (res / "mask.bed").write_text("NC_000962.3\t0\t100\n")
    return res


def test_make_inputs_json_pairs_samples_in_order(tmp_path):
    d = make_fastq_dir(tmp_path, ["B_R2_001.fastq.gz", "A_1.fq.gz", "B_R1_001.fastq.gz", "A_2.fq.gz", "C.R1.fastq.gz", "C.R2.fastq.gz"])
    res = resources(tmp_path)
    out = tmp_path / "inputs.json"
    r = run("make_inputs_json.py", "--fastq-dir", d, "--adapters", res / "adapters.fa", "--reference", res / "H37Rv.gb",
            "--mask", res / "mask.bed", "--set", "use_gubbins=true", "--set", "max_cpus=8", "--output", out)
    assert r.returncode == 0, r.stderr
    data = json.loads(out.read_text())
    assert [Path(p).name for p in data["rMAP_TB.read1"]] == ["A_1.fq.gz", "B_R1_001.fastq.gz", "C.R1.fastq.gz"]
    assert [Path(p).name for p in data["rMAP_TB.read2"]] == ["A_2.fq.gz", "B_R2_001.fastq.gz", "C.R2.fastq.gz"]
    assert data["rMAP_TB.use_gubbins"] is True and data["rMAP_TB.max_cpus"] == 8
    assert data["rMAP_TB.snp_mask_bed"].endswith("mask.bed")


def test_make_inputs_json_refuses_incomplete_pairs(tmp_path):
    d = make_fastq_dir(tmp_path, ["A_R1.fastq.gz", "A_R2.fastq.gz", "B_R1.fastq.gz"])
    res = resources(tmp_path)
    r = run("make_inputs_json.py", "--fastq-dir", d, "--adapters", res / "adapters.fa", "--reference", res / "H37Rv.gb", "--output", tmp_path / "x.json")
    assert r.returncode != 0
    assert "incomplete" in r.stderr
    r = run("make_inputs_json.py", "--fastq-dir", d, "--adapters", res / "adapters.fa", "--reference", res / "H37Rv.gb",
            "--output", tmp_path / "x.json", "--allow-incomplete")
    assert r.returncode == 0
    assert len(json.loads((tmp_path / "x.json").read_text())["rMAP_TB.read1"]) == 1


def test_validate_inputs_accepts_generated_file(tmp_path):
    d = make_fastq_dir(tmp_path, [f"S{i}_R{m}.fastq.gz" for i in range(3) for m in (1, 2)])
    res = resources(tmp_path)
    out = tmp_path / "inputs.json"
    assert run("make_inputs_json.py", "--fastq-dir", d, "--adapters", res / "adapters.fa", "--reference", res / "H37Rv.gb",
               "--mask", res / "mask.bed", "--output", out).returncode == 0
    r = run("validate_inputs.py", out, "--wdl", WDL_PATH)
    assert r.returncode == 0, r.stdout
    assert "0 error(s)" in r.stdout


def test_validate_inputs_reports_problems(tmp_path):
    res = resources(tmp_path)
    bad = {
        "rMAP_TB.read1": ["/nonexistent/A_R1.fastq.gz", "/nonexistent/B_R1.fastq.gz"],
        "rMAP_TB.read2": ["/nonexistent/A_R2.fastq.gz"],
        "rMAP_TB.adapters": str(res / "adapters.fa"),
        "rMAP_TB.mtbc_reference_genbank": str(res / "H37Rv.gb"),
        "rMAP_TB.use_gubbins": "yes",
        "rMAP_TB.not_a_parameter": 1,
        "other.read1": [],
    }
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(bad))
    r = run("validate_inputs.py", p, "--wdl", WDL_PATH)
    assert r.returncode == 1
    for msg in ("read1 has 2 files but read2 has 1", "unknown workflow input 'not_a_parameter'",
                "use_gubbins: expected Boolean", "not in the rMAP_TB. namespace", "file not found"):
        assert msg in r.stdout, msg


def test_example_inputs_validate_structurally():
    # Paths are placeholders, so only structural errors are disallowed.
    for name in ("inputs.example.json", "inputs.minimal.example.json"):
        r = run("validate_inputs.py", REPO_ROOT / "examples" / name, "--wdl", WDL_PATH)
        assert "unknown workflow input" not in r.stdout, name
        assert "expected" not in r.stdout, name


def test_collect_outputs_from_cromwell_log(tmp_path):
    exec_dir = tmp_path / "cromwell-executions"
    exec_dir.mkdir()
    files = {}
    for name in ("integrated_tb_amr_mtbc_phylogenomics_report.html", "tbprofiler_summary.tsv", "phylogenetic_tree.svg", "iqtree_status.txt"):
        f = exec_dir / name
        f.write_text("x")
        files[name] = str(f)
    outputs = {
        "outputs": {
            "rMAP_TB.final_merged_html_report": files["integrated_tb_amr_mtbc_phylogenomics_report.html"],
            "rMAP_TB.tbprofiler_summary_tsv": files["tbprofiler_summary.tsv"],
            "rMAP_TB.tree_svg": files["phylogenetic_tree.svg"],
            "rMAP_TB.iqtree_status": files["iqtree_status.txt"],
            "rMAP_TB.gubbins_status": None,
            "rMAP_TB.fastqc_reports": [files["tbprofiler_summary.tsv"]],
        },
        "id": "abc123",
    }
    log = tmp_path / "run.log"
    log.write_text("[info] workflow started\n" + json.dumps(outputs, indent=2) + "\n[info] done\n")
    dest = tmp_path / "results"
    r = run("collect_outputs.py", "--cromwell-log", log, "--dest", dest)
    assert r.returncode == 0, r.stdout + r.stderr
    assert (dest / "abc123" / "report" / "integrated_tb_amr_mtbc_phylogenomics_report.html").is_file()
    assert (dest / "abc123" / "tables" / "tbprofiler_summary.tsv").is_file()
    assert (dest / "abc123" / "figures" / "phylogenetic_tree.svg").is_file()
    assert (dest / "abc123" / "logs" / "iqtree_status.txt").is_file()
    assert "manifest.tsv" in {p.name for p in (dest / "abc123").iterdir()}


def test_check_wdl_runs_clean():
    r = run("check_wdl.py", WDL_PATH)
    assert r.returncode == 0, r.stdout
