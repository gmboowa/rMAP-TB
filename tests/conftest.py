"""Shared fixtures: a small synthetic MTBC panel that exercises every report section.

Nothing here touches the network or Docker. The fixtures mimic the files the
upstream tasks produce (TB-Profiler summary, canonical resistance table,
IQ-TREE Newick with SH-aLRT/UFBoot labels, SNP matrix, routing table, ...)
closely enough that the embedded report and renderer code runs unmodified.
"""

from __future__ import annotations

import csv
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

SAMPLES = [
    "CONTAM_ERR038255_ERR12861768_10pctNTM_100x",
    "CONTAM_ERR038255_ERR12861768_50pctNTM_100x",
    "ERR1203067", "ERR14187598", "ERR14187918", "ERR266109", "ERR3804031",
    "ERR970410", "SRR25473322", "SRR29440878", "SRR3500411", "SRR35607801", "SRR9157804",
]
NTM_SAMPLES = ["DRR360723", "ERR115016", "SRR12599345"]
LINEAGES = ["lineage4.8.2", "lineage4.8.2", "lineage2.2.1", "lineage4.1.2", "lineage4.1.2", "La1",
            "lineage4", "lineage3", "lineage1.2.1", "lineage4.1", "lineage4.3", "lineage2.2", "M.canetti"]
CATEGORIES = ["Sensitive", "Monoresistance", "MDR-TB", "RR-TB", "Pre-XDR-TB", "Sensitive", "XDR-TB",
              "Hr-TB", "Sensitive", "Resistance not determined by TB-Profiler", "Polyresistance", "Sensitive", "Unknown"]
DRUGS = ["", "ethambutol", "isoniazid,rifampicin", "rifampicin", "rifampicin,isoniazid,levofloxacin", "",
         "rifampicin,isoniazid,moxifloxacin,bedaquiline", "isoniazid", "", "", "ethambutol,pyrazinamide", "", ""]

NEWICK = (
    "((((0000_CONTAM_ERR038255_ERR12861768_10pctNTM_100x:0.0001,0001_CONTAM_ERR038255_ERR12861768_50pctNTM_100x:0.00005)100/100:0.012,"
    "0002_ERR1203067:0.011)98.2/97:0.03,((0003_ERR14187598:0.009,0004_ERR14187918:0.0095)100/100:0.004,"
    "(0009_SRR29440878:0.013,0011_SRR35607801:0.024)88/91:0.002)75.1/64:0.03)66/58:0.02,"
    "(0005_ERR266109:0.35,Reference:0.30)100/100:0.01,((0006_ERR3804031:0.03,0007_ERR970410:0.028)100/100:0.005,"
    "(0008_SRR25473322:0.03,(0010_SRR3500411:0.027,0012_SRR9157804:0.04)45/51:0.003)92/96:0.004)100/100:0.01);\n"
)


def write_tsv(path: Path, header: list[str], rows: list[list]) -> Path:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar="\\")
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    return path


@pytest.fixture(scope="session")
def panel(tmp_path_factory) -> dict:
    """Synthetic 13-MTBC + 3-NTM panel on disk; returns a dict of paths."""
    d = tmp_path_factory.mktemp("panel")
    rnd = random.Random(7)
    p: dict[str, Path] = {"dir": d}

    p["tree"] = d / "tree.nwk"
    p["tree"].write_text(NEWICK)

    tb_rows, res_rows, mut_rows, nonsyn_rows = [], [], [], []
    for i, (s, lin, cat, drugs) in enumerate(zip(SAMPLES, LINEAGES, CATEGORIES, DRUGS)):
        main = lin.split(".")[0]
        key_muts = ", ".join(f"Rv{rnd.randint(1, 3999):04d} c.{rnd.randint(1, 2000)}C>T" for _ in range(rnd.choice([0, 3, 12])))
        tb_rows.append([f"{i:04d}_{s}", "Mycobacterium tuberculosis complex (supported by Kraken2 species typing; TB-Profiler lineage reported)",
                        main, lin, cat, drugs, "", key_muts, f"{s}.json", "YES", "MTBC supported by Kraken2",
                        "Possible mixed infection" if i == 1 else "No evidence of mixed infection",
                        "lineage fraction 0.82" if i == 1 else "lineage fraction 1.00", "ok"])
        res_rows.append([s, cat, drugs, "", key_muts, "ok"])
        for dr in (drugs.split(",") if drugs else []):
            mut_rows.append([f"{i:04d}_{s}", dr, rnd.choice(["rpoB", "katG", "gyrA", "embB"]), f"p.Ser{rnd.randint(100, 500)}Leu",
                             f"c.{rnd.randint(100, 1500)}C>T", "Assoc w R", dr])
        for _ in range(rnd.randint(0, 3)):
            nonsyn_rows.append([f"{i:04d}_{s}", rnd.choice(["rpoB", "katG", "pncA"]), "missense_variant",
                                f"p.Ala{rnd.randint(10, 900)}Val", f"c.{rnd.randint(10, 2500)}C>T", "product"])
    p["tb"] = write_tsv(d / "tbprofiler_summary.tsv",
                        ["sample", "species", "main_lineage", "sub_lineage", "dr_type", "resistant_drugs", "resistance_mutations",
                         "key_mutations", "json_file", "mtbc_selected", "mtbc_selection_reason", "mixed_infection_status",
                         "mixed_infection_evidence", "status"], tb_rows)
    p["res"] = write_tsv(d / "resistance_profile_summary.tsv",
                         ["sample_id", "resistance_profile", "resistant_drugs", "resistance_mutations", "key_mutations", "status"], res_rows)
    p["mut"] = write_tsv(d / "mutation_evidence.tsv", ["sample", "drug", "gene", "mutation", "change", "confidence", "evidence"], mut_rows)
    p["nonsyn"] = write_tsv(d / "nonsyn.tsv", ["sample", "gene", "effect", "aa_change", "nt_change", "product"], nonsyn_rows)

    # A TB-Profiler table with one field far above csv's default 131072-byte
    # limit, reproducing the failure seen on real data.
    big_rows = [r[:6] + ["Rv0001 c.1A>G, " * 12000] + r[7:] for r in tb_rows]
    p["tb_big"] = write_tsv(d / "tbprofiler_summary_bigfield.tsv",
                            ["sample", "species", "main_lineage", "sub_lineage", "dr_type", "resistant_drugs", "resistance_mutations",
                             "key_mutations", "json_file", "mtbc_selected", "mtbc_selection_reason", "mixed_infection_status",
                             "mixed_infection_evidence", "status"], big_rows)

    all_samples = SAMPLES + NTM_SAMPLES
    p["species"] = write_tsv(d / "species_typing.tsv", ["Sample_ID", "Species_Identified", "Evidence"],
                             [[s, "Mycobacterium tuberculosis" if s in SAMPLES else "Mycobacterium abscessus", "Top species-level assignment"] for s in all_samples])
    routing_hdr = ["Sample_ID", "Most_Probable_Mycobacterium_Species", "Mycobacteria_Category", "MTBC_Supported", "MTBC_Reads",
                   "MTBC_Percent", "Evidence", "Workflow_Decision", "Selection_Basis"]
    routing = [[s, "Mycobacterium tuberculosis" if s in SAMPLES else "Mycobacterium abscessus",
                "MTBC" if s in SAMPLES else "Non-MTBC Mycobacterium", "YES" if s in SAMPLES else "NO",
                800000 if s in SAMPLES else 120, "97.5" if s in SAMPLES else "0.05", "Kraken2",
                "Retained for MTBC workflow" if s in SAMPLES else "Reported and excluded", "Kraken2"] for s in all_samples]
    p["routing"] = write_tsv(d / "routing.tsv", routing_hdr, routing)
    p["ntm"] = write_tsv(d / "ntm.tsv", routing_hdr, [r for r in routing if r[2] != "MTBC"])
    p["routing_status"] = d / "routing_status.txt"
    p["routing_status"].write_text(f"{len(SAMPLES)} MTBC-supported samples retained; {len(NTM_SAMPLES)} non-MTBC samples excluded.\n")

    # Symmetric SNP matrix with one identical pair and one distant outlier.
    n = len(SAMPLES)
    mat = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dist = 0 if {i, j} == {0, 1} else (rnd.randint(14000, 17000) if 5 in (i, j) else rnd.randint(400, 2300))
            mat[i][j] = mat[j][i] = dist
    p["matrix"] = d / "snp_matrix.tsv"
    with p["matrix"].open("w") as fh:
        fh.write("snp-dists\t" + "\t".join(SAMPLES) + "\tReference\n")
        for i, s in enumerate(SAMPLES):
            fh.write(s + "\t" + "\t".join(str(v) for v in mat[i]) + "\t1500\n")
        fh.write("Reference\t" + "\t".join(["1500"] * n) + "\t0\n")
    pairs, clusters = [], []
    for i in range(n):
        for j in range(i + 1, n):
            dist = mat[i][j]
            interp = ("Genomically close; review epidemiological linkage" if dist <= 5 else
                      "Intermediate SNP distance; review metadata" if dist <= 12 else "Not clustered by SNP threshold")
            pairs.append([SAMPLES[i], SAMPLES[j], dist, 3800000, interp, "close" if dist <= 5 else "distant"])
            if dist <= 5:
                clusters.append([f"Cluster_{len(clusters) + 1}", SAMPLES[i], SAMPLES[j], dist, interp])
    p["pairs"] = write_tsv(d / "pairs.tsv", ["sample1", "sample2", "snp_distance", "comparable_sites", "interpretation", "cluster_class"], pairs)
    p["clusters"] = write_tsv(d / "clusters.tsv", ["cluster_id", "sample1", "sample2", "snp_distance", "interpretation"], clusters)

    p["depth"] = write_tsv(d / "mean_depth_summary.tsv", ["Sample", "MeanDepth"], [[s, round(rnd.uniform(30, 120), 2)] for s in SAMPLES])
    p["included"] = write_tsv(d / "included_in_iqtree.tsv", ["sample", "alignment_length", "acgt_count", "missing_count", "missing_fraction"],
                              [[s, 120000, 118000, 2000, 0.016] for s in SAMPLES])
    p["excluded_empty"] = write_tsv(d / "excluded_from_iqtree_empty.tsv",
                                    ["sample", "alignment_length", "acgt_count", "missing_count", "missing_fraction", "threshold", "reason", "exclusion_note"], [])
    p["coverage_qc"] = write_tsv(d / "coverage_qc.tsv", ["sample", "mean_depth", "breadth_10x", "coverage_qc", "reason"],
                                 [["SRR9157804", 16.52, 0.8412, "FAIL", "mean_depth<15 or breadth_10x<0.90; excluded from core-SNP alignment"],
                                  ["ERR1203067", 88.1, 0.99, "PASS", "mean depth and 10x breadth above thresholds"]])
    p["metadata"] = write_tsv(d / "sample_metadata.tsv", ["sample", "country", "collection_year"],
                              [[s, "Uganda", 2026] for s in all_samples])
    p["iqstatus"] = d / "iqtree_status.txt"
    p["iqstatus"].write_text("success\n")
    p["filtsum"] = d / "alignment_filtering_summary.txt"
    p["filtsum"].write_text("Alignment filtering summary: 0 samples excluded.\n")
    p["iqreport"] = d / "iqtree.report"
    p["iqreport"].write_text("IQ-TREE 2.3.4\nModel of substitution: GTR+G\n")
    return p
