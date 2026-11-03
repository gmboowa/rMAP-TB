# Benchmarking rMAP-TB against MTBseq, MAGMA and UShER

This directory holds everything needed to reproduce the head-to-head comparison
reported in the rMAP-TB manuscript (Table 7, Table 8 and Figure 6), and the
scripts that constructed the two validation panels.

Every number in the manuscript's benchmarking section is computed here, from
each workflow's own primary output files. No value is taken from any summary
document shipped alongside those outputs.

## Reproducing the comparison

Python 3.9 or later. `compare.py`, `pairs.py` and `tables.py` need only the
standard library; `fig_bench.py` additionally needs matplotlib.

```bash
cd benchmark
bash run_all.sh
```

or, equivalently:

```bash
python3 -I scripts/compare.py    data results   # -> results/comparison.json
python3 -I scripts/pairs.py      data results   # -> results/pairs.json
python3 -I scripts/tables.py     data results   # -> results/table_per_library.md
python3 -I scripts/fig_bench.py  results        # -> results/Figure6_benchmark.{png,pdf}
```

The files already in `results/` are the ones the manuscript reports. Re-running
overwrites them with identical content.

## Layout

```
benchmark/
├── README.md
├── run_all.sh
├── scripts/
│   ├── compare.py        per-library comparison across the four workflows
│   ├── pairs.py          pair-level analysis; the three unmixed control pairs
│   ├── tables.py         the per-library table (manuscript Table 8)
│   └── fig_bench.py      Figure 6, three panels
├── data/                 primary output files of the four workflows
│   ├── rmap/             rMAP-TB v22
│   ├── mtbseq/           MTBseq, 18-library cohort and the 2-library control run
│   ├── magma/            MAGMA
│   └── usher/            UShER
├── results/              outputs of the scripts above, as reported
└── panels/               how the two validation panels were built
    ├── panel1_mixtures/
    └── panel2_mtbc_ntm/
```

### What is in `data/`

The text outputs the scripts read: distance matrices, cluster and group
assignments, per-library statistics, classification tables, the MTBseq joint and
amended position tables, the alignment FASTAs, and the UShER audit file. Large
binary intermediates (BAM files, compressed archives, spreadsheets) are not
included, because no reported value is derived from them. The retained set
reproduces `results/` byte for byte.

The MTBseq tree contains two runs that differ only in cohort composition:

- the 18-library cohort, as submitted; and
- `MTBseq_pure_only/`, the same two unmixed parent libraries run alone at
  identical parameters and the same MTBseq version.

That pair is the internal control for the cohort-coupling result discussed in
the manuscript, and the two runs are what `compare.py` contrasts.

### What is not in `data/`

The FASTQ inputs, the BAM files and the full MTBseq working directories are too
large for version control. The FASTQs are public and can be fetched with
`scripts/fetch_reads.py` in the repository root `scripts/` directory; the
mixtures are then rebuilt deterministically by `panels/panel1_mixtures/`. The
complete, untrimmed output archives of all four workflows are deposited
alongside the tagged release.

## Rebuilding the panels

### Panel 1 — in-silico two-strain mixtures

```bash
python3 ../scripts/fetch_reads.py --panel 1 --out-dir RAW
python3 ../scripts/make_mixtures.py --fastq-dir RAW --out-dir MIX
```

Five parent libraries are downloaded and verified against ENA's published MD5s,
then combined at designed minor-strain fractions. Reads are sampled by
sequenced bases rather than by read count, because the parent libraries have
different read lengths and the designed fraction is a fraction of genome
coverage. Sampling uses a fixed random seed (default 20261101), so the panel is
reproducible to the read. `MIXTURE_PROTOCOL.md` states the design; `STEPS.md` is
the operational record of the run reported in the manuscript.

### Panel 2 — MTBC, NTM and contamination

```bash
python3 ../scripts/fetch_reads.py --panel 2 --out-dir RAW_PANEL2
python3 ../scripts/build_ntm_arm.py \
    --candidates ntm_mtbc_candidates.tsv \
    --ntm-dir RAW_PANEL2 \
    --mtbc-dir RAW \
    --mtbc-accession ERR038255
```

The pure MTBC and NTM libraries need no processing; this step builds only the
two contamination mixtures and writes `ntm_arm_manifest.tsv`, which states the
expected result for every sample before the workflow is run. `--mtbc-dir`
points at the Panel 1 download, so the *M. tuberculosis* side of the mixtures
costs no additional fetch. `--no-contamination` skips the mixtures.

`find_validation_isolates.py` is the ENA/SRA search that selected the 18
accessions and wrote `ntm_mtbc_candidates.tsv`. It is retained so that the
selection is auditable rather than asserted; the accessions it chose are
hard-coded in `fetch_reads.py --panel 2`, so it does not need to be re-run to
rebuild the panel.

## A note on UShER

UShER was run as a clustering method on the core alignment produced by rMAP-TB,
converted to VCF, not as an independent end-to-end pipeline. Its branch-mutation
counts are not equivalent to core-SNP counts. Both points are recorded in
`data/usher/.../audit_summary.txt` and are stated wherever its results appear in
the manuscript.

## Licence

MIT, as for the rest of rMAP-TB.

## Supplementary figures

`scripts/gen_states.py` rebuilds the three report states shown in Additional
file 2, Figure S1. It calls the same fixtures and the same `build_report()`
helper as `tests/test_merge_report.py`, extracting the report program from
`rMAP_TB.wdl` exactly as the tests do, so the HTML it writes is what the
test-suite asserts; it re-checks each test's assertions before writing. A
minimal `pytest` stand-in under `scripts/_pytest_stub/` lets those fixtures run
without installing pytest. With pytest installed, `pytest tests/` checks the
same three states.

```bash
python3 scripts/gen_states.py results          # the three report states, as HTML
python3 scripts/fig_s1.py     results          # Figure S1 (needs playwright)
python3 scripts/fig_s2.py     data/rmap/panel2/panel2_report_tables.json results
python3 scripts/build_af2.py  results          # Additional file 2, both pages
```

`fig_s1.py` screenshots each state with headless Chromium through Playwright;
`fig_s2.py` and `build_af2.py` need only matplotlib. Figure S2 is drawn from
the Panel 2 pairwise distance table in
`data/rmap/panel2/panel2_report_tables.json`, which is the report's own table
data, so the heatmap and Additional file 1 Table S4 come from one source.
