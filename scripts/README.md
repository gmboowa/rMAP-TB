# Scripts

Helper programs for running the workflow and for rebuilding the two validation
panels. Python 3.9 or later; standard library only, with one exception noted
below. All are run from the repository root.

## Running the workflow

| Script | Purpose |
|---|---|
| `make_inputs_json.py` | Creates a Cromwell inputs JSON from a directory of paired FASTQ files, detecting R1/R2 pairs and writing absolute paths |
| `validate_inputs.py` | Checks an inputs JSON against the workflow's declared inputs: unknown keys, types, pairing, file existence, gzip integrity, mask presence |
| `check_wdl.py` | Static checks on `rMAP_TB.wdl`: Bash and Python syntax of every embedded program, Python 3.11 f-string compatibility, call and output wiring |
| `collect_outputs.py` | Gathers the principal outputs of a Cromwell run into `results/<run-id>/{report,tables,figures,alignments,qc,logs}` with a manifest |
| `summarize_report.py` | Prints or exports (JSON) the headline figures of an integrated HTML report: counts, routing, resistance and lineage distributions, tree status, annotated tips, close pairs, exclusions |
| `archive_report.py` | Copies a report into `docs/reports/` under the naming convention `rMAP-TB_<version>_run<stamp>_<n>-sample_<panel>_integrated_report.html`, creating or appending to `docs/reports/README.md` |

Typical sequence:

```bash
python3 scripts/check_wdl.py rMAP_TB.wdl
python3 scripts/make_inputs_json.py --fastq-dir /data/run42/fastq \
    --adapters resources/adapters/adapters.fa \
    --reference resources/reference/H37Rv_NC_000962.3.gb \
    --mask resources/mask/H37Rv_NC_000962.3_PE_PPE_repetitive.bed --output inputs.json
python3 scripts/validate_inputs.py inputs.json
java -jar cromwell.jar run rMAP_TB.wdl --inputs inputs.json | tee run.log
python3 scripts/collect_outputs.py --cromwell-log run.log --dest results/
python3 scripts/summarize_report.py results/<run-id>/report/integrated_tb_amr_mtbc_phylogenomics_report.html
```

The two reports published under `docs/reports/panel1-mixtures/` and
`docs/reports/panel2-mtbc-ntm/` were placed there by hand, under directory
names that match the panels as the manuscript numbers them, rather than through
`archive_report.py`. Use `archive_report.py` for new reports where its naming
convention is wanted.

## Rebuilding the validation panels

These four reproduce the panels used for the validation and benchmarking
reported in the manuscript. The protocols and the analysis scripts are in
`benchmark/`; see `benchmark/README.md`.

| Script | Purpose |
|---|---|
| `fetch_reads.py` | Downloads the public sequencing libraries of either panel from ENA and verifies each file against ENA's published MD5. `--panel 1` is the five mixture parents; `--panel 2` is the MTBC, NTM and contamination panel. Resumable; falls back to ENA's deterministic FTP paths when the portal API is unavailable |
| `make_mixtures.py` | Builds the Panel 1 in-silico two-strain mixtures at defined minor-strain fractions and depths, sampling by sequenced bases rather than read count, under a fixed random seed. Writes a manifest and a Terra data table |
| `find_validation_isolates.py` | Finds candidate ENA runs for Panel 2: resolves each species to a taxon id through ENA's taxonomy service, queries for paired Illumina WGS runs and ranks them by proximity to a target depth. Writes a candidate table for review; downloads nothing |
| `build_ntm_arm.py` | Assembles Panel 2 from the fetched isolates, builds the two MTBC + NTM contamination mixtures, and writes a manifest stating the expected result for every sample before the workflow is run |

```bash
# Panel 1: the in-silico mixtures
python3 scripts/fetch_reads.py --panel 1 --check            # sizes only
python3 scripts/fetch_reads.py --panel 1 --out-dir RAW
python3 scripts/make_mixtures.py --fastq-dir RAW --out-dir MIX --dry-run
python3 scripts/make_mixtures.py --fastq-dir RAW --out-dir MIX

# Panel 2: MTBC, NTM and contamination
python3 scripts/fetch_reads.py --panel 2 --out-dir RAW_PANEL2
python3 scripts/build_ntm_arm.py \
    --candidates ntm_mtbc_candidates.tsv \
    --ntm-dir RAW_PANEL2 --mtbc-dir RAW --mtbc-accession ERR038255
```

The accessions for both panels are built into `fetch_reads.py` and are the same
list given in Additional file 1, Table S1. `find_validation_isolates.py` is the
search that selected the Panel 2 accessions; it is kept so that the selection is
auditable, and does not need to be re-run to rebuild the panel.

## Dependencies

Everything here is Python standard library, with one exception: `fetch_reads.py`
invokes `curl` for resumable downloads. No scripts require network access except
`fetch_reads.py` and `find_validation_isolates.py`, which query ENA.
