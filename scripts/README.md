# Scripts

All scripts use only the Python standard library (Python 3.9 or later) and run from the repository root.

| Script | Purpose |
|---|---|
| `build_snp_mask.py` | Builds the PE/PPE and repetitive-region BED mask from the H37Rv GenBank annotation (see `resources/mask/README.md`) |
| `make_inputs_json.py` | Creates a Cromwell inputs JSON from a directory of paired FASTQ files, detecting R1/R2 pairs and writing absolute paths |
| `validate_inputs.py` | Checks an inputs JSON against the workflow's declared inputs: unknown keys, types, pairing, file existence, gzip integrity, mask presence |
| `check_wdl.py` | Static checks on `rMAP_TB.wdl`: Bash and Python syntax of every embedded program, Python 3.11 f-string compatibility, call and output wiring |
| `collect_outputs.py` | Gathers the principal outputs of a Cromwell run into `results/<run-id>/{report,tables,figures,alignments,qc,logs}` with a manifest |
| `summarize_report.py` | Prints or exports (JSON) the headline figures of an integrated HTML report: counts, routing, resistance and lineage distributions, tree status, annotated tips, close pairs, exclusions |
| `archive_report.py` | Copies a report into `docs/reports/` under the naming convention `rMAP-TB_<version>_run<stamp>_<n>-sample_<panel>_integrated_report.html` and updates `docs/reports/README.md` |

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
python3 scripts/archive_report.py results/<run-id>/report/integrated_tb_amr_mtbc_phylogenomics_report.html \
    --version v22 --panel NTM-contamination-panel --description "..."
```
