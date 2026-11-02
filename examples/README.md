# Example inputs

| File | Use |
|---|---|
| `inputs.minimal.example.json` | The four required inputs plus the SNP mask; every other parameter takes its WDL default. Start here. |
| `inputs.example.json` | Every workflow-level parameter listed with its default value, for runs that need explicit control. |
| `inputs.terra.example.json` | Terra/Dockstore configuration. `read1`/`read2` are taken from a sample set (`this.samples.<column>`); resource files are copies of `resources/` in the workspace bucket. |
| `sample_metadata.example.tsv` | Optional epidemiological metadata. The identifier column (`sample`, `sample_id`, `isolate`, `isolate_id`, `accession`, `run_accession` or `sra_accession`) must match the sample names derived from the FASTQ file names; all other columns are carried into the surveillance metadata table with a `user_` prefix. |

Generate an inputs file from a directory of FASTQ files and check it before running:

```bash
python3 scripts/make_inputs_json.py --fastq-dir /data/run42/fastq \
    --adapters resources/adapters/adapters.fa \
    --reference resources/reference/H37Rv_NC_000962.3.gb \
    --mask resources/mask/H37Rv_NC_000962.3_PE_PPE_repetitive.bed \
    --output inputs.json
python3 scripts/validate_inputs.py inputs.json
java -jar cromwell.jar run rMAP_TB.wdl --inputs inputs.json | tee run.log
python3 scripts/collect_outputs.py --cromwell-log run.log --dest results/
```

Paths in the JSON files must be absolute. Three samples is the minimum for the phylogenomics branch to start; at least three must also pass the coverage gate and the IQ-TREE alignment-quality filter for a tree to be inferred.
