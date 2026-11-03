# MTBseq × rMAP-TB | Reproducible 18-library benchmark package

This package preserves SHA-256-identified copies of the original 18-library MTBseq outputs uploaded to ChatGPT (except the 21.6 GB full archive), reproduces the 153 pairwise MTBseq distances, evaluates the two pure controls separately, and summarizes the mixture-detection results using sample-label ground truth.

**Start with:** `reports/BENCHMARK_ASSESSMENT.md` and `analysis/checked_summary.json`. Do not treat original 237-site MTBseq distances and rMAP-TB whole-core distances as matched-site measurements.

## To do remaining native MTBseq pure-only analysis on your Mac

1. Preserve checksum of existing archive (reads but does not copy it):

   `bash scripts/preserve_original_terra_run.sh "$HOME/Downloads/MTBSEQ/mtbseq_out"`

2. Selectively extract only the two pure controls' Called/Position_Tables files from the already downloaded 21.6 GB archive:

   `python3 scripts/prepare_pure_only_from_tar.py --archive "$HOME/Downloads/MTBSEQ/mtbseq_out/mtbseq_outputs.tar.gz" --output "$HOME/Downloads/MTBseq_pure_only"`

3. With Docker or Colima running, produce new native two-sample MTBseq group-distance output:

   `bash scripts/run_pure_only_mtbseq.sh "$HOME/Downloads/MTBseq_pure_only"`

   The pure-only result is independent of the mixed samples' **joint callability requirement**, but it reuses the original individual-sample variant/position calls instead of remapping reads. This is a new `TBjoin/TBamend/TBgroups` subset analysis, not a full reprocessing of the reads. **It has not been executed in this audit environment.**

## To finish matched-callability testing

1. Verify and download the rMAP-TB alignment from the **18-library** run (not the UShER 10-isolate run):

   `bash scripts/fetch_rmap18_alignment.sh`

2. Execute:

   `python3 scripts/compare_matched_callable_sites.py --rmap-fasta "$HOME/Downloads/MTBseq_rMAP_matched_sites/core.full.aln" --audit-zip baseline_original_outputs/MTBseq_Extra_Audit.zip --amend-zip baseline_original_outputs/MTBseq_Amend_Audit.zip --outdir "$HOME/Downloads/MTBseq_rMAP_matched_sites/analysis"`

3. The script refuses to proceed if the reference length, H37Rv coordinate matches or 18 sample names are inconsistent. When successful, it reports 153 comparisons on the **same jointly callable coordinates** for three predeclared site masks.

## Folder guide

- `baseline_original_outputs/` — original uploaded files, hashes and baseline metadata (not the full 21.6 GB tar archive)
- `analysis/` — reproducible filtering and pure-site calculations, 153-pair descriptive report, mixture evaluation
- `inputs/` — original two pure-control FASTQ URI manifest and MTBseq two-column sample list
- `scripts/` — preservation, selective extraction, subset joint-analysis, reference-matched recalculation
- `reports/` — limitations, scientific findings, interpretive rules

**Scientific restriction:** The uploaded MTBseq matrix's 0/4 differences do not support 18-sample transmission clusters. The pure-only 1,278 candidate differences and the rMAP 1,148 whole-core SNP differences have **different denominators**, so neither is a validated accuracy benchmark for the other.
