# rMAP-TB vs MTBseq — frozen baseline and phased benchmark audit
**Audit date:** 2026-11-03  
**Panel:** 18 synthetic libraries: 15 `MIX_` between distinct lineages, one `CTRL_SAMELIN_` same-lineage mixture, two `PURE_` controls.  
**Analyzed here:** previously uploaded original MTBseq v1.1.0 and rMAP-TB results. No cloud execution was performed for this audit.

## 1. Baseline frozen and recorded, not modified

Original MTBseq submission: `d11eb528-1d56-496c-901b-3213d38af7de`; workflow `e76bf8e8-bc07-46c5-a8d6-0a390c1cc569`. Workflow `MTBSEQ_BENCHMARK`, task `MTBSEQ_FULL`, image `quay.io/biocontainers/mtbseq:1.1.0--hdfd78af_1`.

**Exact logged parameters:** `--step TBfull`, `--project rMAPTB_benchmark`, `--ref M._tuberculosis_H37Rv_2015-11-13`, `--threads 8`, `--minbqual 13`, `--mincovf 4`, `--mincovr 4`, `--minphred 4`, `--minfreq 75`, `--unambig 95`, `--window 12`, `--distance 12`.

The entire original run archive is **not** contained in this small package; it remains on the user's Mac (`~/Downloads/MTBSEQ/mtbseq_out/mtbseq_outputs.tar.gz`, 21,598,164,132 bytes from screenshot/Terminal) and in the original Terra GCS output. Uploaded original logs and relevant upstream and downstream tables are included, with SHA-256 checksums and metadata in `baseline_original_outputs/`. **To preserve the 21.6 GB local archive**, run `scripts/preserve_original_terra_run.sh` (does not duplicate 21.6 GB). Cloud object immutability/version retention has not been enabled or independently confirmed.

## 2. Strict original 18-sample filtering (fully reproduced)

| Original MTBseq stage | Candidate positions |
|:--|--:|
| Amended joint candidate table | 3,148 |
| All 18 have unambiguous base at site | 291 |
| Official phylogenetic eligibility filter | 246 |
| Official 12-base window filter | 237 |
| Variable positions among final 18 sequences | **4** |

With 18 samples and `--unambig 95`, the minimum passing number is 18 (17/18 = 94.44%). A further 45 positions were excluded by official TBamend phylogenetic eligibility, and nine by the window filter. Of the last 237, 233 were constant among the libraries.

All **153/153** original pairwise distances were independently reconstructed from the original 237-position filtered FASTA: **108 comparisons at 0 and 45 at 4**. This validates *internal consistency* of the supplied MTBseq final output, not biological accuracy. The actual site filtering is the primary reason most genomic differences are absent from the matrix; it is not a demonstrated failure of the underlying MTBseq variant caller.

Reference: [MTBseq original manual](https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md), sections TBjoin/TBamend/TBgroups.

## 3. Pure-control analysis (post-hoc, not an independent MTBseq rerun)

Pure controls: `PURE_ERR038255_100x` (rMAP lineage L4) and `PURE_ERR13260175_100x` (L2).

Applying the original run's per-sample genotype requirements to **these two controls only** at joint candidate sites produced **2,201** callable SNP/reference positions. Of those, **1,278** carried different A/C/G/T genotype calls. We enforce both genotype type (`none` or `SNP`) and allele type (`ACGT`), with forward/reverse counts >=4, quality-20 read count >=4, and major allele frequency >=75%. Compared with the less restrictive per-allele-only count of 2,226, 25 positions of other variant call types are rightly excluded; the number of differing acceptable SNP/wild-type calls remains 1,278.

All **1,278** of those variable positions are absent from the 18-library `_u95_phylo` alignment. That proves the 0-SNP pure-control distance in the original 18-library final matrix is an artefact of the **all-library site mask**.

The rMAP-TB full-core alignment reports **1,148 differences over 3,944,469 pairwise-comparable genomic positions** for these two pure libraries. The 1,278 and 1,148 use different regions/quality masks: **do not subtract, divide, or use either as the other's ground truth**. This analysis demonstrates rescue of informative SNP calls when the other 16 libraries no longer determine callability. It is **not** an independent new MTBseq TBjoin run.

To obtain an actual native pure-only MTBseq result without rerunning BWA/GATK, the package contains: `prepare_pure_only_from_tar.py` (extracts just the two controls' original Called and Position_Tables files from the local 21.6 GB archive) and `run_pure_only_mtbseq.sh` (runs the official `TBjoin`, `TBamend`, and `TBgroups` stages in the same MTBseq 1.1.0 Biocontainer, with the baseline settings). Those commands require the Mac's local archive and a functioning Docker/Colima runtime and **have not been executed here**. `TBjoin` still depends on the variant and position table inputs being present and correctly named; the extraction script checks this.

**Statistical limitation:** two pure controls constitute only one independent pure–pure pair; no diagnostic accuracy, sensitivity, or specificity can be claimed for pure-strain SNP distances from this n=1 pair without a broader pure-only panel and validated truth set.

## 4. Mixed-infection recognition separated from SNP distance

**Sample-label truth:** 15 distinct-lineage `MIX_` samples, 1 same-lineage mixture labeled `CTRL_SAMELIN_`, 2 pure-strain controls.

**rMAP-TB/TB-Profiler output:** all 15 `MIX_` libraries are flagged `Possible mixed infection` (15/15); 14/15 have explicit two-lineage grouping (L2+L4 or L3+L4), while `MIX_ERR038255_ERR13260175_50pctB_100x` has only L2 as its major-lineage assignment despite the mixed-infection flag. Both pure controls are unflagged. The same-lineage control is also unflagged, which is expected as a challenge for lineage-based identification; report it separately, not as a simple clinical false negative. These are **descriptive recovery fractions on synthetic samples**, not generalizable diagnostic sensitivity or specificity.

**MTBseq output:** `mtbseq_classification.tsv` supplies a *single* Coll/Homolka lineage per library; `mtbseq_resistance.tsv` is empty of annotated resistance variants in the supplied extraction. Neither table is an explicit binary mixed-infection classifier. **Do not score MTBseq mixed-infection sensitivity from the single-lineage assignment or the zero resistance rows**.

**Orthogonal allele-frequency check for the 12 ERR038255 + ERR13260175 mixtures:** using 1,278 sites at which the corresponding pure parents differ and are high-confidence, examine the supplied MTBseq joint-table `Freq` (reported dominant allele frequency). Across positions with forward/reverse >=4, Qual20 >=4 and depth >=20, the descriptive proxy `100 - dominant_allele_frequency` tracks the nominal mixture labels: median estimates range from ~4.85% (5% nominal) to ~44.21% (50% nominal); the across-12 absolute error of medians is ~0.98 percentage points. This **is not direct alternate-allele read counting**, is not a validated minority-variant caller, and the meaning of A/B orientation is unverified. Restriction: the `ERR038266 + ERR13260106` mixtures lack both parental *pure* controls in the 18-sample panel; no equivalent parental-frequency score is reported for them. The same-lineage control is not scored on the L2/L4 marker panel.

See `analysis/mixture_detection_and_marker_frequencies.tsv` for every library, nominal percentage, reporting distinctions, and marker metrics.

## 5. Matched callable-site rMAP-TB vs MTBseq comparison — prepared, not executed

**What we can verify now:** source samples are matched across the 18 libraries; for 153 pairs the original filtered 237-site MTBseq distances and the rMAP whole-core distances are reproduced/recorded. They are **descriptive only**, because their denominators and mask policies differ.

**What is still missing:** the 18-library rMAP-TB `SNIPPY_CORE_MTBC/core.full.aln` is **not** part of the uploaded files. An unrelated 10-isolate alignment from UShER is NOT a substitute. Without the actual reference-coordinate genomic alignment for these same 18 libraries, matched-site agreement and any fair SNP-distance accuracy or clustering comparison cannot be calculated.

The package supplies:
- `scripts/fetch_rmap18_alignment.sh` to verify and download the expected 18-library GCS object to your Mac;
- `scripts/compare_matched_callable_sites.py` to verify all 18 sample IDs, reference coordinate length, and reference-base identity and then compare SNP distances and genotype concordance on **exactly the same positions and mutually callable alleles**, under (i) the original 237 sites, (ii) original joint-table candidate positions with zero insertion index and matching H37Rv coordinates (up to 2,864 of 3,148 rows), and (iii) the two-pure-parent candidate panel. Its fail-fast checks prevent pretending there is concordance when alignments do not share reference coordinates.

**Mandatory criterion for interpretation:** when all sites are jointly callable and the same reference coordinate scheme is established, compare per-pair distances; do not claim either caller is more accurate solely because one produces more SNPs. Mixed-strain samples should be marked **uninterpretable for direct transmission-clustering calls**. Analyze pure isolate pairs separately. All original runs should remain intact.

## 6. Predefined benchmark endpoints and decision rules

| Endpoint | Analysis population | Output | Status |
|---|---|---|---|
| Original MTBseq filtering and SNP matrix reproducibility | 18 libraries | 3,148 → 291 → 246 → 237; 153 exact matrix values | **COMPLETE** |
| Pure-control candidate-site separation | 2 pure controls | 2,201 pass sites; 1,278 distinct calls | **COMPLETE, post hoc** |
| Native MTBseq pure-only pipeline | Same 2 pure controls | Official two-sample TBjoin/TBamend/TBgroups matrix | **PREPARED; NOT RUN** |
| Mixed-infection descriptive recognition by rMAP-TB | 15 labeled MIX, 1 same-lineage control, 2 pure | Flags and lineage groups (15/15, 14/15; controls separated) | **COMPLETE** |
| MTBseq diagnostic-frequency proxy | 12 mixtures with both pure parent controls | Marker-level major-allele frequency and nominal-proportion concordance | **COMPLETE, exploratory** |
| True per-site matched rMAP-TB vs MTBseq | 18 libraries | Pairwise differences under identical 237/joint/pure-panel callable coordinates | **BLOCKED: 18-sample rMAP full alignment** |
| Whole-core SNP/clustering accuracy | Validated pure-only panel | Agreement with independent truth using harmonized masks | **NOT ESTIMABLE FROM CURRENT MATERIALS** |

## 7. Source provenance and caution

Raw outputs in `baseline_original_outputs/` are copies of user-provided original files. SHA-256 hashes are recorded. The full 21.6-GB archive is not copied, and no cloud changes have been performed. The locally-generated tables, summary and scripts are new audit artifacts. Exact output count/flag assertions are in `scripts/reproduce_baseline_audit.py`.

The MTBseq implementation and filtering-stage interpretation are documented in the original project manual: <https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md>.
