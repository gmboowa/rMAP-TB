# NTM and non-tuberculosis MTBC validation arm — step by step

Run everything in your own Terminal. `/usr/bin/python3` is enough; no installs.
Work in a new folder so this arm never mixes with the 18 mixtures:

```bash
mkdir -p ~/ntm_validation && cd ~/ntm_validation
# copy find_validation_isolates.py, fetch_reads.py and build_ntm_arm.py here
```

---

## What this arm tests, and why each species is in it

Two claims in the manuscript have never been tested with a single sample:
that rMAP-TB discriminates MTBC from NTM, and that its species layer has
zoonotic relevance. This panel tests both, and is deliberately built so that
a failure tells you *which* component failed.

| Group | Species | In the Kraken2 DB? | What a correct result looks like |
|---|---|---|---|
| `mtbc_other` | *M. bovis* ×2 | no | routed **MTBC**; TB-Profiler La1 or M. bovis; Kraken2 says *M. tuberculosis* |
| `mtbc_other` | *M. africanum* ×2 | no | routed **MTBC**; TB-Profiler **L5 or L6**; Kraken2 says *M. tuberculosis* |
| `mtbc_other` | *M. canettii* ×1 | **yes** | routed **MTBC**; Kraken2 **names it correctly** |
| `ntm_in_db` | *M. avium*, *M. intracellulare*, *M. kansasii* | yes | routed **NTM**, named correctly, TB-Profiler skipped |
| `ntm_absent` | *M. abscessus* ×2, *M. fortuitum* ×1 | **no** | **unknown — this is the measurement** |
| `contamination` | *M. tuberculosis* + NTM at 10 % and 50 % | — | stays MTBC, species and lineage unchanged |

*M. canettii* is the positive control for species naming: it is the only
non-tuberculosis MTBC member in your database, so if Kraken2 cannot name it,
the problem is the database build and not the missing genomes.

The `ntm_absent` group is the one worth watching. Your router takes the top
species-level row by read count with **no minimum reads and no minimum
percentage**, so an organism absent from the database may be given a
confident but wrong species name from a handful of spurious reads. Whether
that happens, and how often, is the result. The category should still come
out NTM, so routing may well survive even when the name is wrong — record
both separately.

---

## Step 1 — find candidates (a minute, downloads nothing)

```bash
python3 find_validation_isolates.py
```

For each species it resolves the taxon id through ENA's taxonomy service —
nothing is hard-coded, so no accession is guessed — then queries the portal
for paired Illumina WGS runs, drops anything below 60× or above 3 GB, and
ranks what is left by closeness to 150×. Rows marked `*` are auto-selected
and written to `ntm_mtbc_candidates.tsv`.

Read the output before moving on. If a species reports **SKIP**, try:

```bash
python3 find_validation_isolates.py --min-depth 40 --max-gb 5
python3 find_validation_isolates.py --subtaxa        # e.g. BCG under M. bovis
```

You can edit `ntm_mtbc_candidates.tsv` freely — swap an accession, delete a
row, add one. Everything downstream reads that file, so it is the single
place where the panel is decided.

## Step 2 — download and verify

```bash
python3 fetch_reads.py --from-table ntm_mtbc_candidates.tsv --out-dir RAW_NTM
```

Same fetcher as before: MD5-verified against ENA, resumable, safe to re-run.
Expect roughly 10–15 GB for eleven isolates. Use `--check` first if you want
the sizes without committing.

## Step 3 — build the arm

The pure isolates need no processing. This step builds only the two
contamination mixtures and writes the manifest and Terra table:

```bash
python3 build_ntm_arm.py \
  --candidates ntm_mtbc_candidates.tsv \
  --ntm-dir RAW_NTM \
  --mtbc-dir ~/CRyPTIC/RAW \
  --mtbc-accession ERR038255
```

`--mtbc-dir` points at the M. tuberculosis isolate you already downloaded, so
there is nothing new to fetch for the mixtures. The contaminant defaults to
the first `ntm_in_db` run; override with `--contaminant ERRxxxxxxx`.

Outputs: `ntm_arm_manifest.tsv` (the expected result for every sample) and
`terra_ntm_table.tsv`. Add `--no-contamination` to skip the mixtures.

## Step 4 — upload and rewrite the paths

This is the step that broke the last run. Do both halves.

```bash
BUCKET=gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_ntm

gsutil -m cp RAW_NTM/*.fastq.gz MIX_NTM/*.fastq.gz $BUCKET/

python3 - <<'PY'
import os
B = os.environ["BUCKET"].rstrip("/")
src, out = "terra_ntm_table.tsv", "terra_ntm_table_gs.tsv"
lines = open(src).read().splitlines()
with open(out, "w") as g:
    g.write(lines[0] + "\n")
    for l in lines[1:]:
        if not l.strip():
            continue
        s, r1, r2 = l.split("\t")
        g.write(f"{s}\t{B}/{os.path.basename(r1)}\t{B}/{os.path.basename(r2)}\n")
print("wrote", out)
PY

grep -c 'gs://' terra_ntm_table_gs.tsv      # must equal the number of samples
```

Import `terra_ntm_table_gs.tsv`, **not** `terra_ntm_table.tsv`. Confirm the
Data tab shows `gs://` paths before launching.

## Step 5 — run

Launch `rMAP_TBv11.wdl` with your existing `inputs.json`, unchanged. Build a
set from all rows and use the expression that matches whatever Terra named
the table: if the entity type is `tuberculosis_ntm`, the inputs are
`this.tuberculosis_ntms.read1` and `this.tuberculosis_ntms.read2`, unquoted.

Expect some samples to be routed away from MTBC analyses — that is the
correct behaviour, not a failure. A run where every NTM sample went through
TB-Profiler would mean routing is broken.

## Step 6 — score it

Join the integrated report to `ntm_arm_manifest.tsv` on sample name and
record, per sample:

1. the Kraken2 species call,
2. the routing category (MTBC / NTM / unresolved / non-*Mycobacterium*),
3. the TB-Profiler lineage, where one was run,
4. whether that matches the `expected` column.

The four reportable results are: **routing accuracy** (does every NTM go to
the NTM branch and every MTBC member stay MTBC), **species-naming accuracy
within the database**, **behaviour on organisms absent from the database**,
and **contamination tolerance** (does 10 % and 50 % NTM change the MTBC
call). Send me that table and I will work out the numbers and draft the text.

---

## If something breaks

- **A species reports SKIP at step 1** — relax `--min-depth` / `--max-gb`, or
  try `--subtaxa`. Do not substitute a species from a different group to make
  the panel look complete.
- **An `ntm_absent` sample gets a confident wrong species name** — that is a
  result, not a bug. Record it. It is the strongest argument for either
  rebuilding the database or adding a read-count floor to the species call.
- **An MTBC member is routed to NTM** — that is a real defect. Send me the
  routing TSV and the Kraken2 report for that sample.
- **Everything is "no species-level call"** — likely a database or report
  path problem rather than the samples. Send one Kraken2 report.
