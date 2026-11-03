# Mixture experiment — step by step

Run everything in your own Terminal, in `/Users/gmboowa/CRyPTIC`.
Nothing here needs conda, pip, or any install. `/usr/bin/python3` is enough.

Note: this folder is named CRyPTIC, but this experiment does **not** use CRyPTIC
data. It uses five isolates already in your benchmark set. Consider renaming the
folder to `mixture_validation` so the two experiments don't get confused later.

---

## Step 1 — see what you are about to download (30 seconds, downloads nothing)

```bash
cd /Users/gmboowa/CRyPTIC
python3 fetch_reads.py --check
```

Prints the five isolates with their sizes and read counts. Expect roughly
10–25 GB in total. You have about 150 GB free, which is enough.

If this errors with "could not reach ENA", you have a network problem — fix that
before going further. Nothing has been downloaded at this point.

## Step 2 — download and verify (the slow step)

```bash
python3 fetch_reads.py
```

Downloads into `RAW/` and checks every file against the MD5 ENA publishes. It is
resumable: if it drops out, re-run the same command and it picks up where it
stopped. Files that are already complete and correct are skipped, so re-running
is always safe.

Leave it going. On a normal connection this is the only part measured in hours.

## Step 3 — check the mixing plan before building (seconds)

```bash
python3 make_mixtures.py --fastq-dir RAW --out-dir MIX --dry-run
```

Read two things in the output:

1. **Available depth per isolate.** Every planned mixture must fit inside it.
   Anything that doesn't is marked SKIP, with the reason. If the primary pair
   gets skipped, stop and tell me — the design needs changing, not forcing.
2. **Read lengths.** The ERR038xxx isolates are old submissions and the
   ERR13260xxx ones are recent, so these will probably differ. That is fine and
   expected — mixing is done on sequenced bases, not read counts, so the minor
   fraction stays exact. Write down what it reports; it goes in the Methods.

## Step 4 — build the mixtures (minutes)

```bash
python3 make_mixtures.py --fastq-dir RAW --out-dir MIX
```

Writes 18 sample pairs plus two tables:

- `MIX/mixture_manifest.tsv` — the expected result for every sample. This is
  what you score the run against afterwards. Keep it.
- `MIX/terra_tuberculosis_table.tsv` — the Terra data table, with local paths.

The subsampling is seeded, so this is reproducible. Default seed 20261101 —
state it in the Methods.

## Step 5 — upload

```bash
# replace with your bucket path
BUCKET=gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures

gsutil -m cp MIX/*.fastq.gz $BUCKET/
python3 - <<'PY'
import os
B = os.environ["BUCKET"]
src = "MIX/terra_tuberculosis_table.tsv"
out = "MIX/terra_tuberculosis_table_gs.tsv"
with open(src) as f, open(out, "w") as g:
    g.write(next(f))
    for line in f:
        s, r1, r2 = line.rstrip("\n").split("\t")
        g.write(f"{s}\t{B}/{os.path.basename(r1)}\t{B}/{os.path.basename(r2)}\n")
print("wrote", out)
PY
```

`gsutil` was not on the PATH when I looked, so if it is missing, activate the
conda environment that has it, or install the Google Cloud SDK.

## Step 6 — run

1. Import `MIX/terra_tuberculosis_table_gs.tsv` as the `tuberculosis` table.
2. Build a `tuberculosis_set` from all 18 rows.
3. Launch the current `rMAP_TB.wdl` with your existing `inputs.json` — **no
   changes to either**. The mixtures are ordinary FASTQs.
4. Inputs stay unquoted: `this.tuberculosiss.read1`, `this.tuberculosiss.read2`.

A few hours on the 16 CPU / 64 GB configuration.

## Step 7 — score it

From the integrated report, pull `mixed_infection_status` and
`mixed_infection_evidence` for all 18 samples and join them to
`MIX/mixture_manifest.tsv` on sample name. Send me that join and I will work out
the detection limits and write the Results text.

The four numbers that come out of it are in `MIXTURE_PROTOCOL.md`, section 4.

---

## If something breaks

- **Step 2 stalls or fails repeatedly** — re-run it; it resumes. Verified files
  are never refetched.
- **"MD5 mismatch"** — the file is corrupt. Delete that one file and re-run.
- **Step 3 marks the primary pair SKIP** — do not lower the depths to force it.
  Tell me what it reported and I will revise the design.
- **Any mixture fails in the workflow** — send the task log. Do not drop the
  sample silently; a failure at a specific mixing fraction is itself a result.
