# rMAP-TB mixed-infection validation — in-silico mixture protocol

Purpose: replace "no sensitivity or specificity is claimed for this flag" with a
measured detection limit, using only isolates already in the benchmark set.

## 1. Why these isolates

Chosen from the ten benchmark isolates using lineage, resistance profile and mean
depth from the verified run (20261101_053453_UTC).

| Role | Isolate | Lineage | Resistance | Mean depth |
|---|---|---|---|---|
| Primary, major | ERR038255 | L4.8.2 | Sensitive | 164.1× |
| Primary, minor | ERR13260175 | L2.2.1 | RR-TB (rifampicin) | 188.9× |
| Secondary, major | ERR038266 | L3.1.1 | Sensitive | 152.4× |
| Secondary, minor | ERR13260106 | L4.1.2.1 | RR-TB (rifampicin) | 236.7× |
| Control | ERR038254 + ERR038255 | both L4.8.2 | both sensitive | 151.3× / 164.1× |

The primary pair tests **both** arms of the flag at once: L4 against L2 is the most
divergent pairing available and the most clinically common mixture, so barcode SNPs
appear at the minor fraction; and because one isolate is susceptible and the other
RR-TB, the rifampicin allele also appears at the minor fraction, which is exactly
the heteroresistance condition.

The secondary pair (L3 against L4) is an independent replicate with a different
lineage distance, so a result is not an artefact of one pair.

The control pair is the important one. ERR038254 and ERR038255 are the same
sub-lineage and one SNP apart. Mixing them 50:50 **must not** raise the
mixed-infection flag. Without this, a positive result on the other mixtures is
uninterpretable, because a flag that fires on everything is not a detector.

Excluded: ERR13260124 (29× — too shallow to subsample) and ERR046825 (80× — would
cap the depth series).

## 2. Running it

```
python3 make_mixtures.py --fastq-dir RAW --out-dir MIX --dry-run   # check depths
python3 make_mixtures.py --fastq-dir RAW --out-dir MIX             # build
```

`RAW` holds the untrimmed paired FASTQs named `ERR038255_1.fastq.gz` and so on
(`_R1`/`_R2` and uncompressed also work). Mixing happens before trimming, so the
workflow is run completely unmodified — that is the point.

No dependencies. The subsampling is seeded (`--seed`, default 20261101), so the
mixtures are reproducible; state the seed in the Methods.

It writes 18 sample pairs, `mixture_manifest.tsv` with the expected result for each,
and `terra_tuberculosis_table.tsv`. Replace the local paths in that table with
`gs://` URIs after upload, then import it as the `tuberculosis` table and build a
`tuberculosis_set` from it.

## 3. The design

- **Minor-fraction ladder**, primary pair, 100× total: 5, 10, 20, 30, 50 % minor.
  This is the curve that gives you the detection limit.
- **Reciprocal ladder**, 5–30 % with the susceptible isolate as minor. Detection is
  not symmetric, and a reviewer will ask.
- **Depth series** at a fixed 10 % minor: 40, 70, 100, 140× total. This separates
  "the flag missed it because the fraction was low" from "because the depth was
  low" — two different claims, and only this series tells them apart.
- **Pure endpoints** at 100×, which must come back clean.
- **Secondary pair** at 10, 20, 50 %.
- **Same-lineage control** at 50:50, which must come back clean.

## 4. Reading the result

For each sample take `mixed_infection_status` and `mixed_infection_evidence` from
the TB-Profiler summary and compare with the `expected` column of the manifest.

The reportable numbers are:

1. **Lineage-mixture detection limit** — the lowest minor fraction at which
   "Possible mixed infection" is returned, at 100× and at each depth.
2. **Heteroresistance detection limit** — the lowest minor fraction at which the
   rifampicin variant is called in the 0.10–0.90 window.
3. **Specificity** — the pure samples and the same-lineage control must all return
   "No mixed-infection signal detected". If any of them does not, that is a false
   positive and it is the first thing to fix.
4. **Observed vs expected minor fraction** — plot the TB-Profiler lineage fraction
   against the known mixing fraction. If the relationship is linear, the flag is
   also quantitative and you can say so.

Expect the 5 % points to fail and the 30–50 % points to succeed. Where the
transition falls is the result. A limit around 10–20 % would be a respectable and
publishable finding for a flag built from fields TB-Profiler already reports; it
would be honestly below what fastlin and TBtypeR achieve, and saying so is
stronger than claiming nothing.

Also check what the mixtures do to **SNP distances**, since that is the argument
for having the flag at all. If a flagged mixture's pairwise distances shift
materially against the pure endpoints, that is a second result and it justifies the
module far better than the detection limit alone.

## 5. Scale

18 samples. On the 16 CPU / 64 GB configuration used for the ten-isolate run, this
is a few hours, and the mixtures are smaller than the originals.
