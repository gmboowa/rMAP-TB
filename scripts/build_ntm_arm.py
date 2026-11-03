#!/usr/bin/env python3
"""
Assemble the rMAP-TB NTM / non-tuberculosis MTBC validation arm.

Takes the pure isolates fetched into RAW_NTM, optionally builds two
MTBC + NTM contamination mixtures from an M. tuberculosis isolate already in
RAW, and writes a manifest stating the expected result for every sample plus
a Terra data table.

  python3 build_ntm_arm.py --candidates ntm_mtbc_candidates.tsv
  python3 build_ntm_arm.py --candidates ... --no-contamination

Standard library only.
"""
import argparse, csv, gzip, os, random, sys

GENOME = 4_411_532


def op(p, mode="rt"):
    return gzip.open(p, mode) if str(p).endswith(".gz") else open(p, mode)


def fq_paths(d, acc):
    for a, b in ((f"{acc}_1.fastq.gz", f"{acc}_2.fastq.gz"),
                 (f"{acc}_R1.fastq.gz", f"{acc}_R2.fastq.gz"),
                 (f"{acc}_1.fastq", f"{acc}_2.fastq")):
        p1, p2 = os.path.join(d, a), os.path.join(d, b)
        if os.path.exists(p1) and os.path.exists(p2):
            return p1, p2
    return None, None


def survey(path):
    n = bases = 0
    with op(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                n += 1
                bases += len(line.rstrip("\n"))
    return n, bases


def take(s1, s2, keep, o1, o2):
    with op(s1) as a, op(s2) as b:
        i = 0
        while True:
            r1 = [a.readline() for _ in range(4)]
            r2 = [b.readline() for _ in range(4)]
            if not r1[0] or not r2[0]:
                break
            if i in keep:
                o1.write("".join(r1))
                o2.write("".join(r2))
            i += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="ntm_mtbc_candidates.tsv")
    ap.add_argument("--ntm-dir", default="RAW_NTM")
    ap.add_argument("--mtbc-dir", default="RAW",
                    help="where the M. tuberculosis isolate lives")
    ap.add_argument("--mtbc-accession", default="ERR038255")
    ap.add_argument("--contaminant",
                    help="accession to mix with the M. tuberculosis isolate; "
                         "default is the first ntm_in_db run in the table")
    ap.add_argument("--out-dir", default="MIX_NTM")
    ap.add_argument("--depth", type=float, default=100.0)
    ap.add_argument("--no-contamination", action="store_true")
    ap.add_argument("--seed", type=int, default=20261101)
    a = ap.parse_args()

    with open(a.candidates) as fh:
        cands = list(csv.DictReader(fh, delimiter="\t"))
    if not cands:
        sys.exit("ERROR: %s is empty" % a.candidates)

    rows = []
    missing = []
    for c in cands:
        acc = c["run_accession"]
        p1, p2 = fq_paths(a.ntm_dir, acc)
        if not p1:
            missing.append(acc)
            continue
        rows.append(dict(sample=f"{c['group'].upper()}_{acc}",
                         group=c["group"], composition=f"{acc}=1.00 (pure)",
                         species_submitted=c["species"],
                         expected=c["expected"], read1=p1, read2=p2))
    if missing:
        print("Not found in %s, skipped: %s" % (a.ntm_dir, ", ".join(missing)))
    if not rows:
        sys.exit("ERROR: none of the candidate FASTQs are in %s" % a.ntm_dir)

    if not a.no_contamination:
        cont = a.contaminant
        if not cont:
            inn = [c for c in cands if c["group"] == "ntm_in_db"]
            cont = inn[0]["run_accession"] if inn else None
        m1, m2 = fq_paths(a.mtbc_dir, a.mtbc_accession)
        c1, c2 = fq_paths(a.ntm_dir, cont) if cont else (None, None)
        if not (m1 and c1):
            print("Skipping contamination mixtures: need %s in %s and %s in %s"
                  % (a.mtbc_accession, a.mtbc_dir, cont, a.ntm_dir))
        else:
            os.makedirs(a.out_dir, exist_ok=True)
            rng = random.Random(a.seed)
            stats = {}
            for acc, (x1, x2) in ((a.mtbc_accession, (m1, m2)), (cont, (c1, c2))):
                n, b1 = survey(x1)
                _, b2 = survey(x2)
                stats[acc] = dict(p1=x1, p2=x2, n=n,
                                  depth=(b1 + b2) / GENOME,
                                  bpp=(b1 + b2) / max(n, 1))
                print(f"  {acc:<14} {stats[acc]['depth']:7.1f}x available")
            for pct in (10, 50):
                name = (f"CONTAM_{a.mtbc_accession}_{cont}_{pct:02d}pctNTM"
                        f"_{int(a.depth)}x")
                o1 = os.path.join(a.out_dir, name + "_1.fastq.gz")
                o2 = os.path.join(a.out_dir, name + "_2.fastq.gz")
                print("building", name, "...", flush=True)
                with gzip.open(o1, "wt") as w1, gzip.open(o2, "wt") as w2:
                    for acc, frac in ((a.mtbc_accession, 1 - pct / 100),
                                      (cont, pct / 100)):
                        s = stats[acc]
                        want = int(round(a.depth * frac * GENOME / s["bpp"]))
                        want = min(want, s["n"])
                        take(s["p1"], s["p2"],
                             set(rng.sample(range(s["n"]), want)), w1, w2)
                rows.append(dict(
                    sample=name, group="contamination",
                    composition=f"{a.mtbc_accession}={1-pct/100:.2f}; {cont}={pct/100:.2f}",
                    species_submitted="M. tuberculosis + NTM",
                    expected=("MTBC routing must survive %d%% NTM contamination; "
                              "species call and lineage must stay M. tuberculosis"
                              % pct),
                    read1=o1, read2=o2))

    cols = ["sample", "group", "composition", "species_submitted",
            "expected", "read1", "read2"]
    with open("ntm_arm_manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open("terra_ntm_table.tsv", "w") as fh:
        fh.write("entity:tuberculosis_id\tread1\tread2\n")
        for r in rows:
            fh.write(f"{r['sample']}\t{r['read1']}\t{r['read2']}\n")

    print(f"\n{len(rows)} samples in the NTM arm")
    for g in sorted({r["group"] for r in rows}):
        print(f"  {g:<16} {sum(1 for r in rows if r['group']==g)}")
    print("\n  ntm_arm_manifest.tsv   expected result per sample")
    print("  terra_ntm_table.tsv    upload after replacing paths with gs:// URIs")


if __name__ == "__main__":
    main()
