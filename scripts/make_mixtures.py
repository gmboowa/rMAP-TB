#!/usr/bin/env python3
"""
Build in-silico mixed-strain FASTQ sets for rMAP-TB mixed-infection validation.

Takes paired FASTQs for a set of MTBC isolates and writes read mixtures at
defined minor-strain fractions and total depths, plus a manifest and a Terra
data table. No third-party dependencies.

  python3 make_mixtures.py --fastq-dir RAW --out-dir MIX
  python3 make_mixtures.py --fastq-dir RAW --out-dir MIX --dry-run
"""
import argparse, gzip, io, os, random, sys, csv

GENOME = 4_411_532          # H37Rv NC_000962.3

# ----------------------------------------------------------------- design
# (label, isolate, lineage) pairs. "A" is the major strain at low minor pct.
PRIMARY   = (("ERR038255", "L4.8.2", "sensitive"),
             ("ERR13260175", "L2.2.1", "RR-TB"))
SECONDARY = (("ERR038266", "L3.1.1", "sensitive"),
             ("ERR13260106", "L4.1.2.1", "RR-TB"))
CONTROL   = (("ERR038254", "L4.8.2", "sensitive"),
             ("ERR038255", "L4.8.2", "sensitive"))

BASE_DEPTH   = 100
DEPTH_SERIES = [40, 70, 100, 140]


def design():
    """Yield dicts describing every mixture to build."""
    out = []
    a, b = PRIMARY

    # pure endpoints
    for iso, role in ((a, "A"), (b, "B")):
        out.append(dict(name=f"PURE_{iso[0]}_{BASE_DEPTH}x",
                        comps=[(iso[0], 1.0)], depth=BASE_DEPTH,
                        minor_pct=0, minor=None, series="pure",
                        expect="no mixed-infection signal"))

    # minor = B (resistant, different lineage) -- the main ladder
    for pct in (5, 10, 20, 30, 50):
        out.append(dict(name=f"MIX_{a[0]}_{b[0]}_{pct:02d}pctB_{BASE_DEPTH}x",
                        comps=[(a[0], 1 - pct / 100), (b[0], pct / 100)],
                        depth=BASE_DEPTH, minor_pct=pct, minor=b[0],
                        series="primary_minorB",
                        expect="possible mixed infection + heteroresistance"))

    # minor = A (susceptible), i.e. the reciprocal direction
    for pct in (5, 10, 20, 30):
        out.append(dict(name=f"MIX_{a[0]}_{b[0]}_{pct:02d}pctA_{BASE_DEPTH}x",
                        comps=[(a[0], pct / 100), (b[0], 1 - pct / 100)],
                        depth=BASE_DEPTH, minor_pct=pct, minor=a[0],
                        series="primary_minorA",
                        expect="possible mixed infection"))

    # depth series at a fixed 10 % minor, to separate fraction from depth
    for d in DEPTH_SERIES:
        if d == BASE_DEPTH:
            continue
        out.append(dict(name=f"MIX_{a[0]}_{b[0]}_10pctB_{d}x",
                        comps=[(a[0], 0.9), (b[0], 0.1)], depth=d,
                        minor_pct=10, minor=b[0], series="depth",
                        expect="possible mixed infection + heteroresistance"))

    # second, independent lineage pair
    c, e = SECONDARY
    for pct in (10, 20, 50):
        out.append(dict(name=f"MIX_{c[0]}_{e[0]}_{pct:02d}pctB_{BASE_DEPTH}x",
                        comps=[(c[0], 1 - pct / 100), (e[0], pct / 100)],
                        depth=BASE_DEPTH, minor_pct=pct, minor=e[0],
                        series="secondary",
                        expect="possible mixed infection + heteroresistance"))

    # specificity control: two near-identical same-sub-lineage isolates
    f, g = CONTROL
    out.append(dict(name=f"CTRL_SAMELIN_{f[0]}_{g[0]}_50pct_{BASE_DEPTH}x",
                    comps=[(f[0], 0.5), (g[0], 0.5)], depth=BASE_DEPTH,
                    minor_pct=50, minor=g[0], series="control",
                    expect="NO mixed-infection signal (specificity control)"))
    return out


# ----------------------------------------------------------------- fastq io
def op(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def fq_paths(d, acc):
    for p in ((f"{acc}_1.fastq.gz", f"{acc}_2.fastq.gz"),
              (f"{acc}_R1.fastq.gz", f"{acc}_R2.fastq.gz"),
              (f"{acc}_1.fq.gz", f"{acc}_2.fq.gz"),
              (f"{acc}_1.fastq", f"{acc}_2.fastq")):
        a, b = os.path.join(d, p[0]), os.path.join(d, p[1])
        if os.path.exists(a) and os.path.exists(b):
            return a, b
    raise SystemExit(f"ERROR: no paired FASTQ found for {acc} in {d}")


def survey(path):
    """Return (n_records, total_bases, modal_read_length) for one FASTQ."""
    n = bases = 0
    lens = {}
    with op(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                L = len(line.rstrip("\n"))
                n += 1
                bases += L
                lens[L] = lens.get(L, 0) + 1
    modal = max(lens, key=lens.get) if lens else 0
    return n, bases, modal


def take(src1, src2, keep, out1, out2):
    """Stream both mates, writing records whose index is in `keep`."""
    with op(src1) as a, op(src2) as b:
        idx = 0
        while True:
            r1 = [a.readline() for _ in range(4)]
            r2 = [b.readline() for _ in range(4)]
            if not r1[0] or not r2[0]:
                break
            if idx in keep:
                out1.write("".join(r1))
                out2.write("".join(r2))
            idx += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--seed", type=int, default=20261101)
    ap.add_argument("--genome-size", type=int, default=GENOME)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    plan = design()
    accs = sorted({c for m in plan for c, _ in m["comps"]})

    print("Surveying inputs ...")
    stats = {}
    for acc in accs:
        p1, p2 = fq_paths(a.fastq_dir, acc)
        n, b1, L1 = survey(p1)
        _, b2, L2 = survey(p2)
        depth = (b1 + b2) / a.genome_size
        stats[acc] = dict(p1=p1, p2=p2, n=n, depth=depth, rlen=max(L1, L2),
                          bp_per_pair=(b1 + b2) / max(n, 1))
        print(f"  {acc:<14} {n:>10,} pairs   {depth:7.1f}x available   "
              f"{L1}/{L2} bp reads")

    rl = {s["rlen"] for s in stats.values()}
    if len(rl) > 1:
        print(f"\n  NOTE: read lengths differ across isolates ({sorted(rl)} bp).")
        print("        Mixing is done on sequenced bases, not read counts, so the")
        print("        minor-strain allele fraction is still exact. Say so in the")
        print("        Methods, because a reviewer will check.")

    print("\nPlanned mixtures:")
    ok = True
    for m in plan:
        need = []
        for acc, frac in m["comps"]:
            want = m["depth"] * frac
            need.append(f"{acc}@{want:.0f}x")
            if want > stats[acc]["depth"]:
                ok = False
                m["skip"] = (f"needs {want:.0f}x of {acc} but only "
                             f"{stats[acc]['depth']:.0f}x available")
        flag = "  SKIP: " + m["skip"] if m.get("skip") else ""
        print(f"  {m['name']:<46} {' + '.join(need):<34}{flag}")
    if not ok:
        print("\nSome mixtures exceed available depth; they will be skipped.")
    if a.dry_run:
        return

    os.makedirs(a.out_dir, exist_ok=True)
    rng = random.Random(a.seed)
    rows = []
    for m in plan:
        if m.get("skip"):
            continue
        o1 = os.path.join(a.out_dir, m["name"] + "_1.fastq.gz")
        o2 = os.path.join(a.out_dir, m["name"] + "_2.fastq.gz")
        print(f"building {m['name']} ...", flush=True)
        with gzip.open(o1, "wt") as w1, gzip.open(o2, "wt") as w2:
            for acc, frac in m["comps"]:
                s = stats[acc]
                want_pairs = int(round(m["depth"] * frac * a.genome_size
                                       / s["bp_per_pair"]))
                want_pairs = min(want_pairs, s["n"])
                keep = set(rng.sample(range(s["n"]), want_pairs))
                take(s["p1"], s["p2"], keep, w1, w2)
        rows.append(dict(sample=m["name"], series=m["series"],
                         total_depth=m["depth"], minor_pct=m["minor_pct"],
                         minor_isolate=m["minor"] or "",
                         composition="; ".join(f"{c}={f:.2f}" for c, f in m["comps"]),
                         expected=m["expect"], read1=o1, read2=o2))

    with open(os.path.join(a.out_dir, "mixture_manifest.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(a.out_dir, "terra_tuberculosis_table.tsv"), "w", newline="") as fh:
        fh.write("entity:tuberculosis_id\tread1\tread2\n")
        for r in rows:
            fh.write(f"{r['sample']}\t{r['read1']}\t{r['read2']}\n")

    print(f"\n{len(rows)} mixtures written to {a.out_dir}")
    print("  mixture_manifest.tsv        expected result per sample")
    print("  terra_tuberculosis_table.tsv  upload after replacing paths with gs:// URIs")


if __name__ == "__main__":
    main()
