#!/usr/bin/env python3
"""
Download the public sequencing libraries used by the rMAP-TB validation panels
from the European Nucleotide Archive, and verify each file against the MD5 that
ENA publishes.

    # Panel 1 - the five parent isolates for the in-silico mixtures
    python3 fetch_reads.py --panel 1 --check        # report sizes, download nothing
    python3 fetch_reads.py --panel 1 --out-dir RAW

    # Panel 2 - the MTBC, NTM and contamination panel
    python3 fetch_reads.py --panel 2 --out-dir RAW_PANEL2

    # any accessions you like
    python3 fetch_reads.py ERR038255 ERR13260175 --out-dir RAW

    # or a TSV produced by find_validation_isolates.py
    python3 fetch_reads.py --from-table candidates.tsv --out-dir RAW_PANEL2

Downloads are resumable: re-running skips files that are already present and
verify. Only the Python standard library and curl are required.

Part of rMAP-TB (https://github.com/gmboowa/rMAP-TB), MIT licence.
"""
import argparse
import csv
import hashlib
import os
import subprocess
import sys
import time
import urllib.request

# --------------------------------------------------------------------- panels
# Panel 1: parents of the in-silico two-strain mixtures.
PANEL1 = [
    "ERR038254",    # L4.8.2, susceptible - same-strain control partner
    "ERR038255",    # L4.8.2, susceptible - primary parent
    "ERR038266",    # L3.1.1, susceptible - secondary parent
    "ERR13260106",  # L4.1.2.1, RR-TB     - secondary parent
    "ERR13260175",  # L2.2.1, RR-TB       - primary parent
]

# Panel 2: MTBC isolates spanning human- and animal-adapted lineages, NTM, and
# the M. avium library used to build the contamination mixtures.
PANEL2 = [
    "ERR1203067", "ERR3804031",                 # human-adapted L4, L6
    "ERR14187598", "ERR14187918", "SRR29440878",  # La1 (M. bovis, BCG)
    "SRR35607801", "SRR9157804",                # La2 (M. caprae), La3 (M. orygis)
    "ERR266109",                                # M. canettii
    "ERR970410", "SRR25473322", "SRR3500411", "SRR1239339",  # other animal MTBC
    "ERR12861768",                              # M. avium (also the contaminant)
    "SRR12599345", "SRR28906166",               # M. kansasii, M. intracellulare
    "ERR115016", "DRR360723", "SRR39696440",    # M. abscessus x2, M. fortuitum
]

API = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
       "?accession={acc}&result=read_run"
       "&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count"
       "&format=tsv&download=false")
UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/125.0 Safari/537.36")


def get(url, attempts=4, debug=False):
    """GET with a real User-Agent and backoff. Returns text, or raises."""
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": UA, "Accept": "text/plain"})
            with urllib.request.urlopen(req, timeout=90) as r:
                text = r.read().decode()
            if debug:
                print("    [debug] %s\n    [debug] %d bytes: %r"
                      % (url, len(text), text[:200]))
            return text
        except Exception as e:                      # noqa: BLE001
            last = e
            if i < attempts - 1:
                time.sleep(2 * (i + 1))
    raise last


def ena_path(acc, mate):
    """ENA mirrors FASTQ at a deterministic path, used when the portal API is
    unavailable. Validated against 740 known URLs: the top level is the first
    six characters, then for accessions with more than six digits a
    subdirectory of the surplus digits right-justified to three with zeros,
    then the accession directory."""
    pre, dig = acc[:3], acc[3:]
    sub = "" if len(dig) == 6 else dig[6 - len(dig):].rjust(3, "0") + "/"
    return ("ftp.sra.ebi.ac.uk/vol1/fastq/%s%s/%s%s/%s_%d.fastq.gz"
            % (pre, dig[:3], sub, acc, acc, mate))


def head_size(url):
    """Content-Length via HEAD, or None if the file is not there."""
    req = urllib.request.Request("https://" + url, method="HEAD",
                                 headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except Exception:                               # noqa: BLE001
        return None


def by_path(acc):
    """Fallback record built without the portal API. No MD5 is published on
    this route, so such downloads are verified by size only."""
    pair = {}
    for mate in (1, 2):
        u = ena_path(acc, mate)
        n = head_size(u)
        if n is None:
            return None
        pair[mate] = (u, None, n)
    return dict(acc=acc, pair=pair, reads="0", bases="0")


def query(accs, debug=False):
    """Resolve accessions to FASTQ URLs, sizes and MD5s.

    ENA is queried one accession at a time: the portal API is unreliable with
    comma-separated lists, and one bad accession fails the whole batch.
    """
    out, failed = [], []
    for acc in accs:
        url = API.format(acc=acc)
        try:
            text = get(url, debug=debug)
        except Exception as e:                      # noqa: BLE001
            failed.append((acc, "request failed: %s" % e, url))
            continue
        rows = [l.split("\t") for l in text.strip().splitlines()]
        if len(rows) < 2:
            failed.append((acc, "ENA returned a header but no data row", url))
            continue
        d = dict(zip(rows[0], rows[1]))
        ftp = [u for u in d.get("fastq_ftp", "").split(";") if u]
        if not ftp:
            failed.append((acc, "record exists but lists no FASTQ files", url))
            continue
        md5 = d["fastq_md5"].split(";")
        byt = d["fastq_bytes"].split(";")
        pair = {}
        for u, m, b in zip(ftp, md5, byt):
            base = u.rsplit("/", 1)[-1]
            if base.endswith("_1.fastq.gz"):
                pair[1] = (u, m, int(b))
            elif base.endswith("_2.fastq.gz"):
                pair[2] = (u, m, int(b))
        if 1 not in pair or 2 not in pair:
            names = [u.rsplit("/", 1)[-1] for u in ftp]
            failed.append((acc, "no _1/_2 pair; files are %s" % names, url))
            continue
        out.append(dict(acc=d["run_accession"], pair=pair,
                        reads=d.get("read_count", "0"),
                        bases=d.get("base_count", "0")))
        print("  resolved %s" % acc)

    if failed:
        print("\nThe ENA portal API did not answer for %d accession(s). "
              "Trying the direct FTP path instead." % len(failed))
        still = []
        for acc, why, url in failed:
            rec = by_path(acc)
            if rec:
                n = rec["pair"][1][2] + rec["pair"][2][2]
                print("  %s: found by path, %.2f GB (size-verified only, "
                      "no MD5)" % (acc, n / 1e9))
                out.append(rec)
            else:
                still.append((acc, why, url))
        if still:
            print("\nCould not resolve:")
            for acc, why, url in still:
                print("  %s: %s" % (acc, why))
                print("    %s" % url)
                print("    tried https://%s" % ena_path(acc, 1))
            print("\nPaste a URL into a browser to see what ENA returns.")
        if not out:
            sys.exit("\nNothing resolved; nothing downloaded.")
        print()
    return sorted(out, key=lambda r: r["acc"])


def md5_of(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def fetch(url, dest, want_md5, want_bytes):
    """Download one file, resuming if possible, and verify it."""
    if os.path.exists(dest) and os.path.getsize(dest) == want_bytes:
        if want_md5 is None:
            print("    present, size correct: %s" % os.path.basename(dest))
            return
        print("    present, verifying %s ..." % os.path.basename(dest),
              end="", flush=True)
        if md5_of(dest) == want_md5:
            print(" ok")
            return
        print(" MD5 MISMATCH, refetching")
        os.remove(dest)

    src = url if url.startswith("http") else "https://" + url
    cmd = ["curl", "-L", "-C", "-", "--retry", "5", "--retry-delay", "5",
           "-o", dest, src]
    print("    downloading %s (%.2f GB) ..."
          % (os.path.basename(dest), want_bytes / 1e9), flush=True)
    if subprocess.call(cmd) != 0:
        sys.exit("ERROR: download failed for %s. Re-run to resume." % url)

    if want_md5 is None:
        got = os.path.getsize(dest)
        if got != want_bytes:
            sys.exit("\nERROR: %s is %d bytes, expected %d."
                     % (dest, got, want_bytes))
        print("    size ok (no MD5 published for this route)")
        return
    print("    verifying ...", end="", flush=True)
    if md5_of(dest) != want_md5:
        sys.exit("\nERROR: MD5 mismatch for %s. Delete it and re-run." % dest)
    print(" ok")


def accessions_from_table(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or "run_accession" not in rows[0]:
        sys.exit("ERROR: %s has no run_accession column" % path)
    seen, accs = set(), []
    for r in rows:
        acc = (r["run_accession"] or "").strip()
        if acc and acc not in seen:
            seen.add(acc)
            accs.append(acc)
    return accs


def main():
    ap = argparse.ArgumentParser(
        description="Download rMAP-TB validation-panel reads from ENA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("accessions", nargs="*",
                    help="run accessions to fetch; overrides --panel")
    ap.add_argument("--panel", choices=["1", "2"],
                    help="1: the five mixture parents; "
                         "2: the MTBC/NTM/contamination panel")
    ap.add_argument("--from-table",
                    help="TSV with a run_accession column")
    ap.add_argument("--out-dir", default="RAW")
    ap.add_argument("--check", action="store_true",
                    help="resolve and report sizes, download nothing")
    ap.add_argument("--debug", action="store_true",
                    help="print each ENA request and response")
    a = ap.parse_args()

    if a.accessions:
        accs, label = a.accessions, "command line"
    elif a.from_table:
        accs, label = accessions_from_table(a.from_table), a.from_table
    elif a.panel == "1":
        accs, label = PANEL1, "Panel 1 (mixture parents)"
    elif a.panel == "2":
        accs, label = PANEL2, "Panel 2 (MTBC, NTM and contamination)"
    else:
        ap.error("give accessions, or --panel 1, --panel 2, or --from-table")

    print("%d accession(s) from %s\n" % (len(accs), label))
    recs = query(accs, debug=a.debug)

    total = 0
    print("\n%-14s %14s %14s   %s" % ("Accession", "Reads", "Bases", "Size"))
    for r in recs:
        n = r["pair"][1][2] + r["pair"][2][2]
        total += n
        print("%-14s %14s %14s   %.2f GB"
              % (r["acc"], r["reads"], r["bases"], n / 1e9))
    print("%-14s %44.2f GB" % ("TOTAL", total / 1e9))

    if a.check:
        print("\n--check: nothing downloaded. Re-run without it to fetch into %s/"
              % a.out_dir)
        return

    os.makedirs(a.out_dir, exist_ok=True)
    for r in recs:
        print("\n%s" % r["acc"])
        for mate in (1, 2):
            url, md5, nbytes = r["pair"][mate]
            dest = os.path.join(a.out_dir, "%s_%d.fastq.gz" % (r["acc"], mate))
            fetch(url, dest, md5, nbytes)

    print("\n%d isolate(s) verified in %s/" % (len(recs), a.out_dir))
    if a.panel == "1" or (not a.panel and not a.from_table):
        print("Next:  python3 make_mixtures.py --fastq-dir %s --out-dir MIX "
              "--dry-run" % a.out_dir)
    elif a.from_table:
        print("Next:  python3 build_ntm_arm.py --candidates %s --ntm-dir %s"
              % (a.from_table, a.out_dir))


if __name__ == "__main__":
    main()
