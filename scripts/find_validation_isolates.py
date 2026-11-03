#!/usr/bin/env python3
"""
Find candidate ENA runs for the rMAP-TB NTM / non-tuberculosis MTBC arm.

Resolves each species to a taxon id through ENA's taxonomy service (so no
taxon id is hard-coded or guessed), queries the portal for paired Illumina
WGS runs, ranks them by how close they sit to a target depth, and writes a
candidate table for review.

  python3 find_validation_isolates.py                 # all groups
  python3 find_validation_isolates.py --group ntm_in_db
  python3 find_validation_isolates.py --per-species 5

Nothing is downloaded. Standard library only.
"""
import argparse, json, re, sys, time, urllib.parse, urllib.request

GENOME = 4_411_532          # H37Rv, used only to express depth comparably
UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/125.0 Safari/537.36")

TAXONOMY = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{name}"
ANYNAME = "https://www.ebi.ac.uk/ena/taxonomy/rest/any-name/{name}"
SUGGEST = "https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-search/{name}"

# Several of these species have been renamed. ENA's taxonomy has no exact
# scientific-name match for the older binomial, so each entry carries the
# names actually in use. The 2018 split of Mycobacterium moved the rapid
# growers into Mycobacteroides and Mycolicibacterium, and the animal-adapted
# MTBC members are now recorded as variants of M. tuberculosis.
ALIASES = {
    "Mycobacterium bovis": [
        "Mycobacterium tuberculosis variant bovis",
        "Mycobacterium bovis",
    ],
    "Mycobacterium africanum": [
        "Mycobacterium tuberculosis variant africanum",
        "Mycobacterium africanum",
    ],
    "Mycobacterium canettii": [
        "Mycobacterium canettii",
        "Mycobacterium tuberculosis variant canettii",
    ],
    "Mycobacterium abscessus": [
        "Mycobacteroides abscessus",
        "Mycobacterium abscessus",
    ],
    "Mycobacterium fortuitum": [
        "Mycolicibacterium fortuitum",
        "Mycobacteroides fortuitum",
        "Mycobacterium fortuitum",
    ],
    "Mycobacterium bovis BCG": [
        "Mycobacterium tuberculosis variant bovis BCG",
        "Mycobacterium bovis BCG",
    ],
    "Mycobacterium caprae": [
        "Mycobacterium tuberculosis variant caprae",
        "Mycobacterium caprae",
    ],
    "Mycobacterium microti": [
        "Mycobacterium tuberculosis variant microti",
        "Mycobacterium microti",
    ],
    "Mycobacterium orygis": [
        "Mycobacterium orygis",
        "Mycobacterium tuberculosis variant orygis",
    ],
    "Mycobacterium pinnipedii": [
        "Mycobacterium tuberculosis variant pinnipedii",
        "Mycobacterium pinnipedii",
    ],
    "Mycobacterium mungi": [
        "Mycobacterium tuberculosis variant mungi",
        "Mycobacterium mungi",
    ],
    "Mycobacterium suricattae": [
        "Mycobacterium tuberculosis variant suricattae",
        "Mycobacterium suricattae",
    ],
}
SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"

# Each entry: species, group, why it is in the panel, expected rMAP-TB result.
PANEL = [
    # --- MTBC members that are not M. tuberculosis -----------------------
    ("Mycobacterium bovis", "mtbc_other", 2,
     "zoonotic MTBC; absent from the Kraken2 database",
     "routed MTBC; TB-Profiler La1 / M. bovis; Kraken2 will say M. tuberculosis"),
    ("Mycobacterium africanum", "mtbc_other", 2,
     "MTBC lineages 5 and 6; absent from the Kraken2 database",
     "routed MTBC; TB-Profiler L5 or L6; Kraken2 will say M. tuberculosis"),
    ("Mycobacterium canettii", "mtbc_other", 1,
     "MTBC and IS in the Kraken2 database",
     "routed MTBC; Kraken2 should name M. canettii correctly"),
    ("Mycobacterium bovis BCG", "mtbc_other", 1,
     "vaccine strain; exercises the BCG label in the lineage fix",
     "routed MTBC; TB-Profiler La1 or BCG; must not be called a human lineage"),
    ("Mycobacterium caprae", "mtbc_other", 1,
     "animal-adapted MTBC; absent from the Kraken2 database",
     "routed MTBC; TB-Profiler La2 or M. caprae"),
    ("Mycobacterium orygis", "mtbc_other", 1,
     "animal-adapted MTBC, increasingly reported in humans",
     "routed MTBC; TB-Profiler La3 or M. orygis"),
    ("Mycobacterium microti", "mtbc_other", 1,
     "animal-adapted MTBC; rare in humans",
     "routed MTBC; TB-Profiler M. microti"),
    ("Mycobacterium pinnipedii", "mtbc_other", 1,
     "animal-adapted MTBC; few public genomes, may be skipped",
     "routed MTBC; TB-Profiler M. pinnipedii"),
    ("Mycobacterium mungi", "mtbc_other", 1,
     "animal-adapted MTBC; very few public genomes, likely skipped",
     "routed MTBC if sequenced data exist"),
    ("Mycobacterium suricattae", "mtbc_other", 1,
     "animal-adapted MTBC; very few public genomes, likely skipped",
     "routed MTBC if sequenced data exist"),

    # --- NTM present in the Kraken2 database ------------------------------
    ("Mycobacterium avium", "ntm_in_db", 1,
     "MAC; in the database",
     "routed NTM; Kraken2 should name M. avium; TB-Profiler skipped"),
    ("Mycobacterium intracellulare", "ntm_in_db", 1,
     "MAC; in the database",
     "routed NTM; Kraken2 should name M. intracellulare"),
    ("Mycobacterium kansasii", "ntm_in_db", 1,
     "slow grower; in the database",
     "routed NTM; Kraken2 should name M. kansasii"),

    # --- NTM absent from the Kraken2 database -----------------------------
    ("Mycobacterium abscessus", "ntm_absent", 2,
     "commonest clinical rapid grower; NOT in the database",
     "unknown: either no species-level call, or a confident WRONG species. "
     "This is the point of including it."),
    ("Mycobacterium fortuitum", "ntm_absent", 1,
     "rapid grower; NOT in the database",
     "as above"),
]


def get(url, attempts=4):
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA,
                                                       "Accept": "*/*"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode()
        except Exception as e:
            last = e
            if i < attempts - 1:
                time.sleep(2 * (i + 1))
    raise last


def _lookup(url, name):
    try:
        data = json.loads(get(url.format(name=urllib.parse.quote(name)), attempts=2))
    except Exception:
        return None, None
    if isinstance(data, dict):
        data = [data]
    if not data:
        return None, None
    for rec in data:
        if rec.get("scientificName", "").lower() == name.lower() and rec.get("taxId"):
            return rec["taxId"], rec["scientificName"]
    rec = data[0]
    if rec.get("taxId"):
        return rec["taxId"], rec.get("scientificName", name)
    return None, None


def taxid_for(name):
    """Resolve a species to a taxon id, trying the names actually in use.

    Nothing is hard-coded: the ids come from ENA. Only the NAMES are listed,
    because several of these species have been reclassified and the older
    binomial no longer matches ENA's scientific-name index.
    """
    tried = []
    for candidate in ALIASES.get(name, [name]):
        tried.append(candidate)
        for url in (TAXONOMY, ANYNAME):
            tid, sci = _lookup(url, candidate)
            if tid:
                return tid, None if sci.lower() == name.lower() else \
                    "resolved as %r" % sci
    return None, "no taxon id for any of: %s" % ", ".join(tried)


def runs_for(taxid, subtaxa, limit):
    q = ('tax_tree({t})' if subtaxa else 'tax_eq({t})').format(t=taxid)
    q += (' AND instrument_platform="ILLUMINA"'
          ' AND library_layout="PAIRED"'
          ' AND library_strategy="WGS"')
    params = urllib.parse.urlencode({
        "result": "read_run",
        "query": q,
        "fields": ("run_accession,sample_accession,scientific_name,"
                   "instrument_model,read_count,base_count,fastq_ftp,"
                   "fastq_bytes,fastq_md5"),
        "format": "tsv",
        "limit": str(limit),
    })
    txt = get(SEARCH + "?" + params)
    rows = [l.split("\t") for l in txt.strip().splitlines()]
    if len(rows) < 2:
        return []
    hdr = rows[0]
    out = []
    for r in rows[1:]:
        d = dict(zip(hdr, r))
        ftp = [u for u in d.get("fastq_ftp", "").split(";") if u]
        names = [u.rsplit("/", 1)[-1] for u in ftp]
        if not any(n.endswith("_1.fastq.gz") for n in names):
            continue
        if not any(n.endswith("_2.fastq.gz") for n in names):
            continue
        try:
            bases = int(d.get("base_count") or 0)
        except ValueError:
            continue
        if bases <= 0:
            continue
        d["_depth"] = bases / GENOME
        try:
            d["_gb"] = sum(int(b) for b in d["fastq_bytes"].split(";") if b) / 1e9
        except ValueError:
            d["_gb"] = 0.0
        out.append(d)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-species", type=int, default=3,
                    help="candidates to show per species (default 3)")
    ap.add_argument("--group", choices=["mtbc_other", "ntm_in_db", "ntm_absent"],
                    help="restrict to one group")
    ap.add_argument("--target-depth", type=float, default=150.0)
    ap.add_argument("--min-depth", type=float, default=60.0)
    ap.add_argument("--max-gb", type=float, default=3.0,
                    help="skip runs larger than this per run (default 3 GB)")
    ap.add_argument("--subtaxa", action="store_true",
                    help="include sub-taxa (e.g. BCG under M. bovis)")
    ap.add_argument("--out", default="ntm_mtbc_candidates.tsv")
    a = ap.parse_args()

    panel = [p for p in PANEL if not a.group or p[1] == a.group]
    chosen, notes = [], []

    for species, group, want, why, expect in panel:
        print(f"\n=== {species}  [{group}]  need {want}")
        print(f"    {why}")
        taxid, note = taxid_for(species)
        if not taxid:
            print(f"    SKIP: {note}")
            try:
                sug = json.loads(get(SUGGEST.format(
                    name=urllib.parse.quote(species.split()[-1])), attempts=2))
                names = [r.get("displayName") or r.get("scientificName")
                         for r in (sug if isinstance(sug, list) else [sug])][:6]
                if names:
                    print("    ENA suggests: " + "; ".join(n for n in names if n))
            except Exception:
                pass
            notes.append((species, note))
            continue
        print(f"    taxid {taxid}" + (f"  ({note})" if note else ""))
        try:
            runs = runs_for(taxid, a.subtaxa, 500)
        except Exception as e:
            print(f"    SKIP: portal query failed: {e}")
            notes.append((species, str(e)))
            continue
        runs = [r for r in runs
                if r["_depth"] >= a.min_depth and r["_gb"] <= a.max_gb]
        if not runs:
            print("    SKIP: no run met the depth and size limits")
            notes.append((species, "no run met depth/size limits"))
            continue
        runs.sort(key=lambda r: abs(r["_depth"] - a.target_depth))
        shown = runs[:max(a.per_species, want)]
        for i, r in enumerate(shown):
            mark = "*" if i < want else " "
            print(f"   {mark} {r['run_accession']:<14} {r['_depth']:6.0f}x  "
                  f"{r['_gb']:5.2f} GB  {r['scientific_name']:<34}"
                  f"{r.get('instrument_model','')}")
            if i < want:
                chosen.append(dict(species=species, group=group,
                                   run_accession=r["run_accession"],
                                   sample_accession=r.get("sample_accession", ""),
                                   depth="%.0f" % r["_depth"],
                                   gb="%.2f" % r["_gb"],
                                   model=r.get("instrument_model", ""),
                                   why=why, expected=expect))

    if not chosen:
        sys.exit("\nNothing selected. Nothing to download.")

    cols = ["species", "group", "run_accession", "sample_accession",
            "depth", "gb", "model", "why", "expected"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in chosen:
            fh.write("\t".join(str(c[k]) for k in cols) + "\n")

    total = sum(float(c["gb"]) for c in chosen)
    print(f"\n{len(chosen)} runs selected, {total:.1f} GB total -> {a.out}")
    print("Rows marked * were auto-selected. Edit the file to swap any of them,")
    print("then:  python3 fetch_reads.py --from-table %s --out-dir RAW_NTM" % a.out)
    if notes:
        print("\nSpecies with no candidate:")
        for s, e in notes:
            print(f"  {s}: {e}")
        print("If the reason is a name rather than depth, add the name ENA "
              "suggests to ALIASES at the top of this script and re-run.")
        print("Otherwise try --subtaxa, a lower --min-depth, or a larger --max-gb.")


if __name__ == "__main__":
    main()
