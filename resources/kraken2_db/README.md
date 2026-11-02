# Kraken2 species-typing database

This document describes the Kraken2 database used by the `SPECIES_TYPING` task of rMAP-TB, how it is distributed, what it can and cannot resolve, and how to verify the exact copy used in a given run.

## 1. Location and distribution

The database is not stored in this repository. It is embedded in the species-typing container image and is read by the workflow at a fixed path:

| Item | Value |
|---|---|
| Container image | `gmboowa/mycobacterium-kraken2-bracken:2026.05` (workflow input `species_typing_docker`) |
| Database path inside the image | `/opt/kraken2_db/mycobacterium` |
| Classifier | Kraken2, paired-end mode (`--paired`), default k-mer and confidence settings |
| Bracken | A Bracken k-mer distribution is present in the image but **Bracken is not executed** by rMAP-TB; abundance re-estimation is not performed |

Because the image tag is pinned in the workflow, every rMAP-TB release classifies against a fixed database build. Changing `species_typing_docker` changes the database and should be reported with the results.

## 2. Composition

The database is a curated, *Mycobacterium*-only reference set built from RefSeq assemblies:

| Property | Value |
|---|---|
| Reference sequences | 76 |
| Species represented | 50 |
| Genera | *Mycobacterium* only (including taxa recently reassigned to *Mycobacteroides*, *Mycolicibacterium*, *Mycolicibacillus* and *Mycolicibacter*, which Kraken2 reports under their NCBI taxonomy names) |
| MTBC genomes | 2: *M. tuberculosis* H37Rv (NC_000962.3) and *M. canettii* CIPT 140010059 (NC_015848.1) |
| MTBC complex node | NCBI taxid 77643, reported by Kraken2 at rank `G1` |
| Human sequence | None |
| Other bacterial, viral, fungal or plasmid sequence | None |

The full list of accessions and their taxonomic assignments is written by every workflow run (Section 5) and should be consulted for the authoritative content of the build in use.

## 3. Design rationale

A genus-restricted database was chosen deliberately:

- **Speed and footprint.** The database is small enough to ship inside the container, so no database download or external storage is required on Terra or on local systems.
- **Resolution within the genus.** With the search space limited to mycobacteria, reads are not drawn toward distant taxa, and the top species-level call within the genus is a direct, interpretable statement of the most probable mycobacterium present.
- **Routing, not metagenomics.** The task's purpose is to decide whether a sample is supported as MTBC and, if not, to name the most probable non-tuberculous mycobacterium. It is not intended to profile the complete microbial content of a specimen.

## 4. Limitations and interpretation

These limitations follow directly from the composition above and are stated in the integrated report wherever species calls are shown.

1. **Non-mycobacterial reads are unclassified.** Reads from the host, from other bacteria or from contaminating organisms outside the genus are returned as unclassified. A low percentage of classified reads therefore indicates a low mycobacterial content, not a failed library, and should be read together with the absolute read counts.
2. **Absent taxa cannot be named.** A mycobacterium that is not among the 50 species will be assigned to its nearest represented relative or to a higher node. The percentage of classified reads supporting the call, the per-sample `min_species_call_percent` threshold (default 5 %) and the "low confidence / unresolved" handling in the router exist to limit the impact of this.
3. **No species discrimination within the MTBC.** Only *M. tuberculosis* and *M. canettii* are present. *M. bovis*, *M. bovis* BCG, *M. africanum*, *M. microti*, *M. orygis*, *M. caprae*, *M. pinnipedii* and *M. mungi* are not, so any of these will be called *M. tuberculosis*. A call of *M. tuberculosis* is evidence of MTBC membership; lineage and species within the complex are taken from TB-Profiler downstream.
4. **Complex-node support is used as evidence.** The router reads the read count and percentage assigned at or below taxid 77643 and applies two rules: a majority (`mtbc_node_override_percent`, default 50 %) with at least 1,000 reads routes a sample as MTBC even when a non-MTBC species is the single highest species call, and 10 % with at least 1,000 reads rescues an otherwise unresolved call. In the validation panel, genuine NTM showed 0.01–0.39 % complex-node background, so a majority is not attributable to background k-mer sharing.
5. **Closely related NTM may be interchanged.** Members of species complexes with high genomic similarity (for example the *M. avium* complex or the *M. abscessus* subspecies) can share k-mers; the reported species is the most probable call and should be confirmed by a dedicated method when the distinction matters clinically.
6. **Not a diagnostic.** Species calls are computational predictions produced for research and surveillance use.

## 5. Verifying the database used in a run

The `SPECIES_TYPING` task records the identity of the database it classified against, so the exact build is traceable without access to the image. The following files are produced in the task's `species_typing/` output directory:

| File | Content |
|---|---|
| `kraken2_db_manifest.tsv` | Database path, size and MD5 checksum of `hash.k2d`, MD5 of `taxo.k2d`, number of reference sequences, and whether a Bracken distribution is present |
| `kraken2_db_inspect.txt` | Full output of `kraken2-inspect --db`, listing every taxon in the database with its k-mer counts |
| `kraken2_db_species.tsv` | Species-level (`S`) rows extracted from the inspect output |
| `kraken2_db_accessions.tsv` | Accession-to-taxid map (`seqid2taxid.map`) of the 76 reference sequences |

To reproduce these outside the workflow:

```bash
docker run --rm gmboowa/mycobacterium-kraken2-bracken:2026.05 \
  kraken2-inspect --db /opt/kraken2_db/mycobacterium | awk '$4=="S"' | sort -k2,2nr

docker run --rm gmboowa/mycobacterium-kraken2-bracken:2026.05 \
  cat /opt/kraken2_db/mycobacterium/seqid2taxid.map
```

When publishing results, cite the image tag together with the `hash.k2d` MD5 from `kraken2_db_manifest.tsv`; the pair identifies the database build unambiguously.

## 6. Using a different database

A different Kraken2 database can be used by building an image that places it at `/opt/kraken2_db/mycobacterium` and passing that image as `species_typing_docker`. The routing logic assumes that:

- the database contains the MTBC complex node (taxid 77643) so that complex-level support can be read from the report, and
- at least one MTBC genome is present at species level.

A database that includes host or non-mycobacterial sequence will change the meaning of "percentage of classified reads" and the behaviour of `min_species_call_percent`; thresholds should be re-evaluated before such a database is used for routing.

## 7. Change log

| Image tag | Change |
|---|---|
| `2026.05` | Current build: 76 RefSeq sequences, 50 species, 2 MTBC genomes. Bracken distribution included but unused. |
