# TB-AMR-MTBC-Phylogenomics

A reproducible, containerized WDL/Cromwell workflow for *Mycobacterium tuberculosis* complex (MTBC) antimicrobial-resistance profiling and core-SNP phylogenomics.

The workflow supports paired-end Illumina FASTQ inputs and produces quality-control summaries, TB-Profiler drug-resistance and lineage reports, MTBC-only sample filtering, Snippy-based core-SNP outputs, optional Gubbins recombination filtering, IQ-TREE2 maximum-likelihood phylogeny, tree visualization, and interactive HTML reports for GitHub Pages.

## Recommended repository name

**`TB-AMR-MTBC-Phylogenomics`**

Why this name works:

▪ It clearly communicates tuberculosis, antimicrobial resistance, MTBC, and phylogenomics.  
▪ It matches the WDL workflow name: `TB_AMR_MTBC_Phylogenomics`.  
▪ It is suitable for a GitHub Pages URL such as:  

```text
https://gmboowa.github.io/TB-AMR-MTBC-Phylogenomics/
```

## Workflow overview

```text
Paired-end FASTQ files
⬇
Read trimming
⬇
FastQC + MultiQC
⬇
TB-Profiler AMR, species, and lineage profiling
⬇
MTBC-only sample filtering
⬇
Snippy per-sample variant calling
⬇
Snippy-core alignment
⬇
Optional Gubbins recombination filtering
⬇
IQ-TREE2 maximum-likelihood phylogeny
⬇
Tree visualization + integrated HTML report
```

## Key features

▪ Paired-end FASTQ input support  
▪ Adapter trimming before downstream analysis  
▪ FastQC and MultiQC quality-control reporting  
▪ TB-Profiler-based drug-resistance prediction  
▪ MTBC lineage and sub-lineage reporting  
▪ MTBC-only filtering before phylogenomic reconstruction  
▪ Snippy-based reference-guided variant calling  
▪ Core-SNP alignment generation using Snippy-core  
▪ Optional recombination filtering using Gubbins  
▪ Maximum-likelihood phylogeny using IQ-TREE2  
▪ Optional midpoint-rooted tree visualization  
▪ Interactive HTML reports suitable for GitHub Pages  
▪ Non-synonymous mutation reporting for key TB drug-resistance genes  

## Repository structure

```text
TB-AMR-MTBC-Phylogenomics/
├── README.md
├── TB_5.wdl
├── examples/
│   └── inputs.example.json
└── docs/
    ├── index.html
    ├── DEPLOYMENT.md
    ├── reports/
    │   ├── integrated_tb_amr_mtbc_phylogenomics_report.html
    │   └── tbprofiler_combined_report.html
    └── assets/
        └── .gitkeep
```

## Requirements

| Requirement | Purpose |
|---|---|
| Java | Required to run Cromwell |
| Cromwell | WDL workflow execution engine |
| Docker | Runs containerized bioinformatics tools |
| Paired-end FASTQ files | Input sequencing data |
| Adapter FASTA file | Required for trimming |
| MTBC GenBank reference | Required for Snippy and phylogenomic analysis |

## Main workflow inputs

| Input | Description |
|---|---|
| `input_reads` | Array of paired-end FASTQ files, ordered as R1 followed immediately by matching R2 |
| `adapters` | Adapter FASTA file used by the trimming step |
| `mtbc_reference_genbank` | MTBC reference genome in GenBank format |
| `do_trimming` | Enables read trimming |
| `do_quality_control` | Enables FastQC and MultiQC |
| `do_tb_profiler` | Enables TB-Profiler AMR and lineage analysis |
| `do_phylogeny` | Enables MTBC-only SNP phylogenomics |
| `use_gubbins` | Enables optional recombination filtering |
| `tbprofiler_docker` | Docker image for TB-Profiler |
| `snippy_reference_type` | Reference type; use `genbank` when providing a GenBank reference |
| `iqtree2_model` | IQ-TREE2 substitution model |
| `iqtree2_bootstraps` | Number of bootstrap replicates |
| `min_mtbc_samples_for_tree` | Minimum number of MTBC samples required to build a tree |
| `tb_drug_resistance_genes` | Comma-separated genes used for non-synonymous mutation reporting |

## Example input JSON

An example Cromwell input file is provided here:

```text
examples/inputs.example.json
```

The input FASTQ files must be ordered like this:

```json
{
  "TB_AMR_MTBC_Phylogenomics.input_reads": [
    "/path/to/sample1_1.fastq.gz",
    "/path/to/sample1_2.fastq.gz",
    "/path/to/sample2_1.fastq.gz",
    "/path/to/sample2_2.fastq.gz"
  ]
}
```

## Running the workflow

From the repository root:

```bash
java -jar cromwell-<version>.jar run TB_5.wdl --inputs examples/inputs.example.json
```

For example:

```bash
java -jar cromwell-92.jar run TB_5.wdl --inputs examples/inputs.example.json
```

## Recommended local Docker resources

For small to moderate MTBC datasets on a local workstation:

```text
CPUs:   8
Memory: 16 GB or higher
```

For larger datasets, especially when using Gubbins and IQ-TREE2, consider increasing memory and CPU allocation where possible.

## Main outputs

### Quality control

▪ Trimmed FASTQ files  
▪ FastQC reports  
▪ MultiQC report  
▪ Trimming summary table  

### TB-Profiler and MTBC filtering

▪ TB-Profiler JSON outputs  
▪ Combined TB-Profiler HTML report  
▪ Resistance profile summary  
▪ MTBC-positive sample list  
▪ Non-MTBC or low-confidence excluded sample list  

### SNP phylogenomics

▪ Per-sample Snippy variant calls  
▪ Snippy-core alignment  
▪ Core SNP VCF  
▪ Optional Gubbins recombination-filtered alignment  
▪ IQ-TREE2 maximum-likelihood tree  
▪ Bootstrap-supported tree file  
▪ Tree visualization as PNG or SVG  

### Interactive reports

▪ `integrated_tb_amr_mtbc_phylogenomics_report.html`  
▪ `tbprofiler_combined_report.html`  

## GitHub Pages report site

This repository is prepared for GitHub Pages using the `docs/` folder.

After deployment, the reports will be available at:

```text
https://gmboowa.github.io/TB-AMR-MTBC-Phylogenomics/
```

The `docs/index.html` file provides a card-style dashboard with links to the main reports.

## How to publish the reports on GitHub Pages

1. Push the repository to GitHub.
2. Go to the repository page on GitHub.
3. Open **Settings**.
4. Open **Pages**.
5. Under **Build and deployment**, set:

```text
Source: Deploy from a branch
Branch: main
Folder: /docs
```

6. Save the settings.
7. Wait for GitHub Pages to build the site.
8. Open:

```text
https://gmboowa.github.io/TB-AMR-MTBC-Phylogenomics/
```

## Updating reports after each workflow run

After a successful workflow run, copy the new HTML reports and assets into the `docs/` folder:

```bash
mkdir -p docs/reports docs/assets

cp /path/to/integrated_tb_amr_mtbc_phylogenomics_report.html docs/reports/
cp /path/to/tbprofiler_combined_report.html docs/reports/

# If your reports use external assets, copy them too:
cp -R /path/to/report_assets/* docs/assets/
```

Then commit and push:

```bash
git add README.md TB_5.wdl examples/ docs/
git commit -m "Add TB AMR MTBC phylogenomics workflow and reports site"
git push origin main
```

## Suggested first GitHub commands

```bash
mkdir TB-AMR-MTBC-Phylogenomics
cd TB-AMR-MTBC-Phylogenomics

git init

git add .
git commit -m "Initial TB AMR MTBC phylogenomics workflow repository"

git branch -M main
git remote add origin https://github.com/gmboowa/TB-AMR-MTBC-Phylogenomics.git
git push -u origin main
```

## Interpretation guidance

The integrated report should be interpreted using multiple layers of evidence:

▪ TB-Profiler resistance profile  
▪ MTBC lineage and sub-lineage  
▪ Core-SNP phylogenetic clustering  
▪ Bootstrap support  
▪ Recombination-filtered alignment if Gubbins is enabled  
▪ Metadata such as country, year, collection site, and sample source  

Close clustering alone should not be treated as proof of transmission without epidemiological and sampling context.

## Suggested citation / acknowledgement

If you use this workflow, please cite or acknowledge:

```text
TB-AMR-MTBC-Phylogenomics: a WDL/Cromwell workflow for MTBC antimicrobial-resistance profiling and core-SNP phylogenomics.
```

## License

Add an appropriate license before public release. Recommended options:

▪ MIT License for permissive open-source reuse  
▪ Apache-2.0 License for permissive reuse with explicit patent language  
▪ GPL-3.0 License if derivative workflows should remain open source  

## Maintainer

**Gerald Mboowa**  
Senior Data Engineer, Broad Institute  
Bioinformatics and genomic epidemiology implementation scientist
