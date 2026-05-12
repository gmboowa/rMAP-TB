## rMAP-TB

**rMAP-TB** is a reproducible, Dockerized WDL/Cromwell workflow for public-health-oriented analysis of *Mycobacterium tuberculosis* complex (MTBC) genomic data. It supports paired-end Illumina FASTQ inputs & integrates species typing, TB drug-resistance profiling, lineage interpretation, MTBC-only sample filtering, core-SNP phylogenomics, and interactive surveillance reporting.

The workflow performs read trimming, sequence quality control, Kraken2/Bracken-based Mycobacteria species typing, TB-Profiler resistance and lineage profiling, Snippy-based variant calling, Snippy-core alignment, drug-resistance-associated non-synonymous mutation summarization, pairwise SNP distance estimation, SNP cluster interpretation, optional Gubbins recombination filtering, IQ-TREE2 maximum-likelihood phylogeny & ETE3 tree visualization.

rMAP-TB generates integrated HTML reports & downloadable public-health surveillance outputs, including QC filtering rationale, TB-Profiler mutation-level resistance evidence, lineage distribution summaries, SNP distance heatmaps, pairwise SNP distance tables, SNP cluster summaries & surveillance metadata TSV files.


## Workflow overview
<p align="center">
  <img src="docs/assets/workflow/rMAP_TB.png"
       alt="rMAP-TB workflow"
       width="100%">
</p>

```text
Paired-end FASTQ files
        ⬇
Read trimming with Trimmomatic
        ⬇
Sequence quality control with FastQC
        ⬇
QC aggregation with MultiQC
        ⬇
Mycobacteria species typing with Kraken2 + Bracken
        ⬇
MTBC / non-MTBC Mycobacteria routing
        ⬇

┌──────────────────────────────────────────────────────────────┐
│                 Non-MTBC Mycobacteria / NTM branch           │
│                              ⬇                               │
│                    NTM speciation summary                    │
│                              ⬇                               │
│              Most probable NTM species identified            │
│                              ⬇                               │
│              Species-level evidence and MTBC support         │
│                              ⬇                               │
│              Exclusion from MTBC-specific analysis           │
│                              ⬇                               │
│              Non-MTBC Mycobacteria species summary           │
│                              ⬇                               │
│                    Integrated HTML report                    │
└──────────────────────────────────────────────────────────────┘

        ⬇
MTBC-supported samples only
        ⬇
TB-Profiler species, lineage & AMR profiling
        ⬇
MTBC-only sample filtering
        ⬇
Snippy per-sample variant calling
        ⬇
Mean-depth extraction & variant summary generation
        ⬇
Snippy-core core-genome alignment
        ⬇
Drug-resistance-associated non-synonymous mutation summary
        ⬇
Pairwise SNP distance estimation
        ⬇
SNP cluster interpretation
        ⬇
Lineage distribution summary
        ⬇
SNP distance heatmap generation
        ⬇
QC filtering rationale & surveillance metadata export
        ⬇
Optional Gubbins recombination filtering
        ⬇
IQ-TREE2 maximum-likelihood phylogeny
        ⬇
ETE3 phylogenetic tree visualization
        ⬇
Integrated HTML report with downloadable surveillance outputs
```
## Key features

- **⬤ Paired-end Illumina FASTQ input support**  
- **⬤ Adapter trimming & read preprocessing using Trimmomatic**  
- **⬤ Sequence quality-control assessment using FastQC**  
- **⬤ Aggregated QC reporting using MultiQC**  
- **⬤ Mycobacteria species typing using Kraken2 + Bracken**  
- **⬤ TB-Profiler-based MTBC species, lineage, sub-lineage & drug-resistance profiling**  
- **⬤ WHO-aligned TB drug-resistance classification, including HR-TB, RR-TB, MDR-TB, Pre-XDR-TB & XDR-TB**  
- **⬤ Mutation-level TB-Profiler resistance evidence reporting, including drug, gene, mutation, confidence & evidence fields** 
- **⬤ MTBC-only sample filtering before downstream phylogenomics**  
- **⬤ Snippy-based reference-guided per-sample variant calling**  
- **⬤ Mean-depth extraction & per-sample variant summary reporting**  
- **⬤ Core-genome SNP alignment generation using Snippy-core**  
- **⬤ Non-synonymous mutation reporting for key TB drug-resistance-associated genes**  
- **⬤ Pairwise SNP distance estimation from MTBC core-genome alignments**  
- **⬤ SNP cluster interpretation using configurable SNP-distance thresholds**  
- **⬤ SNP distance heatmap generation for genomic relatedness assessment**  
- **⬤ Lineage distribution summary & visualization**  
- **⬤ Optional recombination filtering using Gubbins**  
- **⬤ Maximum-likelihood phylogenetic inference using IQ-TREE2**  
- **⬤ ETE3-based phylogenetic tree visualization with lineage & resistance metadata**  
- **⬤ Downloadable QC filtering rationale & surveillance metadata TSV outputs**  
- **⬤ Integrated interactive HTML report suitable for GitHub Pages deployment**  
- **⬤ Dockerized modular WDL/Cromwell execution for reproducible analysis**

## Repository structure

```text
rMAP-TB/
├── README.md
├── LICENSE
├── .dockstore.yml
├── .gitignore
├── rMAP_TB.wdl
├── examples/
│   └── inputs.example.json
├── resources/
│   ├── adapters.fa
│   ├── H37Rv.gb
│   └── README.md
└── docs/
    ├── index.html
    ├── DEPLOYMENT.md
    ├── reports/
    │   ├── small_dataset/
    │   ├── medium_dataset/
    │   └── large_dataset/
    └── assets/
        ├── workflow/
        ├── images/
        └── css/
```

## Requirements

| Requirement | Purpose |
|---|---|
| Java | Required to run the Cromwell workflow engine |
| Cromwell | Executes the WDL workflow locally or on supported backends |
| Docker | Runs the containerized bioinformatics tools used by each WDL task |
| Paired-end Illumina FASTQ files | Primary input sequencing data for trimming, QC, species typing, TB-Profiler & variant calling |
| Adapter FASTA file | Required for Trimmomatic adapter trimming |
| MTBC GenBank reference | Required for Snippy reference-guided variant calling & Snippy-core alignment |
| Kraken2/Bracken Mycobacteria database | Required for Mycobacteria species typing; embedded in the workflow Docker image if using the recommended container |
| TB-Profiler database | Required for MTBC lineage & drug-resistance profiling; provided within the TB-Profiler container |
| Sufficient local compute resources | Needed for read processing, variant calling, SNP alignment, recombination filtering, phylogeny & HTML report generation |

## Main workflow inputs

| Input | Description |
|---|---|
| `input_reads` | Array of paired-end Illumina FASTQ files, ordered as R1 followed immediately by the matching R2 file |
| `adapters` | Adapter FASTA file used by Trimmomatic during read trimming |
| `mtbc_reference_genbank` | MTBC reference genome in GenBank format for Snippy variant calling & core-SNP alignment |
| `do_trimming` | Enables adapter trimming & read preprocessing |
| `do_quality_control` | Enables FastQC quality assessment & MultiQC aggregation |
| `do_species_typing` | Enables Mycobacteria species typing using Kraken2 + Bracken |
| `do_tb_profiler` | Enables TB-Profiler-based MTBC species, lineage, sub-lineage & drug-resistance profiling |
| `do_phylogeny` | Enables MTBC-only SNP phylogenomics, including Snippy, Snippy-core, IQ-TREE2 & tree visualization |
| `use_gubbins` | Enables optional recombination filtering before phylogenetic reconstruction |
| `tbprofiler_docker` | Docker image used for TB-Profiler AMR & lineage profiling |
| `species_typing_docker` | Docker image used for Kraken2 + Bracken Mycobacteria species typing |
| `snippy_reference_type` | Reference format used by Snippy; use `genbank` when providing a GenBank reference |
| `iqtree2_model` | IQ-TREE2 nucleotide substitution model used for maximum-likelihood phylogeny |
| `iqtree2_bootstraps` | Number of bootstrap replicates used for phylogenetic support estimation |
| `min_mtbc_samples_for_tree` | Minimum number of MTBC-positive samples required to proceed with tree reconstruction |
| `likely_transmission_snp_threshold` | SNP-distance threshold for identifying genomically close sample pairs requiring epidemiological review |
| `possible_transmission_snp_threshold` | SNP-distance threshold for identifying intermediate-distance sample pairs requiring metadata review |
| `tb_drug_resistance_genes` | Comma-separated list of TB drug-resistance-associated genes used for non-synonymous mutation reporting |
| `tree_title` | Title displayed on the rendered MTBC phylogenetic tree |
| `tree_image_format` | Output format for the ETE3-rendered phylogenetic tree image |

## Example input JSON

An example Cromwell input file is provided here:

```text
examples/inputs.example.json
```

The input FASTQ files must be ordered like this:

```json
  "rMAP_TB.input_reads": [
    "~/sample1_1.fastq.gz",
    "~/sample1_2.fastq.gz",
    "~/sample2_1.fastq.gz",
    "~/sample2_2.fastq.gz"
  ],

  "rMAP_TB.adapters": "~/adapters.fa",
  "rMAP_TB.mtbc_reference_genbank": "~/H37Rv.gb",

  "rMAP_TB.do_trimming": true,
  "rMAP_TB.do_quality_control": true,
  "rMAP_TB.do_species_typing": true,
  "rMAP_TB.do_tb_profiler": true,
  "rMAP_TB.do_phylogeny": true,
  "rMAP_TB.use_gubbins": true,

  "rMAP_TB.tbprofiler_docker": "staphb/tbprofiler:6.6.6",
  "rMAP_TB.species_typing_docker": "gmboowa/mycobacterium-kraken2-bracken:2026.05",
  "rMAP_TB.snippy_reference_type": "genbank",

  "rMAP_TB.iqtree2_model": "GTR+G",
  "rMAP_TB.iqtree2_bootstraps": 1000,
  "rMAP_TB.min_mtbc_samples_for_tree": 3,

  "rMAP_TB.likely_transmission_snp_threshold": 5,
  "rMAP_TB.possible_transmission_snp_threshold": 12,

  "rMAP_TB.report_nonsynonymous_drug_gene_mutations": true,
  "rMAP_TB.tb_drug_resistance_genes": "rpoB,katG,inhA,fabG1,ahpC,embB,pncA,rpsL,rrs,gyrA,gyrB,eis,ethA,ethR,thyA,folC,alr,ddl,gidB,tlyA,rrl,atpE,rv0678,pepQ",

  "rMAP_TB.max_cpus": 8,
  "rMAP_TB.max_memory_gb": 16,
  "rMAP_TB.min_read_length": 50,
  "rMAP_TB.min_mapping_quality": 20,

  "rMAP_TB.tree_title": "MTBC Core-SNP Phylogeny",
  "rMAP_TB.tree_width": 2400,
  "rMAP_TB.tree_height": 1600,
  "rMAP_TB.tree_image_format": "png"
}
```

## Docker Images Used in the Workflow

| Workflow Component | Docker Image | Purpose |
|---|---|---|
| Read trimming | `quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2` | Adapter trimming & read-quality filtering |
| FastQC | `staphb/fastqc:0.11.9` | Per-sample read-level quality-control assessment |
| MultiQC | `ewels/multiqc:latest` | Aggregated QC reporting across samples |
| Species typing | `gmboowa/mycobacterium-kraken2-bracken:2026.05` | *Mycobacterium* species identification using Kraken2 & Bracken |
| TB-Profiler | `staphb/tbprofiler:6.6.6` | MTBC species, lineage, sub-lineage, drug-resistance prediction & mutation-level resistance evidence |
| Snippy | `staphb/snippy:4.6.0` | Reference-guided per-sample SNP calling |
| Snippy-core | `staphb/snippy:4.6.0` | Core-genome SNP alignment generation across MTBC-positive samples |
| Non-synonymous mutation summary | `python:3.11-slim` | Extraction & reporting of non-synonymous mutations in TB drug-resistance-associated genes |
| Pairwise SNP distance & clustering | `python:3.11-slim` | Pairwise SNP distance estimation, reference-sequence exclusion & SNP cluster interpretation |
| Surveillance summary visuals | `python:3.11-slim` | Lineage distribution plots, SNP heatmap generation, QC filtering rationale & surveillance metadata TSV export |
| Gubbins | `staphb/gubbins:3.4.1` | Optional recombination filtering before phylogenetic reconstruction |
| IQ-TREE2 | `gmboowa/iqtree2-python:2.3.4` | Maximum-likelihood phylogenetic inference with bootstrap support |
| Tree visualization | `gmboowa/ete3-render:1.18` | ETE3-based phylogenetic tree rendering with lineage, resistance & bootstrap metadata |
| Report merging | `python:3.11-slim` | Final integrated interactive HTML report generation |

## Running the workflow

From the repository root:

```bash
java -jar cromwell-<version>.jar run rMAP_TB.wdl --inputs ~/inputs.example.json
```

For example:

```bash
java -jar cromwell-92.jar run rMAP_TB.wdl --inputs ~/inputs.example.json
```

## Recommended local Docker resources

For small to moderate MTBC datasets on a local workstation:

```text
CPUs:   8
Memory: 16 GB or higher
```

For larger datasets, especially when using Gubbins & IQ-TREE2, consider increasing memory & CPU allocation where possible.

## Main outputs


### Quality control & read preprocessing

- **⬤ Trimmed paired-end FASTQ files**  
- **⬤ FastQC per-sample HTML reports**  
- **⬤ FastQC ZIP output files**  
- **⬤ MultiQC aggregated quality-control report**  
- **⬤ Trimming summary table**  
- **⬤ QC summary HTML report**  

### Mycobacteria species typing

- **⬤ Kraken2 classification outputs**  
- **⬤ Kraken2 species-level reports**  
- **⬤ Bracken abundance outputs**  
- **⬤ Mycobacteria species typing TSV summary**  
- **⬤ Species typing HTML report**  
- **⬤ Most probable species call per sample**  
- **⬤ Evidence supporting species assignment**  

### TB-Profiler, lineage & MTBC filtering

- **⬤ TB-Profiler JSON outputs**  
- **⬤ TB-Profiler text reports**  
- **⬤ Combined TB-Profiler HTML report**  
- **⬤ TB-Profiler summary TSV**  
- **⬤ MTBC species, lineage & sub-lineage summary**  
- **⬤ WHO-aligned TB drug-resistance profile summary**  
- **⬤ Predicted resistant drugs summary**  
- **⬤ TB-Profiler mutation-level resistance evidence TSV**  
- **⬤ TB-Profiler mutation-level resistance evidence HTML report**  
- **⬤ MTBC-positive sample list**  
- **⬤ MTBC-filtered FASTQ files for downstream phylogenomics**  
- **⬤ MTBC selection/exclusion rationale**  

### Variant calling & core-SNP alignment

- **⬤ Per-sample Snippy variant-calling directories**  
- **⬤ Per-sample VCF files**  
- **⬤ Per-sample aligned FASTA files**  
- **⬤ Per-sample Snippy tabular variant files**  
- **⬤ Snippy logs**  
- **⬤ Variant summary HTML report**  
- **⬤ Mean-depth summary TSV**  
- **⬤ Snippy-core full alignment**  
- **⬤ Snippy-core SNP alignment**  
- **⬤ Core SNP VCF**  
- **⬤ Core SNP tabular summary**  

### Drug-resistance-associated mutation summaries

- **⬤ Non-synonymous mutation TSV summary**  
- **⬤ Non-synonymous mutation HTML report**  
- **⬤ Per-sample collapsible mutation summaries**  
- **⬤ Drug-resistance-associated gene-level mutation reporting**  

### Pairwise SNP distance & cluster interpretation

- **⬤ Pairwise SNP distance matrix TSV**  
- **⬤ Pairwise SNP distance pairs TSV**  
- **⬤ SNP cluster summary TSV**  
- **⬤ SNP distance cluster HTML report**  
- **⬤ Pairwise SNP heatmap PNG**  
- **⬤ Reference/non-sample sequence exclusion log**  
- **⬤ SNP distance task status log**  

### Surveillance summary outputs

- **⬤ Lineage distribution TSV**  
- **⬤ Lineage distribution SVG plot**  
- **⬤ SNP distance heatmap SVG**  
- **⬤ QC filtering rationale TSV**  
- **⬤ Surveillance metadata TSV**  
- **⬤ Surveillance summary HTML report**  
- **⬤ Mean depth, MTBC support, lineage, resistance profile & tree-inclusion metadata**  

### Optional recombination filtering

- **⬤ Gubbins recombination-filtered polymorphic-sites alignment**  
- **⬤ Gubbins recombination-filtered final tree**  
- **⬤ Gubbins log files**  
- **⬤ Recombination-filtering status outputs**  

### Phylogenetic inference & visualization

- **⬤ IQ-TREE2 maximum-likelihood tree file**  
- **⬤ IQ-TREE2 report file**  
- **⬤ IQ-TREE2 log file**  
- **⬤ Bootstrap-supported Newick tree**  
- **⬤ Exportable Newick tree for downstream visualization tools such as iTOL**  
- **⬤ ETE3-rendered MTBC phylogenetic tree image**  
- **⬤ Cleaned tree file used for visualization**  
- **⬤ Tree rendering log**  

### Integrated reports & downloadable public-health outputs

- **⬤ Final integrated interactive HTML report**  
- **⬤ Run metadata file**  
- **⬤ Downloadable TB surveillance metadata TSV**  
- **⬤ Downloadable QC filtering rationale TSV**  
- **⬤ Embedded lineage distribution plot**  
- **⬤ Embedded SNP distance heatmap**  
- **⬤ Embedded MTBC phylogenetic tree**  
- **⬤ GitHub Pages-compatible report outputs**  

Example final report output:

```text
integrated_report.html

```
## GitHub Pages report site

```text
https://gmboowa.github.io/rMAP_TB/
```


## Interpretation guidance

The integrated report should be interpreted using multiple complementary layers of genomic, resistance, quality-control & epidemiological evidence:

- **⬤ Mycobacteria species typing results**  
- **⬤ MTBC selection & filtering rationale**  
- **⬤ TB-Profiler species, lineage & sub-lineage calls**  
- **⬤ TB-Profiler drug-resistance profile**  
- **⬤ Mutation-level resistance evidence, including drug, gene, mutation/change, confidence & evidence fields**  
- **⬤ Non-synonymous mutations in key TB drug-resistance-associated genes**  
- **⬤ Mean depth & sample-level QC indicators**  
- **⬤ Pairwise SNP distances between MTBC isolates**  
- **⬤ SNP cluster interpretation using the configured SNP-distance thresholds**  
- **⬤ SNP distance heatmap for genomic relatedness assessment**  
- **⬤ Core-SNP phylogenetic clustering**  
- **⬤ Bootstrap support values on the phylogenetic tree**  
- **⬤ Recombination-filtered alignment & tree, if Gubbins is enabled**  
- **⬤ Surveillance metadata, including country, year, collection site, sample source, lineage, resistance profile & tree-inclusion status where available**

Close SNP clustering or close placement on a phylogenetic tree should **not** be interpreted as proof of direct transmission on its own. Transmission interpretation should be made only after considering epidemiological linkage, sampling density, collection dates, geography, lineage, resistance profile, sequence quality, SNP distances & bootstrap support.

## Suggested citation / acknowledgement

If you use this workflow, please cite or acknowledge the associated manuscript:

rMAP-TB: a reproducible WDL/Cromwell workflow for *Mycobacterium tuberculosis* complex genomic surveillance and drug-resistance interpretation.

## License

▪ MIT License for permissive open-source reuse  

