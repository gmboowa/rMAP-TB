## TB-AMR-MTBC-Phylogenomics

A reproducible, containerized WDL/Cromwell workflow for *Mycobacterium tuberculosis* complex (MTBC) antimicrobial-resistance profiling & core-SNP phylogenomics.

The workflow supports paired-end Illumina FASTQ inputs & produces quality-control summaries, TB-Profiler drug-resistance & lineage reports, MTBC-only sample filtering, Snippy-based core-SNP outputs, optional Gubbins recombination filtering, IQ-TREE2 maximum-likelihood phylogeny, tree visualization & interactive HTML reports.


## Workflow overview
<div align="center">
<pre>

Paired-end FASTQ files
        ⬇
Read trimming
        ⬇
FastQC + MultiQC
        ⬇
TB-Profiler AMR, species & lineage profiling
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

</pre>
</div>

## Key features

- **⬤ Paired-end FASTQ input support**  
- **⬤ Adapter trimming before downstream analysis**  
- **⬤ FastQC & MultiQC quality-control reporting**  
- **⬤ TB-Profiler-based drug-resistance prediction**  
- **⬤ MTBC lineage & sub-lineage reporting**  
- **⬤ MTBC-only filtering before phylogenomic reconstruction**  
- **⬤ Snippy-based reference-guided variant calling**  
- **⬤ Core-SNP alignment generation using Snippy-core**  
- **⬤ Optional recombination filtering using Gubbins**  
- **⬤ Maximum-likelihood phylogeny using IQ-TREE2**  
- **⬤ Optional midpoint-rooted tree visualization**  
- **⬤ Interactive HTML reports suitable for GitHub Pages**  
- **⬤ Non-synonymous mutation reporting for key TB drug-resistance genes**

## Repository structure

```text
TB-AMR-MTBC-Phylogenomics/
├── README.md
├── TB.wdl
├── examples/
│   └── inputs.example.json
└── docs/
    ├── index.html
    ├── DEPLOYMENT.md
    ├── reports/
    │   ├── integrated_tb_report.html
    │   └── tbprofiler_report.html
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
| MTBC GenBank reference | Required for Snippy & phylogenomic analysis |

## Main workflow inputs

| Input | Description |
|---|---|
| `input_reads` | Array of paired-end FASTQ files, ordered as R1 followed immediately by matching R2 |
| `adapters` | Adapter FASTA file used by the trimming step |
| `mtbc_reference_genbank` | MTBC reference genome in GenBank format |
| `do_trimming` | Enables read trimming |
| `do_quality_control` | Enables FastQC & MultiQC |
| `do_tb_profiler` | Enables TB-Profiler AMR & lineage analysis |
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
    "~/sample1_1.fastq.gz",
    "~/sample1_2.fastq.gz",
    "~/sample2_1.fastq.gz",
    "~/sample2_2.fastq.gz"
  ],

  "TB_AMR_MTBC_Phylogenomics.adapters": "~/adapters.fa",
  "TB_AMR_MTBC_Phylogenomics.mtbc_reference_genbank": "~/H37Rv.gb",

  "TB_AMR_MTBC_Phylogenomics.do_trimming": true,
  "TB_AMR_MTBC_Phylogenomics.do_quality_control": true,
  "TB_AMR_MTBC_Phylogenomics.do_tb_profiler": true,
  "TB_AMR_MTBC_Phylogenomics.do_phylogeny": true,
  "TB_AMR_MTBC_Phylogenomics.use_gubbins": true,

  "TB_AMR_MTBC_Phylogenomics.tbprofiler_docker": "staphb/tbprofiler:6.6.6",
  "TB_AMR_MTBC_Phylogenomics.snippy_reference_type": "genbank",

  "TB_AMR_MTBC_Phylogenomics.iqtree2_model": "GTR+G",
  "TB_AMR_MTBC_Phylogenomics.iqtree2_bootstraps": 1000,
  "TB_AMR_MTBC_Phylogenomics.min_mtbc_samples_for_tree": 3,

  "TB_AMR_MTBC_Phylogenomics.report_nonsynonymous_drug_gene_mutations": true,
  "TB_AMR_MTBC_Phylogenomics.tb_drug_resistance_genes":   "rpoB,katG,inhA,fabG1,ahpC,embB,pncA,rpsL,rrs,gyrA,gyrB,eis,ethA,ethR,thyA,folC,alr,ddl,gidB,tlyA,rrl,atpE,rv0678,pepQ",

  "TB_AMR_MTBC_Phylogenomics.max_cpus": 8,
  "TB_AMR_MTBC_Phylogenomics.max_memory_gb": 16,
  "TB_AMR_MTBC_Phylogenomics.min_read_length": 50,
  "TB_AMR_MTBC_Phylogenomics.min_mapping_quality": 20,

  "TB_AMR_MTBC_Phylogenomics.tree_width": 2400,
  "TB_AMR_MTBC_Phylogenomics.tree_height": 1600,
  "TB_AMR_MTBC_Phylogenomics.tree_image_format": "png"
}
```

## Running the workflow

From the repository root:

```bash
java -jar cromwell-<version>.jar run TB.wdl --inputs ~/inputs.example.json
```

For example:

```bash
java -jar cromwell-92.jar run TB.wdl --inputs ~/inputs.example.json
```

## Recommended local Docker resources

For small to moderate MTBC datasets on a local workstation:

```text
CPUs:   8
Memory: 16 GB or higher
```

For larger datasets, especially when using Gubbins & IQ-TREE2, consider increasing memory & CPU allocation where possible.

## Main outputs

### Quality control

- **⬤ Trimmed FASTQ files**  
- **⬤ FastQC reports**  
- **⬤ MultiQC report**  
- **⬤ Trimming summary table**

### TB-Profiler & MTBC filtering

- **⬤ TB-Profiler JSON outputs**  
- **⬤ Combined TB-Profiler HTML report**  
- **⬤ Resistance profile summary**  
- **⬤ MTBC-positive sample list**  
- **⬤ Non-MTBC or low-confidence excluded sample list**

### SNP phylogenomics

- **⬤ Per-sample Snippy variant calls**  
- **⬤ Snippy-core alignment**  
- **⬤ Core SNP VCF**  
- **⬤ Optional Gubbins recombination-filtered alignment**  
- **⬤ IQ-TREE2 maximum-likelihood tree**  
- **⬤ Bootstrap-supported tree file**  
- **⬤ Tree visualization as PNG or SVG**

### Interactive reports

▪ `integrated_report(s).html`  


## GitHub Pages report site

```text
https://gmboowa.github.io/TB-AMR-MTBC-Phylogenomics/
```


## Interpretation guidance

The integrated report should be interpreted using multiple layers of evidence:

- **⬤ TB-Profiler resistance profile**  
- **⬤ MTBC lineage & sub-lineage**  
- **⬤ Core-SNP phylogenetic clustering**  
- **⬤ Bootstrap support**  
- **⬤ Recombination-filtered alignment (if Gubbins is enabled)**  
- **⬤ Metadata (country, year, collection site, sample source)**

Close clustering alone should not be treated as proof of transmission without epidemiological & sampling context.

## Suggested citation / acknowledgement

If you use this workflow, please cite or acknowledge:

```text
TB-AMR-MTBC-Phylogenomics: a WDL/Cromwell workflow for MTBC antimicrobial-resistance profiling & core-SNP phylogenomics.
```

## License

▪ MIT License for permissive open-source reuse  

