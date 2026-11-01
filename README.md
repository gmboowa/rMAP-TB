rMAP-TB: Rapid Mycobacterial Analysis Pipeline for Tuberculosis & Mycobacterial Genomic Surveillance
---
rMAP-TB is a modular, reproducible, Dockerized WDL/Cromwell bioinformatics workflow for public-health-oriented genomic analysis of Mycobacterium tuberculosis complex (MTBC) and non-MTBC mycobacteria. It supports paired-end Illumina FASTQ inputs and integrates read preprocessing, sequence quality control, species identification, MTBC/NTM sample routing, tuberculosis drug-resistance prediction, species and lineage profiling, MTBC-only sample filtering, core-SNP phylogenomics, SNP distance analysis, optional phylogenetic reconstruction, and interactive genomic surveillance reporting.

The workflow first performs read trimming, sequence quality control, and Kraken2/Bracken-based Mycobacterium species identification. Species typing is used to route samples before TB-Profiler execution: MTBC-supported samples proceed to TB-Profiler for drug-resistance prediction, species identification, and lineage profiling, while non-MTBC mycobacteria are summarized separately through a dedicated nontuberculous mycobacteria (NTM) speciation branch. This approach enables the workflow to report the most probable NTM species without subjecting non-MTBC samples to MTBC-specific resistance and lineage analyses.

For MTBC-supported samples, rMAP-TB integrates reference-based read alignment, variant calling, core-genome SNP identification, and SNP distance analysis to support genomic relatedness assessments and molecular epidemiological investigations. Optional phylogenomic reconstruction enables the characterization of evolutionary relationships among isolates, supporting investigations of potential transmission clusters and population structure.

The pipeline combines antimicrobial resistance prediction, taxonomic identification, phylogenomic analysis, and genomic surveillance reporting within a standardized, portable, and reproducible computational framework. Containerized execution through Docker and workflow orchestration using WDL/Cromwell facilitate consistent analyses across compatible local, institutional, and cloud computing environments.

## Visual Overview of the rMAP-TB Workflow

The rMAP-TB pipeline integrates sequencing quality control, Mycobacterium species identification, MTBC/NTM sample routing, tuberculosis drug-resistance prediction, core-SNP analysis, phylogenomic reconstruction, and genomic surveillance reporting within a reproducible WDL/Cromwell framework.

![rMAP-TB Genomic Surveillance Workflow](docs/assets/workflow/rMAP-TB-workflow.png)

*Figure 1. Overview of the rMAP-TB workflow, showing paired-end Illumina sequencing inputs, read preprocessing and quality control, Kraken2/Bracken-based species identification, MTBC/NTM sample routing, TB-Profiler-based drug-resistance and lineage profiling for MTBC-supported samples, independent NTM species reporting, core-SNP analysis, SNP distance estimation, optional phylogenomic reconstruction, and integrated interactive genomic surveillance reporting.*
