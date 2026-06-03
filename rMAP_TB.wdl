version 1.0

workflow rMAP_TB {
  input {
    Array[File]+ read1
    Array[File]+ read2
    File? sample_metadata_tsv
    File adapters
    File mtbc_reference_genbank

    Boolean do_trimming = true
    Boolean do_quality_control = true
    Boolean do_species_typing = true
    Boolean do_tb_profiler = true
    Boolean do_phylogeny = true
    Boolean use_gubbins = true
    Boolean midpoint_root_tree = true
    Boolean report_nonsynonymous_drug_gene_mutations = true
    Boolean report_snp_distance_clusters = true
    Boolean generate_surveillance_summary_visuals = true

    String trimmomatic_quality_encoding = "phred33"
    String species_typing_docker = "gmboowa/mycobacterium-kraken2-bracken:2026.05"
    String tbprofiler_docker = "staphb/tbprofiler:6.6.6"
    String snippy_reference_type = "genbank"
    String iqtree2_model = "GTR+G"
    Int iqtree2_bootstraps = 1000
    Int min_mtbc_samples_for_tree = 3
    Int likely_transmission_snp_threshold = 5
    Int possible_transmission_snp_threshold = 12

    # Terra runtime controls. Defaults increased for faster execution of TB-Profiler, Snippy, Gubbins, and IQ-TREE.
    # Lower these values if cost control is more important than speed.
    Int max_cpus = 16
    Int max_memory_gb = 64
    Int min_read_length = 50
    Int min_mapping_quality = 20
    Int tree_width = 2400
    Int tree_height = 1600
    String tree_image_format = "png"
    String tb_drug_resistance_genes = "rpoB,katG,inhA,fabG1,ahpC,embB,pncA,rpsL,rrs,gyrA,gyrB,eis,ethA,ethR,thyA,folC,alr,ddl,gidB,tlyA,rrl,atpE,rv0678,pepQ"
  }

  Int cpu_4 = if max_cpus < 4 then max_cpus else 4
  Int cpu_8 = if max_cpus < 8 then max_cpus else 8
  Int cpu_16 = if max_cpus < 16 then max_cpus else 16
  Int min_mtbc_fastq_files_for_tree = min_mtbc_samples_for_tree * 2

  # 0. Prepare paired-end reads supplied from Terra data tables.
  # This follows the same Terra input logic as the WHO_2021 optional-metadata version.
  #
  # For a Terra SET run, use an entity table named TBs and a set table named TBs_set.
  # In the workflow input form, launch from TBs_set and set:
  #   rMAP_TB.read1 = this.TBs.read1
  #   rMAP_TB.read2 = this.TBs.read2
  #
  # The expressions above resolve to Array[File]+ values only when the selected
  # root entity is a set containing multiple TBs members. This task then validates
  # matching R1/R2 counts and creates an interleaved Array[File] for downstream
  # tasks that expect paired reads in R1/R2 order.
  call PREPARE_READ_PAIRS {
    input:
      read1 = read1,
      read2 = read2
  }

  # 1. Trimming is the first executable step. All downstream analysis uses these reads when trimming is enabled.
  if (do_trimming) {
    call TRIMMING {
      input:
        input_reads = PREPARE_READ_PAIRS.paired_reads,
        adapters = adapters,
        trimmomatic_quality_encoding = trimmomatic_quality_encoding,
        cpu = cpu_8,
        min_length = min_read_length
    }
  }

  Array[File] analysis_reads = select_first([TRIMMING.trimmed_reads, PREPARE_READ_PAIRS.paired_reads])

  # 2. Sample QC and Trimming Summary. FastQC runs only after trimming, because it consumes analysis_reads.
  if (do_quality_control) {
    call FASTQC {
      input:
        input_reads = analysis_reads,
        cpu = cpu_8
    }

    call MULTIQC {
      input:
        fastqc_reports = FASTQC.fastqc_reports,
        fastqc_zips = FASTQC.fastqc_zips
    }
  }

  # 3. Species Typing using Kraken2 + Bracken.
  # This step runs after QC/trimming and before TB-Profiler.
  # It uses the custom Mycobacterium-only database image and reports only the most probable species call per sample.
  if (do_species_typing) {
    call SPECIES_TYPING {
      input:
        input_reads = analysis_reads,
        qc_dependency = MULTIQC.multiqc_report,
        docker_image = species_typing_docker,
        cpu = cpu_16,
        memory_gb = max_memory_gb
    }

    call MYCOBACTERIA_SAMPLE_ROUTER {
      input:
        input_reads = analysis_reads,
        species_typing_tsv = SPECIES_TYPING.species_typing_tsv
    }
  }

  Array[File] mtbc_reads = select_first([MYCOBACTERIA_SAMPLE_ROUTER.mtbc_reads, []])
  Int mtbc_sample_count_from_kraken = select_first([MYCOBACTERIA_SAMPLE_ROUTER.mtbc_sample_count, 0])

  # 4. TB-Profiler Resistance, Species, Lineage, and Mutation-Level Evidence Report.
  # TB-Profiler is intentionally chained after Species Typing and receives only Kraken2/Bracken-supported MTBC reads.
  # If the dataset contains only non-MTBC/NTM samples, TB-Profiler is skipped and the final report summarizes NTM speciation.
  if (do_tb_profiler && size(mtbc_reads) >= 2 && mtbc_sample_count_from_kraken > 0) {
    call TB_PROFILER_AND_MTBC_FILTER {
      input:
        input_reads = mtbc_reads,
        qc_dependency = SPECIES_TYPING.species_typing_html,
        species_typing_tsv = SPECIES_TYPING.species_typing_tsv,
        docker_image = tbprofiler_docker,
        cpu = cpu_16,
        memory_gb = max_memory_gb
    }
  }

  # 5. MTBC-only Snippy/core-SNP analysis.
  # If fewer than min_mtbc_samples_for_tree paired MTBC samples are available, this branch is skipped.
  if (do_tb_profiler && do_phylogeny && size(mtbc_reads) >= min_mtbc_fastq_files_for_tree) {
    call SNIPPY_CORE_MTBC {
      input:
        input_reads = mtbc_reads,
        reference_genome = mtbc_reference_genbank,
        reference_type = snippy_reference_type,
        cpu = cpu_16,
        memory_gb = max_memory_gb,
        min_quality = min_mapping_quality
    }
  }

  # 6. Non-synonymous Mutation Summary.
  if (report_nonsynonymous_drug_gene_mutations && defined(SNIPPY_CORE_MTBC.snippy_tab_files)) {
    call TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS {
      input:
        snippy_tab_files = select_first([SNIPPY_CORE_MTBC.snippy_tab_files, []]),
        genes_csv = tb_drug_resistance_genes
    }
  }

  # 7. Pairwise SNP Distance and Cluster Summary.
  # This uses the MTBC-only Snippy core.full.aln before optional Gubbins filtering.
  if (do_phylogeny && report_snp_distance_clusters && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call SNP_DISTANCE_CLUSTERING {
      input:
        core_full_alignment = select_first([SNIPPY_CORE_MTBC.core_full_alignment]),
        likely_transmission_snp_threshold = likely_transmission_snp_threshold,
        possible_transmission_snp_threshold = possible_transmission_snp_threshold,
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

  # 8. Surveillance summary visuals and metadata.
  # This generates lineage distribution, SNP distance heatmap, QC filtering rationale, and surveillance metadata TSV.
  if (do_tb_profiler && generate_surveillance_summary_visuals && defined(TB_PROFILER_AND_MTBC_FILTER.summary_tsv)) {
    call TB_SURVEILLANCE_SUMMARY_VISUALS {
      input:
        tbprofiler_summary_tsv = select_first([TB_PROFILER_AND_MTBC_FILTER.summary_tsv]),
        species_typing_tsv = SPECIES_TYPING.species_typing_tsv,
        pairwise_snp_distance_matrix = SNP_DISTANCE_CLUSTERING.pairwise_snp_distance_matrix,
        mean_depth_tsv = SNIPPY_CORE_MTBC.mean_depth_summary_tsv,
        sample_metadata_tsv = sample_metadata_tsv,
        cpu = cpu_4,
        memory_gb = max_memory_gb
    }
  }

  # 9. Optional recombination filtering. If Gubbins fails internally, its task passes the original alignment forward.
  if (do_phylogeny && use_gubbins && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call GUBBINS_RECOMBINATION {
      input:
        core_full_alignment = select_first([SNIPPY_CORE_MTBC.core_full_alignment]),
        cpu = cpu_16,
        memory_gb = max_memory_gb
    }
  }

  # 10. IQ-TREE uses the Gubbins-filtered alignment when available; otherwise it uses Snippy core.full.aln.
  #     Before IQ-TREE runs, samples with excessive missing, ambiguous, or gap-only
  #     sequence content are removed from the alignment and written to
  #     excluded_from_iqtree.tsv for transparent reporting in the final HTML.
  if (do_phylogeny && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call IQTREE2_PHYLOGENY {
      input:
        alignment = select_first([GUBBINS_RECOMBINATION.filtered_alignment, SNIPPY_CORE_MTBC.core_full_alignment]),
        model = iqtree2_model,
        bootstrap_replicates = iqtree2_bootstraps,
        cpu = cpu_16,
        memory_gb = max_memory_gb,
        midpoint_root_tree = midpoint_root_tree,
        max_missing_fraction_for_tree = 0.50,
        min_non_reference_samples_for_tree = 3
    }
  }

  # 11. MTBC-only Core-SNP Phylogenetic Tree.
  if (do_phylogeny && defined(IQTREE2_PHYLOGENY.final_tree)) {
    call TREE_VISUALIZATION {
      input:
        input_tree = select_first([IQTREE2_PHYLOGENY.final_tree]),
        tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv,
        resistance_profile_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.resistance_profile_summary_tsv,
        iqtree_excluded_samples_tsv = IQTREE2_PHYLOGENY.excluded_from_iqtree,
        width = tree_width,
        height = tree_height,
        image_format = tree_image_format
    }
  }

  # 12. Final merged report.
  # Desired report order:
  # 1. Sample QC and Trimming Summary
  # 2. Species Typing using Kraken2 + Bracken
  # 3. TB-Profiler Resistance, Species, and Lineage Report
  # 4. Resistance Mutation Evidence Summary
  # 5. Non-synonymous Mutation Summary
  # 6. Lineage Distribution Summary
  # 7. Pairwise SNP Distance and Cluster Summary
  # 8. SNP Distance Heatmap
  # 9. MTBC-only Core-SNP Phylogenetic Tree
  # 10. QC Filtering Rationale and Surveillance Metadata
  # 11. Pipeline Provenance and Software Versions
  call MERGE_TB_REPORTS {
  input:
    tbprofiler_html = TB_PROFILER_AND_MTBC_FILTER.combined_html,
    tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv,
    resistance_profile_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.resistance_profile_summary_tsv,
    tbprofiler_mutation_evidence_tsv = TB_PROFILER_AND_MTBC_FILTER.mutation_evidence_tsv,
    tbprofiler_mutation_evidence_html = TB_PROFILER_AND_MTBC_FILTER.mutation_evidence_html,
    mtbc_samples_txt = TB_PROFILER_AND_MTBC_FILTER.mtbc_samples_txt,

    species_typing_html = SPECIES_TYPING.species_typing_html,
    species_typing_tsv = SPECIES_TYPING.species_typing_tsv,
    ntm_species_summary_html = MYCOBACTERIA_SAMPLE_ROUTER.ntm_species_summary_html,
    ntm_species_summary_tsv = MYCOBACTERIA_SAMPLE_ROUTER.ntm_species_summary_tsv,
    mycobacteria_routing_summary_tsv = MYCOBACTERIA_SAMPLE_ROUTER.mycobacteria_routing_summary_tsv,
    mycobacteria_routing_status = MYCOBACTERIA_SAMPLE_ROUTER.routing_status_txt,

    qc_summary_html = MULTIQC.multiqc_report,
    trimming_report_html = TRIMMING.trimming_report,
    variant_summary_html = SNIPPY_CORE_MTBC.variant_summary,

    iqtree_report = IQTREE2_PHYLOGENY.iqtree_report,
    iqtree_excluded_samples_tsv = IQTREE2_PHYLOGENY.excluded_from_iqtree,
    iqtree_included_samples_tsv = IQTREE2_PHYLOGENY.included_in_iqtree,
    iqtree_filtering_summary_txt = IQTREE2_PHYLOGENY.alignment_filtering_summary,
    iqtree_status = IQTREE2_PHYLOGENY.iqtree_status,

    tree_image = TREE_VISUALIZATION.tree_image,
    phylogenetic_tree_newick = TREE_VISUALIZATION.cleaned_tree,

    nonsynonymous_mutations_tsv = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_tsv,
    nonsynonymous_mutations_html = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_html,

    pairwise_snp_distance_matrix = SNP_DISTANCE_CLUSTERING.pairwise_snp_distance_matrix,
    pairwise_snp_distance_pairs = SNP_DISTANCE_CLUSTERING.pairwise_snp_distance_pairs,
    snp_cluster_summary = SNP_DISTANCE_CLUSTERING.snp_cluster_summary,
    snp_distance_cluster_html = SNP_DISTANCE_CLUSTERING.snp_distance_cluster_html,

    lineage_distribution_tsv = TB_SURVEILLANCE_SUMMARY_VISUALS.lineage_distribution_tsv,
    lineage_distribution_svg = TB_SURVEILLANCE_SUMMARY_VISUALS.lineage_distribution_svg,
    snp_distance_heatmap_svg = TB_SURVEILLANCE_SUMMARY_VISUALS.snp_distance_heatmap_svg,
    surveillance_metadata_tsv = TB_SURVEILLANCE_SUMMARY_VISUALS.surveillance_metadata_tsv,
    qc_filtering_rationale_tsv = TB_SURVEILLANCE_SUMMARY_VISUALS.qc_filtering_rationale_tsv,
    surveillance_summary_html = TB_SURVEILLANCE_SUMMARY_VISUALS.surveillance_summary_html
}

  output {
    Array[File]? trimmed_reads = TRIMMING.trimmed_reads
    File? trimming_report = TRIMMING.trimming_report
    File? trimming_summary = TRIMMING.trimming_summary
    Array[File]? trimming_logs = TRIMMING.trimming_logs

    Array[File]? fastqc_reports = FASTQC.fastqc_reports
    Array[File]? fastqc_zips = FASTQC.fastqc_zips
    File? fastqc_summary_html = FASTQC.fastqc_summary_html
    File? fastqc_summary_tsv = FASTQC.fastqc_summary_tsv
    File? fastqc_log = FASTQC.fastqc_log
    File? multiqc_report = MULTIQC.multiqc_report
    File? multiqc_log = MULTIQC.multiqc_log

    File? species_typing_html = SPECIES_TYPING.species_typing_html
    File? species_typing_tsv = SPECIES_TYPING.species_typing_tsv
    Array[File]? species_typing_kraken_reports = SPECIES_TYPING.kraken_reports
    Array[File]? species_typing_kraken_outputs = SPECIES_TYPING.kraken_outputs

    File? mycobacteria_routing_summary_tsv = MYCOBACTERIA_SAMPLE_ROUTER.mycobacteria_routing_summary_tsv
    File? ntm_species_summary_tsv = MYCOBACTERIA_SAMPLE_ROUTER.ntm_species_summary_tsv
    File? ntm_species_summary_html = MYCOBACTERIA_SAMPLE_ROUTER.ntm_species_summary_html
    File? mtbc_species_summary_tsv = MYCOBACTERIA_SAMPLE_ROUTER.mtbc_species_summary_tsv
    File? non_mtbc_species_summary_tsv = MYCOBACTERIA_SAMPLE_ROUTER.non_mtbc_species_summary_tsv
    File? mycobacteria_routing_status = MYCOBACTERIA_SAMPLE_ROUTER.routing_status_txt
    Array[File]? kraken_selected_mtbc_reads = MYCOBACTERIA_SAMPLE_ROUTER.mtbc_reads

    Array[File]? tbprofiler_json = TB_PROFILER_AND_MTBC_FILTER.json_reports
    Array[File]? tbprofiler_txt = TB_PROFILER_AND_MTBC_FILTER.txt_reports
    File? tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv
    File? tbprofiler_combined_html = TB_PROFILER_AND_MTBC_FILTER.combined_html
    File? tbprofiler_mutation_evidence_tsv = TB_PROFILER_AND_MTBC_FILTER.mutation_evidence_tsv
    File? tbprofiler_mutation_evidence_html = TB_PROFILER_AND_MTBC_FILTER.mutation_evidence_html
    File? mtbc_samples = TB_PROFILER_AND_MTBC_FILTER.mtbc_samples_txt
    Array[File]? mtbc_reads_for_phylogeny = TB_PROFILER_AND_MTBC_FILTER.mtbc_reads

    File? variant_summary_html = SNIPPY_CORE_MTBC.variant_summary
    File? core_full_alignment = SNIPPY_CORE_MTBC.core_full_alignment
    File? core_snp_alignment = SNIPPY_CORE_MTBC.core_snp_alignment
    File? snippy_core_vcf = SNIPPY_CORE_MTBC.core_vcf
    Array[File]? snippy_tab_files = SNIPPY_CORE_MTBC.snippy_tab_files
    File? mean_depth_summary_tsv = SNIPPY_CORE_MTBC.mean_depth_summary_tsv

    File? nonsynonymous_drug_gene_mutations_tsv = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_tsv
    File? nonsynonymous_drug_gene_mutations_html = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_html

    File? pairwise_snp_distance_matrix = SNP_DISTANCE_CLUSTERING.pairwise_snp_distance_matrix
    File? pairwise_snp_distance_pairs = SNP_DISTANCE_CLUSTERING.pairwise_snp_distance_pairs
    File? snp_cluster_summary = SNP_DISTANCE_CLUSTERING.snp_cluster_summary
    File? snp_distance_cluster_html = SNP_DISTANCE_CLUSTERING.snp_distance_cluster_html

    File? lineage_distribution_tsv = TB_SURVEILLANCE_SUMMARY_VISUALS.lineage_distribution_tsv
    File? lineage_distribution_svg = TB_SURVEILLANCE_SUMMARY_VISUALS.lineage_distribution_svg
    File? snp_distance_heatmap_svg = TB_SURVEILLANCE_SUMMARY_VISUALS.snp_distance_heatmap_svg
    File? surveillance_metadata_tsv = TB_SURVEILLANCE_SUMMARY_VISUALS.surveillance_metadata_tsv
    File? qc_filtering_rationale_tsv = TB_SURVEILLANCE_SUMMARY_VISUALS.qc_filtering_rationale_tsv
    File? surveillance_summary_html = TB_SURVEILLANCE_SUMMARY_VISUALS.surveillance_summary_html

    File? gubbins_filtered_alignment = GUBBINS_RECOMBINATION.filtered_alignment
    File? gubbins_status = GUBBINS_RECOMBINATION.gubbins_status
    File? iqtree_newick = IQTREE2_PHYLOGENY.final_tree
    File? iqtree_status = IQTREE2_PHYLOGENY.iqtree_status
    File? tree_image = TREE_VISUALIZATION.tree_image
    File? tree_render_status = TREE_VISUALIZATION.render_log

    File final_merged_html_report = MERGE_TB_REPORTS.final_report_html
    File run_metadata = MERGE_TB_REPORTS.run_metadata
  }
}
task PREPARE_READ_PAIRS {
  input {
    Array[File]+ read1
    Array[File]+ read2
    String docker_image = "ubuntu:22.04"
  }

  command <<<
    set -euo pipefail
    mkdir -p prepared

    r1=(~{sep=' ' read1})
    r2=(~{sep=' ' read2})

    if [ ${#r1[@]} -ne ${#r2[@]} ]; then
      echo "ERROR: read1 and read2 must contain the same number of files." >&2
      echo "read1 count: ${#r1[@]}" >&2
      echo "read2 count: ${#r2[@]}" >&2
      exit 1
    fi

    echo -e "sample_index\tsample\tread1\tread2\tprepared_read1\tprepared_read2" > prepared_read_pairs.tsv

    for i in "${!r1[@]}"; do
      R1="${r1[$i]}"
      R2="${r2[$i]}"

      sample=$(basename "$R1")
      sample=$(echo "$sample" | sed -E 's/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//' | sed -E 's/(_R?1|_1|\.R?1|\.1)(_|$).*//')
      prefix=$(printf "%04d_%s" "$i" "$sample")

      out_R1="prepared/${prefix}_R1.fastq.gz"
      out_R2="prepared/${prefix}_R2.fastq.gz"

      ln -s "$R1" "$out_R1"
      ln -s "$R2" "$out_R2"

      echo -e "$i\t$sample\t$R1\t$R2\t$out_R1\t$out_R2" >> prepared_read_pairs.tsv
    done
  >>>

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 20 HDD"
  }

  output {
    Array[File] paired_reads = glob("prepared/*.fastq.gz")
    File prepared_read_pairs = "prepared_read_pairs.tsv"
  }
}

task TRIMMING {
  input {
    String docker_image = "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    Array[File]+ input_reads
    File adapters
    String trimmomatic_quality_encoding = "phred33"
    Int cpu = 8
    Int min_length = 50
  }

  command <<<
    set -uo pipefail
    mkdir -p trimmed logs

    files=(~{sep=' ' input_reads})
    n=${#files[@]}

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: input_reads must contain paired reads in R1/R2 order." >&2
      exit 1
    fi

    if [ ! -s "~{adapters}" ]; then
      echo "WARNING: Adapter file is missing or empty. Raw reads will be passed forward." >&2
      use_trimming="false"
    else
      use_trimming="true"
    fi

    echo -e "sample\tinput_pairs\tpaired_reads_output\tstatus" > trimming_summary.tsv

    for ((i=0; i<n; i+=2)); do
      R1="${files[$i]}"
      R2="${files[$((i+1))]}"
      sample=$(basename "$R1")
      sample=$(echo "$sample" | sed -E 's/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//' | sed -E 's/(_R?1|_1|\.R?1|\.1)(_|$).*//')

      out_R1="trimmed/${sample}_R1_paired.fastq.gz"
      out_R2="trimmed/${sample}_R2_paired.fastq.gz"

      if [ "$use_trimming" = "true" ]; then
        if trimmomatic PE \
          -threads ~{cpu} \
          -~{trimmomatic_quality_encoding} \
          "$R1" "$R2" \
          "$out_R1" "trimmed/${sample}_R1_unpaired.fastq.gz" \
          "$out_R2" "trimmed/${sample}_R2_unpaired.fastq.gz" \
          ILLUMINACLIP:~{adapters}:2:30:10 \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:~{min_length} \
          > "logs/${sample}.trimmomatic.log" 2>&1; then

          status="success"

        else
          status="trimming_failed_raw_reads_used"
          echo "WARNING: Trimmomatic failed for ${sample}; raw reads copied forward." >> "logs/${sample}.trimmomatic.log"
          cp "$R1" "$out_R1"
          cp "$R2" "$out_R2"
        fi
      else
        status="adapter_missing_raw_reads_used"
        cp "$R1" "$out_R1"
        cp "$R2" "$out_R2"
      fi

      if [ ! -s "$out_R1" ] || [ ! -s "$out_R2" ]; then
        echo "ERROR: No usable paired reads produced for ${sample}" >&2
        exit 1
      fi

      echo -e "${sample}\t${R1};${R2}\t${out_R1};${out_R2}\t${status}" >> trimming_summary.tsv
    done

    python3 - <<'PYDEPTH'
import re
from pathlib import Path

out = Path("mean_depth_summary.tsv")
out.write_text("sample\tmean_depth\tdepth_source\n", encoding="utf-8")

def sample_from_path(path):
    name = path.name
    for suffix in [".txt", ".log", ".vcf"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name

def extract_depth(text):
    patterns = [
        r"Mean\s+depth\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
        r"Average\s+depth\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
        r"Read\s+depth\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
        r"coverage\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
        r"mean\s+coverage\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            return m.group(1)
    return "Not available"

seen = set()
with out.open("a", encoding="utf-8") as fh:
    for txt in sorted(Path("snippy_results").glob("*/*.txt")):
        sample = txt.parent.name
        depth = extract_depth(txt.read_text(errors="replace"))
        fh.write(f"{sample}\t{depth}\t{txt}\n")
        seen.add(sample)
    for log in sorted(Path("logs").glob("*.snippy.log")):
        sample = sample_from_path(log)
        if sample in seen:
            continue
        depth = extract_depth(log.read_text(errors="replace"))
        fh.write(f"{sample}\t{depth}\t{log}\n")
        seen.add(sample)
PYDEPTH

    python3 - <<'PY'
import csv, html

rows = list(csv.DictReader(open("trimming_summary.tsv"), delimiter="\t"))

out = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Trimming summary</title>
<style>
body{font-family:Arial;margin:24px;background:#f8fafc;color:#111827}
table{border-collapse:collapse;width:100%;background:white}
th,td{border:1px solid #dbe4ee;padding:8px;vertical-align:top}
th{background:#0f766e;color:white}
.ok{color:#166534;font-weight:bold}
.warn{color:#d97706;font-weight:bold}
</style>
</head>
<body>
<h1>Read trimming summary</h1>
<p>Adapter and quality trimming was attempted with Trimmomatic. If trimming failed for a sample, raw reads were copied forward so downstream analysis could continue.</p>
<table>
<thead>
<tr><th>Sample</th><th>Input read pair</th><th>Output paired reads</th><th>Status</th></tr>
</thead>
<tbody>
"""

for r in rows:
    status = r.get("status", "")
    cls = "ok" if status == "success" else "warn"
    out += (
        f"<tr><td>{html.escape(r['sample'])}</td>"
        f"<td>{html.escape(r['input_pairs'])}</td>"
        f"<td>{html.escape(r['paired_reads_output'])}</td>"
        f"<td class='{cls}'>{html.escape(status)}</td></tr>\n"
    )

out += "</tbody></table></body></html>"
open("trimming_report.html", "w").write(out)
PY
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "16 GB"
    disks: "local-disk 300 HDD"
  }

  output {
    Array[File] trimmed_reads = glob("trimmed/*_paired.fastq.gz")
    Array[File] trimming_logs = glob("logs/*.log")
    File trimming_summary = "trimming_summary.tsv"
    File trimming_report = "trimming_report.html"
  }
}
task FASTQC {
  input {
    # Biocontainers FastQC image has a stable fastqc executable, but may not include python.
    # This task therefore generates its summary using bash only.
    String docker_image = "quay.io/biocontainers/fastqc:0.11.9--0"
    Array[File]+ input_reads
    Int cpu = 8
  }

  command <<<
    set -uo pipefail
    mkdir -p fastqc logs

    echo "Running FastQC on localized input reads..." > logs/fastqc.command.log

    if command -v fastqc >/dev/null 2>&1; then
      FASTQC_BIN="$(command -v fastqc)"
    elif [ -x /opt/conda/bin/fastqc ]; then
      FASTQC_BIN="/opt/conda/bin/fastqc"
    elif [ -x /usr/local/bin/fastqc ]; then
      FASTQC_BIN="/usr/local/bin/fastqc"
    elif [ -x /usr/bin/fastqc ]; then
      FASTQC_BIN="/usr/bin/fastqc"
    else
      FASTQC_BIN=""
      echo "WARNING: FastQC executable was not found. Creating fallback QC summary and continuing." >> logs/fastqc.command.log
    fi

    if [ -n "$FASTQC_BIN" ]; then
      echo "Using FastQC: ${FASTQC_BIN}" >> logs/fastqc.command.log
      if "$FASTQC_BIN" -t ~{cpu} -o fastqc ~{sep=' ' input_reads} >> logs/fastqc.command.log 2>&1; then
        fastqc_status="success"
      else
        fastqc_status="fastqc_failed"
        echo "WARNING: FastQC failed. A summary report will still be generated." >> logs/fastqc.command.log
      fi
    else
      fastqc_status="fastqc_not_found"
    fi

    echo -e "sample\tfastqc_html\tfastqc_zip\tstatus" > fastqc_summary.tsv

    shopt -s nullglob
    html_files=(fastqc/*_fastqc.html)

    if [ ${#html_files[@]} -eq 0 ]; then
      echo -e "NO_FASTQC_OUTPUT\tNA\tNA\t${fastqc_status}" >> fastqc_summary.tsv
    else
      for html_file in "${html_files[@]}"; do
        base=$(basename "$html_file" _fastqc.html)
        zip_file="fastqc/${base}_fastqc.zip"
        if [ -f "$zip_file" ]; then
          echo -e "${base}\t${html_file}\t${zip_file}\t${fastqc_status}" >> fastqc_summary.tsv
        else
          echo -e "${base}\t${html_file}\tNA\twarning_zip_missing" >> fastqc_summary.tsv
        fi
      done
    fi

    cat > fastqc_summary.html <<'HTML_HEAD'
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>FastQC summary</title>
<style>
body{font-family:Arial;margin:24px;background:#f8fafc;color:#111827}
table{border-collapse:collapse;width:100%;background:white}
th,td{border:1px solid #dbe4ee;padding:8px;vertical-align:top}
th{background:#2563eb;color:white}
.ok{color:#166534;font-weight:bold}.warn{color:#d97706;font-weight:bold}.fail{color:#b91c1c;font-weight:bold}
</style>
</head>
<body>
<h1>FastQC per-sample report summary</h1>
<p>FastQC was run as a dedicated task. If FastQC failed or was unavailable, the workflow continued and the status is recorded here.</p>
<table>
<thead><tr><th>Sample/read file</th><th>FastQC HTML</th><th>FastQC ZIP</th><th>Status</th></tr></thead>
<tbody>
HTML_HEAD

    tail -n +2 fastqc_summary.tsv | while IFS=$'\t' read -r sample html zip status; do
      cls="warn"
      if [ "$status" = "success" ]; then cls="ok"; fi
      if echo "$status" | grep -qi "failed\|not_found"; then cls="fail"; fi
      printf '<tr><td>%s</td><td>%s</td><td>%s</td><td class="%s">%s</td></tr>\n' "$sample" "$html" "$zip" "$cls" "$status" >> fastqc_summary.html
    done

    cat >> fastqc_summary.html <<'HTML_TAIL'
</tbody>
</table>
</body>
</html>
HTML_TAIL

    # Always exit successfully so QC never blocks downstream TB-Profiler/phylogenomics.
    exit 0
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    continueOnReturnCode: [0]
  }

  output {
    Array[File] fastqc_reports = glob("fastqc/*_fastqc.html")
    Array[File] fastqc_zips = glob("fastqc/*_fastqc.zip")
    File fastqc_summary_html = "fastqc_summary.html"
    File fastqc_summary_tsv = "fastqc_summary.tsv"
    File fastqc_log = "logs/fastqc.command.log"
  }
}

task MULTIQC {
  input {
    String docker_image = "multiqc/multiqc:v1.25"
    Array[File] fastqc_reports
    Array[File] fastqc_zips
  }

  command <<<
    set -uo pipefail
    mkdir -p fastqc_input multiqc logs

    reports=(~{sep=' ' fastqc_reports})
    zips=(~{sep=' ' fastqc_zips})

    for f in "${reports[@]}"; do
      if [ -n "$f" ] && [ -f "$f" ]; then
        cp "$f" fastqc_input/
      fi
    done

    for f in "${zips[@]}"; do
      if [ -n "$f" ] && [ -f "$f" ]; then
        cp "$f" fastqc_input/
      fi
    done

    if command -v multiqc >/dev/null 2>&1; then
      MULTIQC_BIN="$(command -v multiqc)"
    elif [ -x /opt/conda/bin/multiqc ]; then
      MULTIQC_BIN="/opt/conda/bin/multiqc"
    elif [ -x /usr/local/bin/multiqc ]; then
      MULTIQC_BIN="/usr/local/bin/multiqc"
    else
      echo "ERROR: MultiQC executable was not found in PATH or common locations." >&2
      MULTIQC_BIN=""
    fi

    echo "Using MultiQC: ${MULTIQC_BIN}" > logs/multiqc.command.log
    echo "FastQC input files:" >> logs/multiqc.command.log
    find fastqc_input -maxdepth 1 -type f | sort >> logs/multiqc.command.log || true

    if [ -n "$MULTIQC_BIN" ] && \
       ( ls fastqc_input/*_fastqc.zip >/dev/null 2>&1 || ls fastqc_input/*_fastqc.html >/dev/null 2>&1 ); then

      if "${MULTIQC_BIN}" fastqc_input -o multiqc --force >> logs/multiqc.command.log 2>&1; then
        echo "MultiQC completed successfully." >> logs/multiqc.command.log
      else
        echo "WARNING: MultiQC failed. A fallback report will be generated." >> logs/multiqc.command.log
      fi
    else
      echo "WARNING: No FastQC outputs were available for MultiQC, or MultiQC was not found." >> logs/multiqc.command.log
    fi

    if [ ! -f multiqc/multiqc_report.html ]; then
      cat > multiqc/multiqc_report.html <<'HTML'
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>QC summary</title>
<style>
body{font-family:Arial;margin:24px;background:#f8fafc;color:#111827}
.warn{background:#fef3c7;color:#92400e;padding:12px;border-radius:10px;border-left:5px solid #d97706}
</style>
</head>
<body>
<h1>QC summary</h1>
<p class="warn">MultiQC did not generate a combined report. Check the FastQC task outputs and MultiQC stderr/log files.</p>
</body>
</html>
HTML
    fi
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File multiqc_report = "multiqc/multiqc_report.html"
    Array[File] multiqc_data = glob("multiqc/*")
    File multiqc_log = "logs/multiqc.command.log"
  }
}
task SPECIES_TYPING {
  input {
    Array[File]+ input_reads
    File? qc_dependency
    String docker_image = "gmboowa/mycobacterium-kraken2-bracken:2026.05"
    Int cpu = 16
    Int memory_gb = 64
  }

  command <<<
    set -euo pipefail

    mkdir -p species_typing kraken_reports kraken_outputs

    files=(~{sep=' ' input_reads})
    n=${#files[@]}

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: input_reads must contain paired-end FASTQ files in R1/R2 order." >&2
      exit 1
    fi

    ###########################################################################
    # Species typing output
    #
    # MTBC_Supported is the explicit downstream-selection field.
    #
    # Downstream logic:
    #   - TB_PROFILER_AND_MTBC_FILTER should select samples for Snippy/core-SNP/
    #     IQ-TREE only when MTBC_Supported == YES.
    #   - TB-Profiler species/lineage/resistance are annotations only and should
    #     not determine phylogeny selection.
    ###########################################################################

    echo -e "Sample_ID\tSpecies_Identified\tEvidence\tMTBC_Supported\tMTBC_Reads\tMTBC_Percent\tSelection_Basis" > species_typing/species_typing.tsv

    for ((i=0; i<n; i+=2)); do
      r1="${files[$i]}"
      r2="${files[$((i+1))]}"

      base="$(basename "$r1")"
      sample_id="$base"
      sample_id="${sample_id%_R1_paired.fastq.gz}"
      sample_id="${sample_id%_R2_paired.fastq.gz}"
      sample_id="${sample_id%_R1_paired.fq.gz}"
      sample_id="${sample_id%_R2_paired.fq.gz}"
      sample_id="${sample_id%_R1.fastq.gz}"
      sample_id="${sample_id%_R2.fastq.gz}"
      sample_id="${sample_id%_R1.fq.gz}"
      sample_id="${sample_id%_R2.fq.gz}"
      sample_id="${sample_id%_1.fastq.gz}"
      sample_id="${sample_id%_2.fastq.gz}"
      sample_id="${sample_id%_1.fq.gz}"
      sample_id="${sample_id%_2.fq.gz}"
      sample_id="${sample_id%.fastq.gz}"
      sample_id="${sample_id%.fq.gz}"

      report="kraken_reports/${sample_id}.kraken2.report.txt"
      output="kraken_outputs/${sample_id}.kraken2.output.txt"

      kraken2 \
        --db /opt/kraken2_db/mycobacterium \
        --threads ~{cpu} \
        --paired \
        --report "$report" \
        --output "$output" \
        "$r1" "$r2"

      #########################################################################
      # Kraken2 report columns:
      #   1 = percent
      #   2 = reads assigned to clade
      #   3 = reads assigned directly to taxon
      #   4 = rank code
      #   5 = NCBI taxid
      #   6+ = taxon name
      #
      # TaxID 77643 is the Mycobacterium tuberculosis complex node.
      #########################################################################

      top_species_line="$(awk '$4=="S"' "$report" | sort -k2,2nr | head -1 || true)"
      mtbc_line="$(awk '$4=="G1" && $5=="77643"' "$report" | head -1 || true)"

      if [ -z "$top_species_line" ]; then
        final_call="No species-level Mycobacterium call"
        top_percent="0.00"
        top_reads="0"
        top_taxid="NA"
        evidence="No species-level Kraken2 assignment detected"
      else
        top_percent="$(echo "$top_species_line" | awk '{print $1}')"
        top_reads="$(echo "$top_species_line" | awk '{print $2}')"
        top_taxid="$(echo "$top_species_line" | awk '{print $5}')"
        top_species="$(echo "$top_species_line" | awk '{for(i=6;i<=NF;i++) printf $i (i<NF ? " " : "")}')"
        final_call="$top_species"
      fi

      if [ -n "$mtbc_line" ]; then
        mtbc_percent="$(echo "$mtbc_line" | awk '{print $1}')"
        mtbc_reads="$(echo "$mtbc_line" | awk '{print $2}')"
      else
        mtbc_percent="0.00"
        mtbc_reads="0"
      fi

      #########################################################################
      # Conservative top-species-first MTBC routing rule.
      #
      # The top species-level Kraken2 call is treated as the primary routing
      # evidence. A clear non-MTBC Mycobacterium top call, such as
      # Mycobacterium avium, is NOT routed into the MTBC workflow even when the
      # MTBC-complex node has background/cross-assigned reads.
      #
      # Routing order:
      #   1. If the top species-level call is a recognized MTBC member, route as
      #      MTBC.
      #   2. If the top species-level call is a clear non-MTBC Mycobacterium,
      #      exclude from MTBC-specific downstream analyses.
      #   3. Only when the species-level call is unresolved, allow strong MTBC
      #      complex-node support (>=10% and >=1,000 reads) to route as MTBC.
      #########################################################################

      final_call_lc="$(echo "$final_call" | tr '[:upper:]' '[:lower:]')"

      mtbc_percent_float="$(python3 - <<PYMTBC
try:
    print(float("${mtbc_percent}"))
except Exception:
    print(0.0)
PYMTBC
)"

      mtbc_reads_int="$(python3 - <<PYREADS
try:
    print(int(float("${mtbc_reads}")))
except Exception:
    print(0)
PYREADS
)"

      top_call_is_mtbc="NO"
      if echo "$final_call_lc" | grep -Eq "mycobacterium tuberculosis complex|tuberculosis complex|(^|[^a-z])mtbc([^a-z]|$)|mycobacterium tuberculosis|m\. tuberculosis|mycobacterium africanum|mycobacterium bovis|mycobacterium caprae|mycobacterium canettii|mycobacterium microti|mycobacterium pinnipedii|mycobacterium orygis|mycobacterium mungi|mycobacterium suricattae|dassie bacillus"; then
        top_call_is_mtbc="YES"
      fi

      top_call_is_clear_non_mtbc="NO"
      if [ "$top_call_is_mtbc" = "NO" ] && echo "$final_call_lc" | grep -q "mycobacterium"; then
        if ! echo "$final_call_lc" | grep -Eq "no species-level|unknown|not resolved|unresolved|^na$|^n/a$"; then
          top_call_is_clear_non_mtbc="YES"
        fi
      fi

      strong_mtbc_node_support="$(python3 - <<PYSTRONG
pct = float("${mtbc_percent_float}")
reads = int("${mtbc_reads_int}")
print("YES" if (pct >= 10.0 and reads >= 1000) else "NO")
PYSTRONG
)"

      if [ "$top_call_is_mtbc" = "YES" ]; then
        mtbc_supported="YES"
        selection_basis="Selected for MTBC workflow because the top Kraken2 species-level call is a recognized MTBC member"
      elif [ "$top_call_is_clear_non_mtbc" = "YES" ]; then
        mtbc_supported="NO"
        selection_basis="Not selected for MTBC workflow because the top Kraken2 species-level call is a clear non-MTBC Mycobacterium; MTBC-complex reads are treated as background/cross-assignment"
      elif [ "$strong_mtbc_node_support" = "YES" ]; then
        mtbc_supported="YES"
        selection_basis="Selected for MTBC workflow because the species-level call was unresolved and Kraken2 detected strong MTBC-complex support: ${mtbc_reads} reads; ${mtbc_percent}%"
      else
        mtbc_supported="NO"
        selection_basis="Not selected for MTBC workflow because Kraken2/Bracken species typing did not provide MTBC support; low-level MTBC reads alone are treated as background/cross-assignment"
      fi

      if [ -z "$top_species_line" ]; then
        evidence="No species-level Kraken2 assignment detected; MTBC support: ${mtbc_reads} reads; ${mtbc_percent}%"
      else
        evidence="Top species-level assignment: ${final_call} (${top_reads} reads; ${top_percent}%); MTBC support: ${mtbc_reads} reads; ${mtbc_percent}%"
      fi

      echo -e "${sample_id}\t${final_call}\t${evidence}\t${mtbc_supported}\t${mtbc_reads}\t${mtbc_percent}\t${selection_basis}" >> species_typing/species_typing.tsv
    done

    python3 <<'PY'
import csv
from html import escape

rows = []
with open("species_typing/species_typing.tsv", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        rows.append(row)

def species_badge(label):
    label = label or ""
    text = label.strip().lower()

    if (
        text == "mycobacterium tuberculosis"
        or text == "mycobacterium tuberculosis complex"
        or "mycobacterium tuberculosis" in text
        or "tuberculosis complex" in text
        or text == "mtbc"
    ):
        return (
            '<span style="display:inline-block;'
            'padding:4px 8px;'
            'border-radius:999px;'
            'background:#6F42C1;'
            'color:white;'
            'font-size:12px;'
            'font-weight:700;">'
            f'{escape(label)}</span>'
        )

    if text == "no species-level mycobacterium call":
        return (
            '<span style="display:inline-block;'
            'padding:4px 8px;'
            'border-radius:999px;'
            'background:#9CA3AF;'
            'color:white;'
            'font-size:12px;'
            'font-weight:700;">'
            f'{escape(label)}</span>'
        )

    return (
        '<span style="display:inline-block;'
        'padding:4px 8px;'
        'border-radius:999px;'
        'background:#28A745;'
        'color:white;'
        'font-size:12px;'
        'font-weight:700;">'
        f'{escape(label)}</span>'
    )

def mtbc_badge(value):
    value = (value or "").strip().upper()

    if value == "YES":
        return (
            '<span style="display:inline-block;'
            'padding:4px 8px;'
            'border-radius:999px;'
            'background:#28A745;'
            'color:white;'
            'font-size:12px;'
            'font-weight:700;">YES</span>'
        )

    return (
        '<span style="display:inline-block;'
        'padding:4px 8px;'
        'border-radius:999px;'
        'background:#DC2626;'
        'color:white;'
        'font-size:12px;'
        'font-weight:700;">NO</span>'
    )

html = []
html.append('<section class="card">')
html.append('<h2>2. Species Typing using Kraken2 + Bracken</h2>')
html.append('<p>Species typing was performed using Kraken2 against a custom Mycobacterium-only database embedded in the Docker image <code>gmboowa/mycobacterium-kraken2-bracken:2026.05</code>. The table reports the single most probable species-level call for each sample based on the highest species-level Kraken2 assignment and supporting taxonomic evidence.</p>')
html.append('<p><strong>Phylogeny selection rule:</strong> samples are selected for downstream Snippy/core-SNP/IQ-TREE analysis only when Kraken2/Bracken supports MTBC. TB-Profiler species, lineage, and resistance are reported later as annotations but do not determine phylogeny selection.</p>')
html.append('<div class="table-wrap">')
html.append('<table>')
html.append('<thead><tr>')
html.append('<th>Sample ID</th>')
html.append('<th>Species Identified</th>')
html.append('<th>Evidence Supporting Call</th>')
html.append('<th>MTBC Supported</th>')
html.append('<th>MTBC Reads</th>')
html.append('<th>MTBC Percent</th>')
html.append('<th>Selection Basis</th>')
html.append('</tr></thead>')
html.append('<tbody>')

if rows:
    for row in rows:
        html.append(
            "<tr>"
            f"<td>{escape(row.get('Sample_ID',''))}</td>"
            f"<td>{species_badge(row.get('Species_Identified',''))}</td>"
            f"<td>{escape(row.get('Evidence',''))}</td>"
            f"<td>{mtbc_badge(row.get('MTBC_Supported',''))}</td>"
            f"<td>{escape(row.get('MTBC_Reads',''))}</td>"
            f"<td>{escape(row.get('MTBC_Percent',''))}</td>"
            f"<td>{escape(row.get('Selection_Basis',''))}</td>"
            "</tr>"
        )
else:
    html.append('<tr><td colspan="7">No species typing results were generated.</td></tr>')

html.append('</tbody></table></div></section>')

with open("species_typing/species_typing.html", "w") as out:
    out.write("\n".join(html))
PY
  >>>

  output {
    Array[File] kraken_reports = glob("kraken_reports/*.report.txt")
    Array[File] kraken_outputs = glob("kraken_outputs/*.output.txt")
    File species_typing_tsv = "species_typing/species_typing.tsv"
    File species_typing_html = "species_typing/species_typing.html"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 200 HDD"
    timeout: "24 hours"
  }
}
task MYCOBACTERIA_SAMPLE_ROUTER {
  input {
    Array[File]+ input_reads
    File species_typing_tsv
    String docker_image = "python:3.11-slim"
    Int cpu = 2
    Int memory_gb = 8
  }

  command <<<
    set -euo pipefail

    mkdir -p routed_reads/mtbc routing ntm_species

    python3 - <<'PY'
import csv
import html
import re
import shutil
from pathlib import Path

species_tsv = Path("~{species_typing_tsv}")
input_files = [Path(x) for x in """~{sep='\n' input_reads}""".splitlines() if x.strip()]

outdir = Path("routing")
outdir.mkdir(exist_ok=True)
ntm_dir = Path("ntm_species")
ntm_dir.mkdir(exist_ok=True)
mtbc_dir = Path("routed_reads/mtbc")
mtbc_dir.mkdir(parents=True, exist_ok=True)

def sample_from_name(path):
    name = Path(path).name
    patterns = [
        r"_R1_paired\.fastq\.gz$", r"_R2_paired\.fastq\.gz$",
        r"_R1_paired\.fq\.gz$", r"_R2_paired\.fq\.gz$",
        r"_R1\.fastq\.gz$", r"_R2\.fastq\.gz$",
        r"_R1\.fq\.gz$", r"_R2\.fq\.gz$",
        r"_1\.fastq\.gz$", r"_2\.fastq\.gz$",
        r"_1\.fq\.gz$", r"_2\.fq\.gz$",
        r"\.fastq\.gz$", r"\.fq\.gz$", r"\.fastq$", r"\.fq$"
    ]
    for pat in patterns:
        name = re.sub(pat, "", name)
    return name

pairs = {}
for i in range(0, len(input_files), 2):
    r1 = input_files[i]
    r2 = input_files[i + 1] if i + 1 < len(input_files) else None
    sample = sample_from_name(r1)
    if r2 is not None:
        pairs[sample] = (r1, r2)

rows = []
if species_tsv.exists() and species_tsv.stat().st_size > 0:
    with species_tsv.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="	"))

def clean(v):
    return str(v or "").strip()

def is_mtbc_supported(row):
    return clean(row.get("MTBC_Supported") or row.get("mtbc_supported")).upper() == "YES"

def species_is_mtbc(species):
    text = clean(species).lower()
    mtbc_terms = [
        "mycobacterium tuberculosis complex",
        "tuberculosis complex",
        "mtbc",
        "mycobacterium tuberculosis",
        "m. tuberculosis",
        "mycobacterium africanum",
        "mycobacterium bovis",
        "mycobacterium caprae",
        "mycobacterium canettii",
        "mycobacterium microti",
        "mycobacterium pinnipedii",
        "mycobacterium orygis",
        "mycobacterium mungi",
        "mycobacterium suricattae",
        "dassie bacillus"
    ]
    return any(term in text for term in mtbc_terms)

def species_is_unresolved(species):
    text = clean(species).lower()
    return (not text) or text in ["no species-level mycobacterium call", "no species-level call", "unknown", "na", "n/a", "not resolved", "unresolved"]

def species_is_clear_non_mtbc_mycobacterium(species):
    text = clean(species).lower()
    return ("mycobacterium" in text) and (not species_is_mtbc(text)) and (not species_is_unresolved(text))

def classify(row):
    species = clean(row.get("Species_Identified") or row.get("Species identified") or row.get("species") or row.get("Species"))
    evidence = clean(row.get("Evidence") or row.get("evidence"))

    # Top-species-first routing: a clear non-MTBC Mycobacterium top call
    # overrides any background MTBC-complex read signal.
    if species_is_mtbc(species):
        return "MTBC"
    if species_is_clear_non_mtbc_mycobacterium(species):
        return "NTM"
    if is_mtbc_supported(row):
        return "MTBC"
    if species_is_unresolved(species):
        return "Unresolved Mycobacterium"
    if "mycobacterium" in f"{species} {evidence}".lower():
        return "NTM"
    return "Non-Mycobacterium or low-confidence"

def decision_for(classification):
    if classification == "MTBC":
        return "Retained for TB-Profiler, MTBC AMR interpretation, and MTBC-specific phylogenomics"
    if classification == "NTM":
        return "Reported as non-MTBC Mycobacterium and excluded from TB-Profiler and MTBC-specific phylogenomics"
    if classification == "Unresolved Mycobacterium":
        return "Reported as unresolved Mycobacterium and excluded from MTBC-specific downstream analyses"
    return "Reported and excluded from MTBC-specific downstream analyses"

routing_rows = []
mtbc_count = ntm_count = unresolved_count = non_myco_count = 0

for row in rows:
    sample = clean(row.get("Sample_ID") or row.get("sample") or row.get("sample_id") or row.get("Sample ID"))
    species = clean(row.get("Species_Identified") or row.get("Species identified") or row.get("species") or row.get("Species")) or "Not resolved"
    evidence = clean(row.get("Evidence") or row.get("evidence")) or "No supporting evidence reported"
    mtbc_supported = clean(row.get("MTBC_Supported") or row.get("mtbc_supported")) or "NO"
    mtbc_reads = clean(row.get("MTBC_Reads") or row.get("mtbc_reads")) or "0"
    mtbc_percent = clean(row.get("MTBC_Percent") or row.get("mtbc_percent")) or "0.00"
    selection_basis = clean(row.get("Selection_Basis") or row.get("selection_basis"))

    classification = classify(row)
    decision = decision_for(classification)

    if classification == "MTBC":
        mtbc_count += 1
        if sample in pairs:
            r1, r2 = pairs[sample]
            shutil.copyfile(r1, mtbc_dir / f"{sample}_R1.fastq.gz")
            shutil.copyfile(r2, mtbc_dir / f"{sample}_R2.fastq.gz")
    elif classification == "NTM":
        ntm_count += 1
    elif classification == "Unresolved Mycobacterium":
        unresolved_count += 1
    else:
        non_myco_count += 1

    routing_rows.append({
        "Sample_ID": sample,
        "Most_Probable_Mycobacterium_Species": species,
        "Mycobacteria_Category": classification,
        "MTBC_Supported": mtbc_supported,
        "MTBC_Reads": mtbc_reads,
        "MTBC_Percent": mtbc_percent,
        "Evidence": evidence,
        "Workflow_Decision": decision,
        "Selection_Basis": selection_basis
    })

fields = [
    "Sample_ID", "Most_Probable_Mycobacterium_Species", "Mycobacteria_Category",
    "MTBC_Supported", "MTBC_Reads", "MTBC_Percent", "Evidence",
    "Workflow_Decision", "Selection_Basis"
]

def write_rows(path, selected):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="	")
        writer.writeheader()
        writer.writerows(selected)

write_rows(outdir / "mycobacteria_routing_summary.tsv", routing_rows)
write_rows(ntm_dir / "ntm_species_summary.tsv", [r for r in routing_rows if r["Mycobacteria_Category"] == "NTM"])
write_rows(outdir / "mtbc_species_summary.tsv", [r for r in routing_rows if r["Mycobacteria_Category"] == "MTBC"])
write_rows(outdir / "non_mtbc_species_summary.tsv", [r for r in routing_rows if r["Mycobacteria_Category"] != "MTBC"])

has_mtbc = mtbc_count > 0
has_ntm = ntm_count > 0
all_ntm = mtbc_count == 0 and ntm_count > 0

(outdir / "mtbc_sample_count.txt").write_text(str(mtbc_count) + "\n")
(outdir / "ntm_sample_count.txt").write_text(str(ntm_count) + "\n")
(outdir / "unresolved_mycobacterium_count.txt").write_text(str(unresolved_count) + "\n")
(outdir / "non_mycobacterium_count.txt").write_text(str(non_myco_count) + "\n")
(outdir / "has_mtbc.txt").write_text(("true" if has_mtbc else "false") + "\n")
(outdir / "has_ntm.txt").write_text(("true" if has_ntm else "false") + "\n")
(outdir / "all_ntm.txt").write_text(("true" if all_ntm else "false") + "\n")

if all_ntm:
    status = "All classified Mycobacterium samples were non-MTBC/NTM. TB-Profiler and MTBC-specific phylogenomics should be skipped."
elif has_mtbc and has_ntm:
    status = "Mixed MTBC and NTM dataset. MTBC-supported samples are routed forward; NTM samples are reported and excluded from MTBC-specific analyses."
elif has_mtbc:
    status = "MTBC-supported dataset. MTBC samples are routed forward for TB-Profiler and MTBC-specific analyses."
else:
    status = "No MTBC-supported samples were detected. MTBC-specific downstream analyses should be skipped."

(outdir / "routing_status.txt").write_text(status + "\n")

body = []
for r in routing_rows:
    if r["Mycobacteria_Category"] != "MTBC":
        body.append(
            "<tr>"
            f"<td>{html.escape(r['Sample_ID'])}</td>"
            f"<td><em>{html.escape(r['Most_Probable_Mycobacterium_Species'])}</em></td>"
            f"<td>{html.escape(r['Mycobacteria_Category'])}</td>"
            f"<td>{html.escape(r['MTBC_Reads'])}</td>"
            f"<td>{html.escape(r['MTBC_Percent'])}</td>"
            f"<td>{html.escape(r['Evidence'])}</td>"
            f"<td>{html.escape(r['Workflow_Decision'])}</td>"
            "</tr>"
        )
if not body:
    body.append("<tr><td colspan='7'>No non-MTBC Mycobacterium or NTM samples were detected by Kraken2/Bracken.</td></tr>")

html_text = f"""<section class="card">
<h2>Non-MTBC Mycobacteria Species Summary</h2>
<p><strong>Routing status:</strong> {html.escape(status)}</p>
<p>Kraken2/Bracken species typing was used as the first routing layer. MTBC-supported samples are retained for TB-Profiler and MTBC-specific phylogenomics. Non-MTBC Mycobacterium samples are reported here and excluded from MTBC-specific downstream analyses.</p>
<table>
<thead><tr><th>Sample ID</th><th>Most probable Mycobacterium species</th><th>Category</th><th>MTBC reads</th><th>MTBC percent</th><th>Evidence</th><th>Workflow decision</th></tr></thead>
<tbody>{''.join(body)}</tbody>
</table>
</section>
"""
(ntm_dir / "ntm_species_summary.html").write_text(html_text)
PY
  >>>

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 100 HDD"
    timeout: "24 hours"
  }

  output {
    Array[File] mtbc_reads = glob("routed_reads/mtbc/*.fastq.gz")
    File mycobacteria_routing_summary_tsv = "routing/mycobacteria_routing_summary.tsv"
    File ntm_species_summary_tsv = "ntm_species/ntm_species_summary.tsv"
    File ntm_species_summary_html = "ntm_species/ntm_species_summary.html"
    File mtbc_species_summary_tsv = "routing/mtbc_species_summary.tsv"
    File non_mtbc_species_summary_tsv = "routing/non_mtbc_species_summary.tsv"
    File routing_status_txt = "routing/routing_status.txt"
    Int mtbc_sample_count = read_int("routing/mtbc_sample_count.txt")
    Int ntm_sample_count = read_int("routing/ntm_sample_count.txt")
    Int unresolved_mycobacterium_count = read_int("routing/unresolved_mycobacterium_count.txt")
    Int non_mycobacterium_count = read_int("routing/non_mycobacterium_count.txt")
    Boolean has_mtbc = read_boolean("routing/has_mtbc.txt")
    Boolean has_ntm = read_boolean("routing/has_ntm.txt")
    Boolean all_ntm = read_boolean("routing/all_ntm.txt")
  }
}

task TB_PROFILER_AND_MTBC_FILTER {
  input {
    Array[File]+ input_reads
    File? qc_dependency
    File? species_typing_tsv
    String docker_image = "staphb/tbprofiler:6.6.6"
    Int cpu = 16
    Int memory_gb = 64
  }

  command <<<
    set -uo pipefail
    mkdir -p tbprofiler_results mtbc_reads logs mutation_evidence

    # Cromwell requires declared output files to exist at task completion.
    # Create safe placeholders up front; successful report-rendering steps overwrite them later.
    cat > tbprofiler_combined_report.html <<'EOF_TBPROF_PLACEHOLDER'
<!doctype html><html><head><meta charset="utf-8"><title>TB-Profiler report</title></head><body><h1>TB-Profiler report</h1><p>TB-Profiler report generation started, but no final rendered table was produced.</p></body></html>
EOF_TBPROF_PLACEHOLDER

    cat > mutation_evidence/tbprofiler_mutation_evidence.html <<'EOF_MUT_PLACEHOLDER'
<!doctype html><html><head><meta charset="utf-8"><title>TB-Profiler mutation evidence</title></head><body><h1>Resistance Mutation Evidence Summary</h1><p>No rendered mutation-evidence table was produced.</p></body></html>
EOF_MUT_PLACEHOLDER

    species_tsv="~{if defined(species_typing_tsv) then species_typing_tsv else ""}"

    if command -v tb-profiler >/dev/null 2>&1; then
      TBPROFILER_BIN="$(command -v tb-profiler)"
    elif command -v tb_profiler >/dev/null 2>&1; then
      TBPROFILER_BIN="$(command -v tb_profiler)"
    elif [ -x /opt/conda/bin/tb-profiler ]; then
      TBPROFILER_BIN="/opt/conda/bin/tb-profiler"
    else
      echo "ERROR: TB-Profiler executable not found." >&2
      exit 127
    fi

    echo "Using TB-Profiler: ${TBPROFILER_BIN}" > logs/tbprofiler.command.log

    files=(~{sep=' ' input_reads})
    n=${#files[@]}

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: input_reads must contain paired reads in R1/R2 order." >&2
      exit 1
    fi

    echo -e "sample\tspecies\tmain_lineage\tsub_lineage\tdr_type\tresistant_drugs\tresistance_mutations\tkey_mutations\tjson_file\tmtbc_selected\tmtbc_selection_reason\tstatus" > tbprofiler_summary.tsv
    echo -e "sample_id\tresistance_profile\tresistant_drugs\tresistance_mutations\tkey_mutations\tstatus" > resistance_profile_summary.tsv
    echo -e "sample\tdrug\tgene\tmutation\tchange\tconfidence\tevidence\tsource_json" > mutation_evidence/tbprofiler_mutation_evidence.tsv
    : > mtbc_samples.txt

    for ((i=0; i<n; i+=2)); do
      R1="${files[$i]}"
      R2="${files[$((i+1))]}"
      sample=$(basename "$R1")
      sample=$(echo "$sample" | sed -E 's/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//' | sed -E 's/(_R?1_paired|_R?1|_1_paired|_1|\.R?1|\.1)(_|$).*//')

      echo "Running TB-Profiler for ${sample}" >> logs/tbprofiler.command.log

      if "${TBPROFILER_BIN}" profile \
        -1 "$R1" \
        -2 "$R2" \
        -p "$sample" \
        --dir tbprofiler_results \
        --threads ~{cpu} >> "logs/${sample}.tbprofiler.log" 2>&1; then
        status="success"
      else
        status="tbprofiler_failed"
        echo "WARNING: TB-Profiler failed for ${sample}; workflow will continue." >> logs/tbprofiler.command.log
      fi

      json=$(find tbprofiler_results -name "${sample}*.json" | head -n 1 || true)

      if [ -z "$json" ]; then
        mkdir -p tbprofiler_results/results
        json="tbprofiler_results/results/${sample}.results.json"
        echo "{\"id\":\"${sample}\",\"sample\":\"${sample}\",\"error\":\"TB-Profiler did not generate JSON\"}" > "$json"
      fi

      python3 - "$json" "$sample" "$R1" "$R2" "$status" "$species_tsv" <<'PYTB'
import json, sys, os, shutil, re, csv

json_file, sample, r1, r2, status, species_tsv = sys.argv[1:7]

try:
    data = json.load(open(json_file))
except Exception as e:
    data = {"error": str(e)}

def get_path(obj, paths):
    for path in paths:
        cur = obj
        ok = True
        for part in path.split("."):
            if isinstance(cur, dict) and part in cur:
                cur = cur[part]
            else:
                ok = False
                break
        if ok and cur not in (None, "", [], {}):
            return cur
    return ""

def clean_value(v):
    if v in (None, "", [], {}):
        return ""
    if isinstance(v, str):
        return v.strip()
    if isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, list):
        out = []
        for x in v:
            y = clean_value(x)
            if y:
                out.append(y)
        return "; ".join(out)
    if isinstance(v, dict):
        for key in [
            "name", "drug", "gene", "change", "mutation",
            "original_mutation", "confidence", "source", "evidence"
        ]:
            if key in v and v[key] not in (None, "", [], {}):
                return clean_value(v[key])
        return ""
    return str(v)

def uniq_list(xs):
    seen = []
    seen_lower = set()
    for x in xs:
        x = str(x or "").strip()
        xl = x.lower()
        if not x:
            continue
        if xl in ["none", "none reported", "not reported", "na", "n/a", "unknown"]:
            continue
        if xl not in seen_lower:
            seen.append(x)
            seen_lower.add(xl)
    return seen

def uniq(xs):
    return ", ".join(uniq_list(xs))

def safe_tsv(x):
    return str(x or "").replace("\t", " ").replace("\n", " ").replace("\r", " ").strip()

def normalize_sample_id(name):
    s = str(name or "").strip()
    s = s.split("/")[-1]
    s = s.split("\\")[-1]
    s = re.sub(r"(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq|\.gz)$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"(_R?1_paired|_R?2_paired|_R?1|_R?2|_1_paired|_2_paired|_1|_2|\.R?1|\.R?2|\.1|\.2)$", "", s)
    s = re.sub(r"^\d{3,}[_-](?=(?:[DES]RR|SAM[END]?|ERS|SRS|DRS|SRX|ERX|DRX|[A-Za-z]))", "", s)
    return s.strip()

def normalize_drug_name(x):
    s = str(x or "").strip().lower()
    s = re.sub(r"[_\-]+", " ", s)
    s = re.sub(r"\s+", " ", s)

    aliases = {
        "inh": "isoniazid",
        "isoniazid": "isoniazid",
        "h": "isoniazid",

        "rif": "rifampicin",
        "rmp": "rifampicin",
        "rifampin": "rifampicin",
        "rifampicin": "rifampicin",
        "r": "rifampicin",

        "pza": "pyrazinamide",
        "pyrazinamide": "pyrazinamide",
        "z": "pyrazinamide",

        "emb": "ethambutol",
        "ethambutol": "ethambutol",
        "e": "ethambutol",

        "sm": "streptomycin",
        "str": "streptomycin",
        "streptomycin": "streptomycin",
        "s": "streptomycin",

        "levo": "levofloxacin",
        "levofloxacin": "levofloxacin",
        "lfx": "levofloxacin",

        "moxi": "moxifloxacin",
        "moxifloxacin": "moxifloxacin",
        "mfx": "moxifloxacin",

        "ofx": "ofloxacin",
        "ofloxacin": "ofloxacin",

        "gatifloxacin": "gatifloxacin",
        "gfx": "gatifloxacin",

        "ciprofloxacin": "ciprofloxacin",
        "cfx": "ciprofloxacin",

        "amikacin": "amikacin",
        "amk": "amikacin",

        "kanamycin": "kanamycin",
        "kan": "kanamycin",

        "capreomycin": "capreomycin",
        "cap": "capreomycin",

        "bedaquiline": "bedaquiline",
        "bdq": "bedaquiline",

        "linezolid": "linezolid",
        "lzd": "linezolid",

        "clofazimine": "clofazimine",
        "cfz": "clofazimine",

        "ethionamide": "ethionamide",
        "eto": "ethionamide",

        "prothionamide": "prothionamide",
        "pto": "prothionamide",

        "cycloserine": "cycloserine",
        "cs": "cycloserine",

        "para aminosalicylic acid": "para-aminosalicylic acid",
        "para-aminosalicylic acid": "para-aminosalicylic acid",
        "pas": "para-aminosalicylic acid",

        "delamanid": "delamanid",
        "dlm": "delamanid",

        "pretomanid": "pretomanid",
        "pa": "pretomanid",
    }

    return aliases.get(s, s)

def split_drug_tokens(value):
    if not value:
        return []

    if isinstance(value, list):
        values = value
    else:
        values = [value]

    out = []

    for v in values:
        txt = str(v or "")
        txt = txt.replace(";", ",")
        txt = txt.replace("|", ",")
        txt = txt.replace("/", ",")
        txt = re.sub(r"\band\b", ",", txt, flags=re.IGNORECASE)

        for part in txt.split(","):
            p = normalize_drug_name(part)
            if p and p not in [
                "none", "none reported", "not reported",
                "susceptible", "sensitive", "na", "n/a", "unknown"
            ]:
                out.append(p)

    return out

def extract_drug_like_terms(text):
    known = [
        "isoniazid", "rifampicin", "rifampin", "pyrazinamide", "ethambutol",
        "streptomycin", "levofloxacin", "moxifloxacin", "ofloxacin",
        "gatifloxacin", "ciprofloxacin", "amikacin", "kanamycin",
        "capreomycin", "bedaquiline", "linezolid", "clofazimine",
        "ethionamide", "prothionamide", "cycloserine",
        "para-aminosalicylic acid", "delamanid", "pretomanid"
    ]

    text_lower = str(text or "").lower()
    hits = []

    for d in known:
        if d in text_lower and d not in hits:
            hits.append(d)

    return hits

def is_negative_or_non_resistance_text(value):
    s = str(value or "").lower()
    negative_terms = [
        "susceptible",
        "sensitive",
        "not resistant",
        "no resistance",
        "benign",
        "neutral",
        "synonymous",
        "phylogenetic",
        "lineage marker",
        "lineage_marker",
        "low confidence",
        "low_confidence",
        "uncertain significance",
        "unknown significance",
        "not associated with resistance"
    ]
    return any(t in s for t in negative_terms)

def is_positive_resistance_text(value):
    s = str(value or "").lower()

    if not s:
        return False

    if is_negative_or_non_resistance_text(s):
        return False

    positive_terms = [
        "resistant",
        "resistance",
        "assoc w r",
        "associated with resistance",
        "confers resistance",
        "predicted resistant",
        "high confidence",
        "moderate confidence",
        '"r"',
        "'r'"
    ]

    return any(t in s for t in positive_terms)

def should_count_variant_for_who_classification(block_name, item, drug_values, confidence, evidence, source):
    """
    Count only canonical TB-Profiler resistance calls for classification.

    Important:
      - dr_variants and resistance_variants are treated as canonical resistance evidence.
      - other_variants, candidate variants, lineage markers, synonymous/benign calls,
        and low-confidence calls are not allowed to drive WHO classification.
      - Mutation evidence may still be displayed in the evidence table, but only this
        function determines which drugs enter WHO classification.
    """
    if not drug_values:
        return False

    block_lower = str(block_name or "").lower()

    if block_lower in [
        "other_variants",
        "dr_variants_candidate",
        "candidate_variants",
        "variants",
        "mutations"
    ]:
        return False

    blob = " ".join([
        json.dumps(item, sort_keys=True),
        str(confidence or ""),
        str(evidence or ""),
        str(source or "")
    ])

    if is_negative_or_non_resistance_text(blob):
        return False

    if block_lower in ["dr_variants", "resistance_variants"]:
        return True

    return is_positive_resistance_text(blob)

def classify_who_2021_tb_resistance(resistant_drugs, tbprofiler_drtype):
    """
    Exact WHO 2021+ hierarchy implemented for rMAP-TB.

    The decision order is deliberately strict:

      1. XDR-TB:
         MDR/RR-TB plus resistance to any fluoroquinolone plus resistance to
         bedaquiline or linezolid.

      2. Pre-XDR-TB:
         MDR/RR-TB plus resistance to any fluoroquinolone.

      3. MDR-TB:
         Resistance to at least both isoniazid and rifampicin.

      4. RR-TB:
         Rifampicin-resistant TB without detected isoniazid resistance.

      5. Hr-TB:
         Isoniazid-resistant and rifampicin-susceptible TB.

      6. Monoresistance:
         Resistance to exactly one anti-TB drug, excluding cases already classified
         as RR-TB, MDR-TB, Pre-XDR-TB, XDR-TB, or Hr-TB.

      7. Polyresistance:
         Resistance to more than one anti-TB drug, excluding MDR/RR-TB,
         Pre-XDR-TB, and XDR-TB.

      8. No resistance detected by TB-Profiler.

    The function uses drug-level calls first. TB-Profiler drtype is used only
    when no parseable drug-level calls are available.
    """
    drugs = set()

    for d in resistant_drugs:
        for token in split_drug_tokens(d):
            drugs.add(token)

    drtype_text = str(tbprofiler_drtype or "").strip()
    drtype_lower = drtype_text.lower().replace("_", "-")

    has_inh = "isoniazid" in drugs
    has_rif = "rifampicin" in drugs

    fluoroquinolones = {
        "levofloxacin",
        "moxifloxacin",
        "ofloxacin",
        "gatifloxacin",
        "ciprofloxacin"
    }

    xdr_group_a_non_fq = {
        "bedaquiline",
        "linezolid"
    }

    has_fq = bool(drugs.intersection(fluoroquinolones))
    has_bedaquiline_or_linezolid = bool(drugs.intersection(xdr_group_a_non_fq))

    has_mdr_rr = has_rif

    # WHO 2021+ priority order.
    # This is the critical fix: XDR and Pre-XDR must be evaluated before
    # Hr-TB, monoresistance, or polyresistance.
    if has_mdr_rr and has_fq and has_bedaquiline_or_linezolid:
        return "XDR-TB"

    if has_mdr_rr and has_fq:
        return "Pre-XDR-TB"

    if has_rif and has_inh:
        return "MDR-TB"

    if has_rif and not has_inh:
        return "RR-TB"

    if has_inh and not has_rif:
        return "Hr-TB"

    if len(drugs) == 1:
        return "Monoresistance"

    if len(drugs) > 1:
        return "Polyresistance"

    # Fallback to TB-Profiler drtype only when no parseable drug-level calls exist.
    if "xdr" in drtype_lower and "pre" not in drtype_lower:
        return "XDR-TB"

    if "pre-xdr" in drtype_lower or "pre xdr" in drtype_lower or "prexdr" in drtype_lower:
        return "Pre-XDR-TB"

    if "mdr" in drtype_lower:
        return "MDR-TB"

    if re.search(r"\brr[- ]?tb\b", drtype_lower) or "rifampicin-resistant" in drtype_lower or "rifampicin resistant" in drtype_lower:
        return "RR-TB"

    if re.search(r"\bhr[- ]?tb\b", drtype_lower) or "isoniazid-resistant" in drtype_lower or "isoniazid resistant" in drtype_lower:
        return "Hr-TB"

    if "mono" in drtype_lower:
        return "Monoresistance"

    if "poly" in drtype_lower:
        return "Polyresistance"

    return "No resistance detected by TB-Profiler"

def get_kraken_support(sample, species_tsv):
    kraken_species = ""
    kraken_evidence = ""
    supports_mtbc = False

    normalized_sample = normalize_sample_id(sample)

    if species_tsv and os.path.exists(species_tsv):
        try:
            with open(species_tsv) as sfh:
                reader = csv.DictReader(sfh, delimiter="	")
                for row in reader:
                    sid = (
                        row.get("Sample_ID") or
                        row.get("sample") or
                        row.get("sample_id") or
                        row.get("Sample ID") or
                        ""
                    )

                    if normalize_sample_id(sid) == normalized_sample:
                        kraken_species = (
                            row.get("Species_Identified") or
                            row.get("species") or
                            row.get("Species identified") or
                            row.get("Species") or
                            ""
                        )
                        kraken_evidence = (
                            row.get("Evidence") or
                            row.get("evidence") or
                            row.get("Kraken2_Bracken_Evidence") or
                            ""
                        )

                        # Use the explicit router field first. Do not infer MTBC
                        # support from the word "MTBC" in the evidence sentence,
                        # because every row reports an MTBC-read count.
                        explicit = str(row.get("MTBC_Supported") or row.get("mtbc_supported") or "").strip().upper()
                        if explicit in ["YES", "TRUE", "1"]:
                            supports_mtbc = True
                        elif explicit in ["NO", "FALSE", "0"]:
                            supports_mtbc = False
                        else:
                            # Conservative fallback for legacy species TSVs.
                            # Prioritize the top species call. A clear non-MTBC
                            # Mycobacterium call should not be rescued by
                            # background MTBC-complex node reads.
                            species_lc = kraken_species.lower()
                            mtbc_terms = [
                                "mycobacterium tuberculosis complex",
                                "tuberculosis complex",
                                "mtbc",
                                "mycobacterium tuberculosis",
                                "m. tuberculosis",
                                "mycobacterium africanum",
                                "mycobacterium bovis",
                                "mycobacterium caprae",
                                "mycobacterium canettii",
                                "mycobacterium microti",
                                "mycobacterium pinnipedii",
                                "mycobacterium orygis",
                                "mycobacterium mungi",
                                "mycobacterium suricattae",
                                "dassie bacillus"
                            ]
                            top_call_is_mtbc = any(t in species_lc for t in mtbc_terms)
                            clear_non_mtbc = (
                                "mycobacterium" in species_lc
                                and not top_call_is_mtbc
                                and not any(t in species_lc for t in [
                                    "no species-level", "unknown", "not resolved", "unresolved"
                                ])
                            )
                            try:
                                mtbc_pct = float(row.get("MTBC_Percent") or row.get("mtbc_percent") or 0)
                            except Exception:
                                mtbc_pct = 0.0
                            try:
                                mtbc_reads = int(float(row.get("MTBC_Reads") or row.get("mtbc_reads") or 0))
                            except Exception:
                                mtbc_reads = 0
                            supports_mtbc = bool(top_call_is_mtbc or ((not clear_non_mtbc) and mtbc_pct >= 10.0 and mtbc_reads >= 1000))
                        break
        except Exception:
            kraken_species = ""
            kraken_evidence = ""
            supports_mtbc = False

    return kraken_species, kraken_evidence, supports_mtbc

species = clean_value(get_path(data, [
    "species", "main_species", "taxon", "organism",
    "phylogeny.species", "lineage.species"
]))

main_lineage = clean_value(get_path(data, [
    "main_lineage", "main_lin", "lineage", "lin", "phylogeny.lineage"
]))

sub_lineage = clean_value(get_path(data, [
    "sub_lineage", "sublineage", "sublin", "sub_lin", "phylogeny.sublineage"
]))

dr_type_raw = clean_value(get_path(data, [
    "drtype", "dr_type", "resistance_type", "drug_resistance_type", "prediction.drtype"
]))

# This list is used for WHO classification.
# It must contain only canonical TB-Profiler resistant-drug calls.
resistant_drugs_for_classification = []

# This list is used for display and should mirror classification drugs, not all mutation-associated drugs.
resistant_drugs_display_values = []

resistance_mutations = []
key_mutations = []
mutation_evidence_rows = []

variant_blocks = []

for key in [
    "dr_variants",
    "resistance_variants",
    "other_variants",
    "dr_variants_candidate",
    "candidate_variants",
    "variants",
    "mutations"
]:
    block = data.get(key, [])

    if isinstance(block, dict):
        block = list(block.values())

    if isinstance(block, list):
        for x in block:
            if isinstance(x, dict):
                x["_rmap_tb_source_block"] = key
                variant_blocks.append(x)

for item in variant_blocks:
    source_block = item.get("_rmap_tb_source_block", "")

    gene = clean_value(item.get("gene") or item.get("locus_tag") or item.get("locus") or "")
    change = clean_value(
        item.get("change") or
        item.get("protein_change") or
        item.get("nucleotide_change") or
        item.get("original_mutation") or
        item.get("hgvs") or
        item.get("mutation") or ""
    )

    drug = (
        item.get("drug") or
        item.get("drug_name") or
        item.get("name") or
        item.get("drugs") or
        item.get("drug_resistance") or ""
    )

    confidence = clean_value(
        item.get("confidence") or
        item.get("grade") or
        item.get("prediction") or
        item.get("frequency") or ""
    )

    evidence = clean_value(
        item.get("evidence") or
        item.get("annotation") or
        item.get("effect") or
        item.get("type") or ""
    )

    source = clean_value(
        item.get("source") or
        item.get("catalogue") or
        item.get("database") or source_block or ""
    )

    if gene or change:
        key_mutations.append(" ".join([gene, change]).strip())

    drug_values = []

    if drug:
        if isinstance(drug, list):
            drug_values = [clean_value(d) for d in drug if clean_value(d)]
        else:
            drug_values = [clean_value(drug)]
    else:
        drug_values = extract_drug_like_terms(evidence)

    count_for_who = should_count_variant_for_who_classification(
        source_block,
        item,
        drug_values,
        confidence,
        evidence,
        source
    )

    normalized_drugs = []
    for d in drug_values:
        for token in split_drug_tokens(d):
            if token:
                normalized_drugs.append(token)

    if count_for_who:
        for d in normalized_drugs:
            resistant_drugs_for_classification.append(d)
            resistant_drugs_display_values.append(d)

    drug_for_row = "; ".join(drug_values) if drug_values else "No drug assignment reported"

    if drug_values or gene or change:
        mut = []

        if drug_values:
            mut.append("; ".join(drug_values))

        if gene or change:
            mut.append(" ".join([gene, change]).strip())

        if confidence:
            mut.append(confidence)

        if source:
            mut.append(source)

        resistance_mutations.append(" | ".join([x for x in mut if x]))

        mutation_evidence_rows.append([
            sample,
            drug_for_row,
            gene or "Not reported",
            change or "Not reported",
            change or "Not reported",
            confidence or "Not reported",
            evidence or "Not reported",
            json_file
        ])

# Parse structured drug-level phenotype tables if present.
# These are allowed to drive WHO classification only when the item clearly reports resistance.
for key in ["drug_table", "drugs", "resistance"]:
    block = data.get(key)

    if isinstance(block, dict):
        for drug, val in block.items():
            s = json.dumps(val).lower() if isinstance(val, (dict, list)) else str(val).lower()

            if is_positive_resistance_text(s):
                for token in split_drug_tokens(drug):
                    resistant_drugs_for_classification.append(token)
                    resistant_drugs_display_values.append(token)

    elif isinstance(block, list):
        for item in block:
            if not isinstance(item, dict):
                continue

            drug = clean_value(
                item.get("drug") or
                item.get("name") or
                item.get("Drug") or
                item.get("drug_name") or ""
            )

            s = json.dumps(item).lower()

            if drug and is_positive_resistance_text(s):
                for token in split_drug_tokens(drug):
                    resistant_drugs_for_classification.append(token)
                    resistant_drugs_display_values.append(token)

kraken_species, kraken_evidence, kraken_supports_mtbc = get_kraken_support(sample, species_tsv)

raw_lineage_text_for_interpretation = " ".join([main_lineage, sub_lineage]).lower().strip()
has_tbprofiler_lineage = bool(
    re.search(r"(^|[^a-z0-9])(lineage)?[ _-]?[1-9](\.|$|[^0-9])", raw_lineage_text_for_interpretation)
)

has_tbprofiler_species = bool(species)

meaningful_drtype = str(dr_type_raw or "").strip().lower() not in [
    "", "none", "not reported", "unknown", "na", "n/a", "susceptible", "sensitive"
]

canonical_resistant_drugs = uniq_list(resistant_drugs_for_classification)
has_resistance_evidence = bool(canonical_resistant_drugs) or meaningful_drtype

dr_type = classify_who_2021_tb_resistance(canonical_resistant_drugs, dr_type_raw)

if status != "success":
    dr_type = "Resistance not determined by TB-Profiler"
elif not has_resistance_evidence and kraken_supports_mtbc and not has_tbprofiler_lineage and not has_tbprofiler_species:
    dr_type = "Resistance not determined by TB-Profiler"
elif not has_resistance_evidence:
    dr_type = "No resistance detected by TB-Profiler"

if not mutation_evidence_rows:
    if dr_type == "No resistance detected by TB-Profiler":
        mutation_evidence_rows.append([
            sample,
            "No drug-resistance mutation detected by TB-Profiler",
            "Not reported",
            "Not reported",
            "Not reported",
            "Not reported",
            "TB-Profiler completed without reporting mutation-level drug-resistance evidence",
            json_file
        ])
    elif dr_type == "Resistance not determined by TB-Profiler":
        mutation_evidence_rows.append([
            sample,
            "Resistance not determined by TB-Profiler",
            "Not reported",
            "Not reported",
            "Not reported",
            "Not reported",
            "TB-Profiler did not provide interpretable mutation-level drug-resistance evidence for this sample",
            json_file
        ])
    else:
        mutation_evidence_rows.append([
            sample,
            "No mutation-level evidence extracted",
            "Not reported",
            "Not reported",
            "Not reported",
            "Not reported",
            "No mutation-level resistance evidence was available in the TB-Profiler JSON output",
            json_file
        ])

with open("mutation_evidence/tbprofiler_mutation_evidence.tsv", "a") as me:
    for row in mutation_evidence_rows:
        me.write("\t".join([safe_tsv(x) for x in row]) + "\n")

###########################################################################
# Kraken2/Bracken-based MTBC selection for downstream phylogenomics
#
# IMPORTANT:
#   - This is the only logic that determines mtbc_selected, mtbc_reads,
#     and mtbc_samples.txt.
#   - TB-Profiler species and lineage are retained as report annotations only.
#   - TB-Profiler must not select samples for Snippy/core-SNP/IQ-TREE.
#   - IQTREE2_PHYLOGENY later applies a second-stage alignment-quality filter.
###########################################################################

if kraken_supports_mtbc:
    is_mtbc = True

    if has_tbprofiler_lineage:
        reason = "Selected for phylogeny because Kraken2/Bracken species typing supports MTBC; TB-Profiler lineage reported"
    else:
        reason = "Selected for phylogeny because Kraken2/Bracken species typing supports MTBC; TB-Profiler lineage not resolved"
else:
    is_mtbc = False

    if species_tsv and os.path.exists(species_tsv):
        if kraken_species or kraken_evidence:
            reason = "Not selected for phylogeny because Kraken2/Bracken species typing did not support MTBC"
        else:
            reason = "Not selected for phylogeny because no matching Kraken2/Bracken species-typing result supported MTBC"
    else:
        reason = "Not selected for phylogeny because Kraken2/Bracken species-typing results were not provided"

species_display = species

if kraken_supports_mtbc:
    if has_tbprofiler_lineage:
        if species_display:
            species_display = f"{species_display} (Kraken2/Bracken supports MTBC; TB-Profiler lineage reported)"
        else:
            species_display = "Mycobacterium tuberculosis complex (supported by Kraken2/Bracken species typing; TB-Profiler lineage reported)"
    else:
        if species_display:
            species_display = f"{species_display} (Kraken2/Bracken supports MTBC; TB-Profiler lineage not resolved)"
        else:
            species_display = "Mycobacterium tuberculosis complex (supported by Kraken2/Bracken species typing; TB-Profiler lineage not resolved)"
else:
    if species_display:
        species_display = f"{species_display} (not selected for phylogeny; Kraken2/Bracken did not support MTBC)"
    else:
        species_display = "Non-MTB / not classified as MTBC by Kraken2/Bracken"

if is_mtbc:
    os.makedirs("mtbc_reads", exist_ok=True)
    shutil.copy(r1, f"mtbc_reads/{sample}_R1.fastq.gz")
    shutil.copy(r2, f"mtbc_reads/{sample}_R2.fastq.gz")

    with open("mtbc_samples.txt", "a") as fh:
        fh.write(sample + "\n")

main_lineage_display = main_lineage or "Not resolved by TB-Profiler"
sub_lineage_display = sub_lineage or "Not resolved by TB-Profiler"

if str(main_lineage_display).strip().lower() in ["not reported", "unknown", "none", "na", "n/a"]:
    main_lineage_display = "Not resolved by TB-Profiler"

if str(sub_lineage_display).strip().lower() in ["not reported", "unknown", "none", "na", "n/a"]:
    sub_lineage_display = "Not resolved by TB-Profiler"

resistant_drugs_display = uniq(canonical_resistant_drugs) or "None reported"
resistance_mutations_display = uniq(resistance_mutations) or "None reported"
key_mutations_display = uniq(key_mutations) or "None reported"

line = [
    sample,
    species_display or "Not reported",
    main_lineage_display,
    sub_lineage_display,
    dr_type or "Not reported",
    resistant_drugs_display,
    resistance_mutations_display,
    key_mutations_display,
    json_file,
    "YES" if is_mtbc else "NO",
    reason,
    status
]

with open("tbprofiler_summary.tsv", "a") as out:
    out.write("\t".join([safe_tsv(x) for x in line]) + "\n")

resistance_profile_line = [
    sample,
    dr_type or "Not reported",
    resistant_drugs_display,
    resistance_mutations_display,
    key_mutations_display,
    status
]

with open("resistance_profile_summary.tsv", "a") as rout:
    rout.write("\t".join([safe_tsv(x) for x in resistance_profile_line]) + "\n")
PYTB
    done

    python3 - <<'PYHTML'
import csv, html

rows = list(csv.DictReader(open("tbprofiler_summary.tsv"), delimiter="\t"))

cols = [
    ("sample", "Sample ID"),
    ("species", "Species / MTBC member"),
    ("main_lineage", "Main lineage"),
    ("sub_lineage", "Sub-lineage"),
    ("dr_type", "Resistance profile"),
    ("resistant_drugs", "Predicted resistant drugs"),
    ("resistance_mutations", "Resistance-associated mutations"),
    ("key_mutations", "All key mutations"),
    ("mtbc_selected", "Selected for SNP tree by Kraken2/Bracken"),
    ("mtbc_selection_reason", "Kraken2/Bracken selection reason"),
    ("status", "TB-Profiler status")
]

colors = [
    "#0f172a", "#075985", "#7c2d12", "#4c1d95",
    "#991b1b", "#9f1239", "#d97706", "#b45309",
    "#166534", "#365314", "#374151"
]

css = """
body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#f8fafc;color:#0f172a}
.card{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:18px;margin-bottom:18px;box-shadow:0 1px 4px rgba(15,23,42,.08)}
table{border-collapse:collapse;width:100%;font-size:13px;background:white}
th,td{border:1px solid #e2e8f0;padding:9px;vertical-align:top}
th{color:white}
.yes{background:#dcfce7;color:#166534;font-weight:bold}
.no{background:#fee2e2;color:#991b1b;font-weight:bold}
.ok{color:#166534;font-weight:bold}
.fail{color:#b91c1c;font-weight:bold}
.pill{display:inline-block;border-radius:999px;padding:3px 8px;background:#e0f2fe}
.muted{color:#475569}
"""

out = [
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>TB-Profiler MTBC AMR report</title>",
    "<style>" + css + "</style></head><body>",
    "<div class='card'><h1>TB-Profiler drug-resistance, species and lineage report</h1>",
    "<p class='muted'>This table summarizes TB-Profiler JSON outputs, resistance-associated mutations, and Kraken2/Bracken-based MTBC sample selection for downstream core-SNP phylogenomics. TB-Profiler species and lineage are reported as annotations only and do not determine SNP-tree inclusion.</p>",
    "<p class='muted'><strong>Phylogeny selection rule:</strong> samples are selected for Snippy/core-SNP/IQ-TREE only when Kraken2/Bracken species typing supports MTBC. IQ-TREE may later exclude selected samples from the final tree if their core-SNP alignment has excessive missing, ambiguous, or gap content.</p>",
    "<p class='muted'><strong>WHO 2021+ classification rule used:</strong> XDR-TB = MDR/RR-TB plus fluoroquinolone resistance plus bedaquiline or linezolid resistance; Pre-XDR-TB = MDR/RR-TB plus fluoroquinolone resistance; RR-TB = rifampicin resistance without detected isoniazid resistance; MDR-TB = resistance to at least both isoniazid and rifampicin; MDR/RR-TB is the umbrella term for MDR-TB or RR-TB; Hr-TB = isoniazid resistance without rifampicin resistance.</p>",
    "</div>",
    "<div class='card'><table><thead><tr>"
]

for i, (_, label) in enumerate(cols):
    out.append(f"<th style='background:{colors[i]}'>{html.escape(label)}</th>")

out.append("</tr></thead><tbody>")

for r in rows:
    out.append("<tr>")

    for key, _ in cols:
        val = html.escape(r.get(key, "") or "")

        if key == "mtbc_selected":
            cls = "yes" if val == "YES" else "no"
            out.append(f"<td class='{cls}'>{val}</td>")

        elif key == "status":
            cls = "ok" if val == "success" else "fail"
            out.append(f"<td class='{cls}'>{val}</td>")

        elif key in ("main_lineage", "sub_lineage", "dr_type"):
            out.append(f"<td><span class='pill'>{val}</span></td>")

        else:
            out.append(f"<td>{val}</td>")

    out.append("</tr>")

out.append("</tbody></table></div></body></html>")

open("tbprofiler_combined_report.html", "w").write("\n".join(out))
PYHTML

    if [ ! -s tbprofiler_combined_report.html ]; then
      echo "WARNING: tbprofiler_combined_report.html was not produced; writing fallback report." >> logs/tbprofiler.command.log
      python3 - <<'PYHTML_FALLBACK' || true
import html, pathlib
summary = pathlib.Path("tbprofiler_summary.tsv")
body = ""
if summary.exists():
    body = "<h2>Raw TB-Profiler summary TSV</h2><pre>" + html.escape(summary.read_text()) + "</pre>"
else:
    body = "<p>tbprofiler_summary.tsv was not found.</p>"
pathlib.Path("tbprofiler_combined_report.html").write_text(
    "<!doctype html><html><head><meta charset='utf-8'><title>TB-Profiler MTBC AMR report</title>"
    "<style>body{font-family:Arial,Helvetica,sans-serif;margin:24px}pre{white-space:pre-wrap;background:#f8fafc;border:1px solid #e2e8f0;padding:12px}</style>"
    "</head><body><h1>TB-Profiler drug-resistance, species and lineage report</h1>"
    "<p><strong>Note:</strong> The styled HTML renderer did not complete, so this fallback report preserves the raw summary table for Cromwell output collection.</p>"
    + body + "</body></html>"
)
PYHTML_FALLBACK
    fi

    python3 - <<'PYMUTHTML'
import csv, html

rows = list(csv.DictReader(open("mutation_evidence/tbprofiler_mutation_evidence.tsv"), delimiter="\t"))

css = """
body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#f8fafc;color:#0f172a}
.card{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:18px;margin-bottom:18px;box-shadow:0 1px 4px rgba(15,23,42,.08)}
table{border-collapse:collapse;width:100%;font-size:13px;background:white}
th,td{border:1px solid #e2e8f0;padding:9px;vertical-align:top}
th{background:#7c2d12;color:white}
.badge-drug{display:inline-block;border-radius:999px;padding:4px 8px;background:#fee2e2;color:#991b1b;font-weight:700;font-size:12px}
.badge-gene{display:inline-block;border-radius:999px;padding:4px 8px;background:#e0f2fe;color:#075985;font-weight:700;font-size:12px}
.badge-confidence{display:inline-block;border-radius:999px;padding:4px 8px;background:#dcfce7;color:#166534;font-weight:700;font-size:12px}
.muted{color:#475569}
"""

out = [
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>TB-Profiler mutation-level resistance evidence</title>",
    "<style>" + css + "</style></head><body>",
    "<div class='card'><h1>Resistance Mutation Evidence Summary</h1>",
    "<p class='muted'>This table extracts mutation-level drug-resistance evidence from TB-Profiler JSON outputs, including drug or evidence source, gene, mutation/change, confidence, and evidence fields where reported. Classification is not inferred from candidate, benign, lineage, synonymous, or low-confidence evidence rows.</p></div>",
    "<div class='card'><table><thead><tr>",
    "<th>Sample ID</th>",
    "<th>Drug / Evidence source</th>",
    "<th>Gene</th>",
    "<th>Mutation</th>",
    "<th>Change</th>",
    "<th>Confidence</th>",
    "<th>Evidence / associated drug(s)</th>",
    "</tr></thead><tbody>"
]

for r in rows:
    sample = html.escape(r.get("sample", "") or "")
    drug = html.escape(r.get("drug", "") or "")
    gene = html.escape(r.get("gene", "") or "")
    mutation = html.escape(r.get("mutation", "") or "")
    change = html.escape(r.get("change", "") or "")
    confidence = html.escape(r.get("confidence", "") or "")
    evidence = html.escape(r.get("evidence", "") or "")

    out.append("<tr>")
    out.append(f"<td>{sample}</td>")
    out.append(f"<td><span class='badge-drug'>{drug}</span></td>")
    out.append(f"<td><span class='badge-gene'>{gene}</span></td>")
    out.append(f"<td>{mutation}</td>")
    out.append(f"<td>{change}</td>")
    out.append(f"<td><span class='badge-confidence'>{confidence}</span></td>")
    out.append(f"<td>{evidence}</td>")
    out.append("</tr>")

out.append("</tbody></table></div></body></html>")

open("mutation_evidence/tbprofiler_mutation_evidence.html", "w").write("\n".join(out))
PYMUTHTML

    if [ ! -s mutation_evidence/tbprofiler_mutation_evidence.html ]; then
      echo "WARNING: tbprofiler_mutation_evidence.html was not produced; writing fallback mutation-evidence report." >> logs/tbprofiler.command.log
      python3 - <<'PYMUTHTML_FALLBACK' || true
import html, pathlib
tsv = pathlib.Path("mutation_evidence/tbprofiler_mutation_evidence.tsv")
body = ""
if tsv.exists():
    body = "<h2>Raw mutation evidence TSV</h2><pre>" + html.escape(tsv.read_text()) + "</pre>"
else:
    body = "<p>tbprofiler_mutation_evidence.tsv was not found.</p>"
pathlib.Path("mutation_evidence/tbprofiler_mutation_evidence.html").write_text(
    "<!doctype html><html><head><meta charset='utf-8'><title>TB-Profiler mutation evidence</title>"
    "<style>body{font-family:Arial,Helvetica,sans-serif;margin:24px}pre{white-space:pre-wrap;background:#f8fafc;border:1px solid #e2e8f0;padding:12px}</style>"
    "</head><body><h1>Resistance Mutation Evidence Summary</h1>"
    "<p><strong>Note:</strong> The styled mutation-evidence renderer did not complete, so this fallback report preserves the raw TSV for Cromwell output collection.</p>"
    + body + "</body></html>"
)
PYMUTHTML_FALLBACK
    fi
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 500 HDD"
    timeout: "120 hours"
  }

  output {
    Array[File] json_reports = glob("tbprofiler_results/**/*.json")
    Array[File] txt_reports = glob("tbprofiler_results/**/*.txt")
    Array[File] tbprofiler_logs = glob("logs/*.tbprofiler.log")
    File tbprofiler_command_log = "logs/tbprofiler.command.log"
    File summary_tsv = "tbprofiler_summary.tsv"
    File resistance_profile_summary_tsv = "resistance_profile_summary.tsv"
    File combined_html = "tbprofiler_combined_report.html"
    File mutation_evidence_tsv = "mutation_evidence/tbprofiler_mutation_evidence.tsv"
    File mutation_evidence_html = "mutation_evidence/tbprofiler_mutation_evidence.html"
    File mtbc_samples_txt = "mtbc_samples.txt"
    Array[File] mtbc_reads = glob("mtbc_reads/*.fastq.gz")
  }
}
task SNIPPY_CORE_MTBC {
  input {
    String docker_image = "staphb/snippy:4.6.0"

    # These reads should already be filtered upstream by TB_PROFILER_AND_MTBC_FILTER.
    # Selection for this task is Kraken2/Bracken-based, not TB-Profiler-based.
    # TB-Profiler species/lineage/resistance are annotations only and must not
    # determine which samples enter this task.
    Array[File]+ input_reads

    File reference_genome
    String reference_type = "genbank"
    Int cpu = 16
    Int memory_gb = 64
    Int min_quality = 20
  }

  command <<<
    set -uo pipefail
    mkdir -p snippy_results snippy_core logs

    echo "SNIPPY_CORE_MTBC started." > logs/snippy.command.log
    echo "Selection rule: input_reads are expected to be Kraken2/Bracken-selected MTBC reads from TB_PROFILER_AND_MTBC_FILTER.mtbc_reads." >> logs/snippy.command.log
    echo "TB-Profiler species, lineage, and resistance annotations do not determine entry into this task." >> logs/snippy.command.log
    echo "IQTREE2_PHYLOGENY will later apply alignment-quality filtering before final tree inference." >> logs/snippy.command.log
    echo "" >> logs/snippy.command.log

    if command -v snippy >/dev/null 2>&1; then
      SNIPPY_BIN="$(command -v snippy)"
    elif [ -x /usr/local/bin/snippy ]; then
      SNIPPY_BIN="/usr/local/bin/snippy"
    else
      echo "ERROR: snippy executable not found" >&2
      exit 127
    fi

    if command -v snippy-core >/dev/null 2>&1; then
      SNIPPY_CORE_BIN="$(command -v snippy-core)"
    else
      SNIPPY_CORE_BIN="snippy-core"
    fi

    echo "Using snippy: ${SNIPPY_BIN}" >> logs/snippy.command.log
    echo "Using snippy-core: ${SNIPPY_CORE_BIN}" >> logs/snippy.command.log
    echo "Minimum variant quality: ~{min_quality}" >> logs/snippy.command.log
    echo "CPU threads: ~{cpu}" >> logs/snippy.command.log
    echo "" >> logs/snippy.command.log

    if command -v samtools >/dev/null 2>&1; then
      echo "samtools detected:" > logs/samtools.version.log
      samtools --version >> logs/samtools.version.log 2>&1 || true
    else
      echo "WARNING: samtools not found; mean depth will be reported as NA." > logs/samtools.version.log
    fi

    if [ "~{reference_type}" = "genbank" ]; then
      cp "~{reference_genome}" reference_input.gbk
      ref="reference_input.gbk"
      grep -q '^LOCUS' "$ref" || { echo "ERROR: invalid GenBank reference"; exit 1; }
    else
      cp "~{reference_genome}" reference_input.fa
      ref="reference_input.fa"
      grep -q '^>' "$ref" || { echo "ERROR: invalid FASTA reference"; exit 1; }
    fi

    files=(~{sep=' ' input_reads})
    n=${#files[@]}

    echo "Number of Kraken2/Bracken-selected MTBC FASTQ files received: ${n}" >> logs/snippy.command.log

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: MTBC reads must be paired R1/R2 files." >&2
      echo "ERROR: Odd number of input FASTQ files received by SNIPPY_CORE_MTBC." >> logs/snippy.command.log
      exit 1
    fi

    paired_count=$((n / 2))
    echo "Number of Kraken2/Bracken-selected paired MTBC samples received: ${paired_count}" >> logs/snippy.command.log
    echo "" >> logs/snippy.command.log

    echo -e "sample\tstatus\tvcf\taligned_fasta\tselection_basis" > variant_summary.tsv
    echo -e "Sample\tMeanDepth" > mean_depth_summary.tsv
    echo -e "Sample\tSelectedForSnippyCore\tSelectionBasis" > snippy_core_selection_summary.tsv

    successful_samples=()

    for ((i=0; i<n; i+=2)); do
      R1="${files[$i]}"
      R2="${files[$((i+1))]}"
      sample=$(basename "$R1" | sed -E 's/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//' | sed -E 's/(_R?1|_1|\.R?1|\.1)(_|$).*//')
      outdir="snippy_results/${sample}"

      selection_basis="Selected for Snippy/core-SNP analysis because upstream Kraken2/Bracken species typing supported MTBC"

      echo "Running snippy for ${sample}" >> logs/snippy.command.log
      echo -e "${sample}\tYES\t${selection_basis}" >> snippy_core_selection_summary.tsv

      if "$SNIPPY_BIN" \
        --cpus ~{cpu} \
        --minqual ~{min_quality} \
        --ref "$ref" \
        --R1 "$R1" \
        --R2 "$R2" \
        --outdir "$outdir" \
        --prefix "$sample" \
        --force \
        > "logs/${sample}.snippy.log" 2>&1; then

        status="success"
        successful_samples+=("$outdir")

      else
        status="snippy_failed"
        echo "WARNING: snippy failed for ${sample}" >> logs/snippy.command.log
      fi

      vcf="${outdir}/${sample}.vcf"
      aln="${outdir}/${sample}.aligned.fa"

      if [ ! -f "$vcf" ]; then vcf="NA"; fi
      if [ ! -f "$aln" ]; then aln="NA"; fi

      depth="NA"

      if command -v samtools >/dev/null 2>&1; then
        bam_candidates=(
          "${outdir}/${sample}.bam"
          "${outdir}/snps.bam"
        )

        for bam in "${bam_candidates[@]}"; do
          if [ -f "$bam" ]; then
            depth=$(samtools depth "$bam" 2>/dev/null | awk '{sum+=$3} END {if(NR>0) printf "%.2f", sum/NR; else print "NA"}')
            break
          fi
        done
      fi

      if [ -z "$depth" ]; then
        depth="NA"
      fi

      echo -e "${sample}\t${status}\t${vcf}\t${aln}\t${selection_basis}" >> variant_summary.tsv
      echo -e "${sample}\t${depth}" >> mean_depth_summary.tsv
    done

    core_status="skipped"

    if [ ${#successful_samples[@]} -ge 2 ]; then
      echo "Running snippy-core on ${#successful_samples[@]} successfully processed Kraken2/Bracken-selected MTBC samples" >> logs/snippy.command.log

      if "$SNIPPY_CORE_BIN" --ref "$ref" --prefix snippy_core/core "${successful_samples[@]}" >> logs/snippy_core.log 2>&1; then
        core_status="success"
      else
        core_status="snippy_core_failed"
        echo "WARNING: snippy-core failed" >> logs/snippy.command.log
      fi
    else
      echo "WARNING: fewer than 2 successful Kraken2/Bracken-selected MTBC samples, skipping snippy-core" >> logs/snippy.command.log
    fi

    echo "Snippy-core status: ${core_status}" >> logs/snippy.command.log

    if [ ! -f mean_depth_summary.tsv ]; then
      echo -e "Sample\tMeanDepth" > mean_depth_summary.tsv
    fi

    python3 - <<'PY'
import csv, html

rows = list(csv.DictReader(open("variant_summary.tsv"), delimiter="\t"))

depth_by_sample = {}
try:
    with open("mean_depth_summary.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sample = r.get("Sample") or r.get("sample") or ""
            depth = r.get("MeanDepth") or r.get("mean_depth") or "NA"
            if sample:
                depth_by_sample[sample] = depth
except Exception:
    depth_by_sample = {}

selection_rows = []
try:
    with open("snippy_core_selection_summary.tsv") as fh:
        selection_rows = list(csv.DictReader(fh, delimiter="\t"))
except Exception:
    selection_rows = []

selection_basis_by_sample = {}
for r in selection_rows:
    sample = r.get("Sample") or ""
    basis = r.get("SelectionBasis") or ""
    if sample:
        selection_basis_by_sample[sample] = basis

out = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>MTBC Snippy-core summary</title>
<style>
body{font-family:Arial;margin:24px;background:#f8fafc;color:#111827}
table{border-collapse:collapse;width:100%;background:white}
th,td{border:1px solid #ddd;padding:8px;vertical-align:top}
th{background:#1d4ed8;color:white}
.ok{color:#166534;font-weight:bold}
.fail{color:#b91c1c;font-weight:bold}
.note{background:#eef6ff;border-left:5px solid #2563eb;padding:12px;border-radius:10px;margin:12px 0;line-height:1.45}
.badge-green{display:inline-block;border-radius:999px;padding:4px 8px;background:#28A745;color:white;font-size:12px;font-weight:700}
.badge-red{display:inline-block;border-radius:999px;padding:4px 8px;background:#b91c1c;color:white;font-size:12px;font-weight:700}
</style>
</head>
<body>
<h1>MTBC core-SNP variant-calling summary</h1>

<div class="note">
<strong>Selection rule:</strong> samples entering this Snippy/core-SNP task were selected upstream using Kraken2/Bracken MTBC support from the species-typing step. TB-Profiler species, lineage, and resistance outputs are annotations only and do not determine which samples enter Snippy/core-SNP/IQ-TREE analysis.
<br><br>
<strong>Second-stage tree filtering:</strong> IQTREE2_PHYLOGENY may later exclude selected samples from final IQ-TREE inference if the core-SNP alignment has no usable A/C/G/T bases or excessive missing, ambiguous, or gap content.
</div>

<table>
<thead>
<tr>
<th>Sample</th>
<th>Status</th>
<th>Mean depth</th>
<th>Selected for Snippy/core-SNP by Kraken2/Bracken</th>
<th>Selection basis</th>
<th>VCF</th>
<th>Aligned FASTA</th>
</tr>
</thead>
<tbody>
"""

for r in rows:
    status = r["status"]
    cls = "ok" if status == "success" else "fail"
    sample = r["sample"]
    mean_depth = depth_by_sample.get(sample, "NA")
    selection_basis = r.get("selection_basis") or selection_basis_by_sample.get(sample, "Selected upstream by Kraken2/Bracken MTBC support")

    out += (
        "<tr>"
        f"<td>{html.escape(sample)}</td>"
        f"<td class='{cls}'>{html.escape(status)}</td>"
        f"<td>{html.escape(mean_depth)}</td>"
        "<td><span class='badge-green'>YES</span></td>"
        f"<td>{html.escape(selection_basis)}</td>"
        f"<td>{html.escape(r['vcf'])}</td>"
        f"<td>{html.escape(r['aligned_fasta'])}</td>"
        "</tr>\n"
    )

out += "</tbody></table></body></html>"

open("variant_summary.html","w").write(out)
PY
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 500 HDD"
    timeout: "120 hours"
  }

  output {
    Array[File] vcf_files = glob("snippy_results/*/*.vcf")
    Array[File] aligned_fastas = glob("snippy_results/*/*.aligned.fa")
    Array[File] snippy_tab_files = glob("snippy_results/*/*.tab")
    Array[File] snippy_logs = glob("logs/*.snippy.log")
    File snippy_command_log = "logs/snippy.command.log"
    File samtools_version_log = "logs/samtools.version.log"

    File variant_summary = "variant_summary.html"
    File mean_depth_summary_tsv = "mean_depth_summary.tsv"
    File snippy_core_selection_summary_tsv = "snippy_core_selection_summary.tsv"

    File? core_full_alignment = "snippy_core/core.full.aln"
    File? core_snp_alignment = "snippy_core/core.aln"
    File? core_tab = "snippy_core/core.tab"
    File? core_vcf = "snippy_core/core.vcf"
  }
}

task TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS {
  input {
    String docker_image = "python:3.11-slim"
    Array[File] snippy_tab_files
    String genes_csv = "rpoB,katG,inhA,fabG1,ahpC,embB,pncA,rpsL,rrs,gyrA,gyrB,eis,ethA,ethR,thyA,folC,alr,ddl,gidB,tlyA,rrl,atpE,rv0678,pepQ"
  }

  command <<<
    set -uo pipefail
    mkdir -p nonsynonymous_drug_gene_mutations logs

    echo "Input snippy tab files:" > logs/nonsyn.log
    printf "%s\n" ~{sep='\n' snippy_tab_files} >> logs/nonsyn.log || true

    python3 - <<'PY'
import csv, html, sys
from pathlib import Path

files = [Path(x) for x in """~{sep='\n' snippy_tab_files}""".splitlines() if x.strip()]
genes = {g.strip().lower() for g in "~{genes_csv}".split(',') if g.strip()}

keep_terms = [
    'missense','stop_gained','stop_lost','start_lost',
    'frameshift','inframe_insertion','inframe_deletion',
    'disruptive_inframe','conservative_inframe',
    'protein_altering','coding_sequence_variant'
]
exclude_terms = ['synonymous_variant']

rows = []

if not files:
    rows.append({
        'sample':'NA',
        'gene':'NA',
        'position':'',
        'ref':'',
        'alt':'',
        'type':'',
        'effect':'No Snippy .tab files provided — mutation analysis skipped',
        'aa_change':'',
        'nt_change':'',
        'product':'',
        'evidence':''
    })

for tab in files:
    sample = tab.stem.split('.')[0]

    try:
        if not tab.exists():
            rows.append({
                'sample': sample,
                'gene':'ERROR',
                'position':'',
                'ref':'',
                'alt':'',
                'type':'',
                'effect': f'{tab} missing',
                'aa_change':'',
                'nt_change':'',
                'product':'',
                'evidence':''
            })
            continue

        with open(tab, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')

            for r in reader:
                gene = (r.get('GENE') or r.get('LOCUS_TAG') or '').strip()
                effect = (r.get('EFFECT') or '').strip()
                effect_l = effect.lower()

                if not gene or gene.lower() not in genes:
                    continue

                if any(x in effect_l for x in exclude_terms):
                    continue

                if not any(x in effect_l for x in keep_terms):
                    continue

                rows.append({
                    'sample': sample,
                    'gene': gene,
                    'position': r.get('POS',''),
                    'ref': r.get('REF',''),
                    'alt': r.get('ALT',''),
                    'type': r.get('TYPE',''),
                    'effect': effect,
                    'aa_change': r.get('AA_POS','') or r.get('AA_CHANGE',''),
                    'nt_change': r.get('NT_POS','') or r.get('NT_CHANGE',''),
                    'product': r.get('PRODUCT',''),
                    'evidence': r.get('EVIDENCE','')
                })

    except Exception as e:
        rows.append({
            'sample': sample,
            'gene':'ERROR',
            'position':'',
            'ref':'',
            'alt':'',
            'type':'',
            'effect': f'Parse error: {e}',
            'aa_change':'',
            'nt_change':'',
            'product':'',
            'evidence':''
        })

fields = ['sample','gene','position','ref','alt','type','effect','aa_change','nt_change','product','evidence']

out_dir = Path('nonsynonymous_drug_gene_mutations')
out_dir.mkdir(exist_ok=True)

tsv_file = out_dir / 'nonsynonymous_drug_gene_mutations.tsv'
with open(tsv_file, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader()
    w.writerows(rows)

css = '''
body{font-family:Arial;margin:24px;background:#f8fafc;color:#111827}
.card{background:white;border:1px solid #e5e7eb;border-radius:14px;padding:18px;box-shadow:0 1px 6px rgba(0,0,0,.08)}
table{border-collapse:collapse;width:100%;font-size:13px}
th{color:white;padding:9px;text-align:left}
td{border-bottom:1px solid #e5e7eb;padding:8px}
.gene{font-weight:bold;color:#7c2d12}
.mut{color:#b91c1c;font-weight:bold}
.warn{color:#d97706;font-weight:bold}
'''
colors = ['#0f766e','#7c2d12','#2563eb','#374151','#374151','#4c1d95','#b91c1c','#d97706','#6d28d9','#087f5b','#6b7280']

labels = ['Sample','Gene','Position','REF','ALT','Type','Effect','AA change','NT change','Product','Evidence']

html_out = [
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>Non-synonymous mutations</title>",
    f"<style>{css}</style></head><body>",
    "<div class='card'><h1>Non-synonymous mutations in TB drug-resistance genes</h1>",
    "<table><thead><tr>"
]

for i, lab in enumerate(labels):
    html_out.append(f"<th style='background:{colors[i]}'>{html.escape(lab)}</th>")

html_out.append("</tr></thead><tbody>")

for r in rows:
    html_out.append("<tr>")
    for f in fields:
        val = html.escape(str(r.get(f,"")))
        cls = "gene" if f=="gene" else "mut" if f in ("effect","alt") else "warn" if "error" in val.lower() else ""
        html_out.append(f"<td class='{cls}'>{val}</td>")
    html_out.append("</tr>")

html_out.append("</tbody></table></div></body></html>")

(out_dir / 'nonsynonymous_drug_gene_mutations.html').write_text("\n".join(html_out))
PY
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File nonsynonymous_mutations_tsv = "nonsynonymous_drug_gene_mutations/nonsynonymous_drug_gene_mutations.tsv"
    File nonsynonymous_mutations_html = "nonsynonymous_drug_gene_mutations/nonsynonymous_drug_gene_mutations.html"
    File nonsynonymous_log = "logs/nonsyn.log"
  }
}

task SNP_DISTANCE_CLUSTERING {
  input {
    String docker_image = "python:3.11-slim"
    File core_full_alignment
    Int likely_transmission_snp_threshold = 5
    Int possible_transmission_snp_threshold = 12
    Int cpu = 8
    Int memory_gb = 64
  }

  command <<<
    set -uo pipefail

    mkdir -p snp_distance logs

    cp "~{core_full_alignment}" snp_distance/core.full.aln

    echo "Running fail-safe pure-Python SNP distance clustering on core.full.aln" > logs/snp_distance.command.log
    echo "Reference/non-sample sequences are excluded from sample-level SNP distance and cluster reporting." >> logs/snp_distance.command.log

    echo "Detecting system architecture and CPU features..." >> logs/snp_distance.command.log

    ARCH="$(uname -m)"
    echo "Architecture detected: ${ARCH}" >> logs/snp_distance.command.log

    if command -v lscpu >/dev/null 2>&1; then
      CPU_FLAGS="$(lscpu | grep -i 'Flags' || true)"
    else
      CPU_FLAGS="$(cat /proc/cpuinfo | grep -m1 -i 'flags' || true)"
    fi

    echo "CPU flags detected:" >> logs/snp_distance.command.log
    echo "${CPU_FLAGS}" >> logs/snp_distance.command.log

    apt-get update >/dev/null 2>&1 || true
    apt-get install -y --no-install-recommends \
      python3-pip \
      python3-dev \
      build-essential \
      >/dev/null 2>&1 || true

    python3 -m pip install --upgrade pip setuptools wheel >/dev/null 2>&1 || true

    if [ "${ARCH}" = "x86_64" ]; then
      echo "Using Intel/x86_64-safe Python package installation..." >> logs/snp_distance.command.log

      # Helps avoid illegal-instruction crashes from optimized BLAS/NumPy builds
      export OPENBLAS_CORETYPE=Haswell
      export NUMPY_EXPERIMENTAL_ARRAY_FUNCTION=0

      python3 -m pip install --quiet \
        "numpy<2.0" \
        "pandas<2.2" \
        "matplotlib<3.9" \
        "seaborn<0.13" \
        || {
          echo "WARNING: Standard x86_64 pip install failed. Retrying with no binary for numpy..." >> logs/snp_distance.command.log
          python3 -m pip install --quiet --no-binary=numpy "numpy<2.0" || true
          python3 -m pip install --quiet "pandas<2.2" "matplotlib<3.9" "seaborn<0.13" || true
        }

    elif [ "${ARCH}" = "aarch64" ] || [ "${ARCH}" = "arm64" ]; then
      echo "Using ARM64/aarch64-safe Python package installation..." >> logs/snp_distance.command.log

      python3 -m pip install --quiet \
        "numpy<2.0" \
        "pandas<2.2" \
        "matplotlib<3.9" \
        "seaborn<0.13" \
        || {
          echo "WARNING: ARM64 pip install failed. Retrying with source-compatible installation..." >> logs/snp_distance.command.log
          python3 -m pip install --quiet --no-binary=numpy "numpy<2.0" || true
          python3 -m pip install --quiet "pandas<2.2" "matplotlib<3.9" "seaborn<0.13" || true
        }

    else
      echo "Unknown architecture: ${ARCH}. Using conservative package installation..." >> logs/snp_distance.command.log

      python3 -m pip install --quiet \
        "numpy<2.0" \
        "pandas<2.2" \
        "matplotlib<3.9" \
        "seaborn<0.13" \
        || true
    fi

    echo "Checking Python package imports..." >> logs/snp_distance.command.log

    python3 - <<'PY' >> logs/snp_distance.command.log 2>&1
import sys

packages = ["numpy", "pandas", "matplotlib", "seaborn"]

for pkg in packages:
    try:
        __import__(pkg)
        print(f"OK: {pkg} imported successfully")
    except Exception as e:
        print(f"WARNING: {pkg} could not be imported: {e}")

print("Python package check completed.")
PY

    python3 - <<'PY' >> logs/snp_distance.command.log 2>&1
import csv
import html
import traceback
from pathlib import Path

alignment_file = Path("snp_distance/core.full.aln")
matrix_file = Path("snp_distance/pairwise_snp_distance_matrix.tsv")
pairs_file = Path("snp_distance/pairwise_snp_distance_pairs.tsv")
clusters_file = Path("snp_distance/snp_cluster_summary.tsv")
html_file = Path("snp_distance/snp_distance_cluster_report.html")
status_file = Path("snp_distance/snp_distance_status.txt")
excluded_file = Path("snp_distance/excluded_reference_or_non_sample_sequences.tsv")
heatmap_png = Path("snp_distance/pairwise_snp_heatmap.png")

likely_threshold = int("~{likely_transmission_snp_threshold}")
possible_threshold = int("~{possible_transmission_snp_threshold}")

REFERENCE_NAMES = {
    "reference",
    "ref",
    "h37rv",
    "h37rvsiena",
    "nc_000962",
    "mycobacterium_tuberculosis_h37rv"
}


def is_reference_name(name):
    n = str(name or "").strip().lower()
    return n in REFERENCE_NAMES or n.startswith("reference") or n.startswith("ref_")


def write_fallback_outputs(message):
    matrix_file.write_text("sample\n", encoding="utf-8")
    pairs_file.write_text(
        "sample1\tsample2\tsnp_distance\tcomparable_sites\tinterpretation\tcluster_class\n",
        encoding="utf-8"
    )
    clusters_file.write_text(
        "cluster_id\tsample1\tsample2\tsnp_distance\tinterpretation\n",
        encoding="utf-8"
    )
    excluded_file.write_text("sequence\treason\n", encoding="utf-8")
    status_file.write_text(str(message) + "\n", encoding="utf-8")

    html_file.write_text(
        "<!doctype html>\n"
        "<html><head><meta charset='utf-8'>"
        "<title>Pairwise SNP Distance and Cluster Summary</title>"
        "<style>"
        "body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#f8fafc;color:#0f172a}"
        ".card{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:18px;margin-bottom:18px}"
        ".note{background:#fff7ed;border-left:5px solid #f97316;padding:12px;border-radius:10px}"
        "</style>"
        "</head><body><div class='card'>"
        "<h1>Pairwise SNP Distance and Cluster Summary</h1>"
        f"<div class='note'>{html.escape(str(message))}</div>"
        "</div></body></html>\n",
        encoding="utf-8"
    )


def read_alignment(path):
    text = path.read_text(errors="replace").splitlines()

    records = {}
    name = None
    seq = []

    is_fasta = any(line.startswith(">") for line in text if line.strip())

    if is_fasta:
        for line in text:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).upper()

                name = line[1:].split()[0]
                seq = []

            else:
                seq.append(line)

        if name is not None:
            records[name] = "".join(seq).upper()

        return records

    for line in text:
        line = line.strip()

        if not line:
            continue

        parts = line.split()

        if len(parts) >= 2 and set(parts[1].upper()).issubset(set("ACGTNRYKMSWBDHV-.?")):
            records[parts[0]] = "".join(parts[1:]).upper()

    return records


def snp_distance(seq1, seq2):
    missing = {"N", "-", "?", "."}
    comparable = set("ACGT")

    n = min(len(seq1), len(seq2))

    dist = 0
    compared = 0

    for i in range(n):
        a = seq1[i]
        b = seq2[i]

        if a in missing or b in missing:
            continue

        if a not in comparable or b not in comparable:
            continue

        compared += 1

        if a != b:
            dist += 1

    return dist, compared


def interpret(distance):
    if distance <= likely_threshold:
        return (
            "Genomically close; review epidemiological linkage",
            "close"
        )

    if distance <= possible_threshold:
        return (
            "Intermediate SNP distance; review with epidemiological metadata",
            "intermediate"
        )

    return (
        "Not clustered by SNP threshold",
        "distant"
    )


def write_heatmap_png(matrix, samples):
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        from matplotlib.colors import LinearSegmentedColormap

        matrix_df = pd.DataFrame(
            matrix,
            index=samples,
            columns=samples
        )

        cmap = LinearSegmentedColormap.from_list(
            "snp_scale",
            ["#d1fae5", "#fdba74", "#fecaca"]
        )

        plt.figure(figsize=(12, 10))

        sns.heatmap(
            matrix_df,
            annot=True,
            fmt=".0f",
            cmap=cmap,
            linewidths=0.5,
            cbar=False,
            square=True
        )

        plt.title(
            "Pairwise SNP Distance Heatmap",
            fontsize=34,
            fontweight="bold",
            pad=45
        )

        plt.xticks(
            rotation=60,
            ha="left",
            fontsize=18
        )

        plt.yticks(
            rotation=0,
            fontsize=18
        )

        plt.figtext(
            0.5,
            0.02,
            "Reference sequence excluded. Lower values indicate closer genomic relatedness.",
            ha="center",
            fontsize=16,
            color="#475569"
        )

        plt.tight_layout(rect=[0, 0.06, 1, 0.94])

        plt.savefig(
            heatmap_png,
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

    except Exception as exc:
        print("WARNING: Heatmap PNG generation failed:", str(exc))
        heatmap_png.write_bytes(b"")


try:
    records_all = read_alignment(alignment_file)

    excluded = []
    records = {}

    for name, seq in records_all.items():
        if is_reference_name(name):
            excluded.append((
                name,
                "Reference/non-sample sequence excluded from sample-level SNP distance analysis"
            ))
        else:
            records[name] = seq

    with excluded_file.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sequence", "reason"])
        writer.writerows(excluded)

    samples = list(records.keys())

    if len(records_all) == 0:
        write_fallback_outputs(
            "No sequences were found in the core alignment; SNP distance outputs were generated as empty placeholders."
        )
        raise SystemExit(0)

    if len(samples) == 0:
        write_fallback_outputs(
            "Only reference/non-sample sequences were found after filtering; no sample-level SNP distances were generated."
        )
        raise SystemExit(0)

    if len(samples) == 1:
        sample = samples[0]

        matrix_file.write_text(
            "sample\t" + sample + "\n" + sample + "\t0\n",
            encoding="utf-8"
        )

        pairs_file.write_text(
            "sample1\tsample2\tsnp_distance\tcomparable_sites\tinterpretation\tcluster_class\n",
            encoding="utf-8"
        )

        clusters_file.write_text(
            "cluster_id\tsample1\tsample2\tsnp_distance\tinterpretation\n",
            encoding="utf-8"
        )

        status_file.write_text(
            "Only one sample sequence was present after reference filtering; no pairwise SNP distances were calculated.\n",
            encoding="utf-8"
        )

        html_file.write_text(
            "<!doctype html><html><head><meta charset='utf-8'>"
            "<title>Pairwise SNP Distance and Cluster Summary</title>"
            "<style>"
            "body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#f8fafc;color:#0f172a}"
            ".card{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:18px;margin-bottom:18px}"
            ".note{background:#eff6ff;border-left:5px solid #2563eb;padding:12px;border-radius:10px}"
            "</style>"
            "</head>"
            "<body>"
            "<div class='card'>"
            "<h1>Pairwise SNP Distance and Cluster Summary</h1>"
            "<div class='note'>Only one sample sequence was present after reference filtering; pairwise SNP distances could not be generated.</div>"
            "</div>"
            "</body></html>\n",
            encoding="utf-8"
        )

        raise SystemExit(0)

    lengths = {len(v) for v in records.values()}

    if len(lengths) > 1:
        print(
            "WARNING: Alignment sequences have unequal lengths; distances calculated over shared aligned positions."
        )

    matrix = []
    comparable_matrix = []

    for s1 in samples:
        row = []
        comparable_row = []

        for s2 in samples:
            d, c = snp_distance(records[s1], records[s2])

            row.append(d)
            comparable_row.append(c)

        matrix.append(row)
        comparable_matrix.append(comparable_row)

    with matrix_file.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")

        writer.writerow(["sample"] + samples)

        for sample, row in zip(samples, matrix):
            writer.writerow([sample] + row)

    pairs = []

    for i, sample1 in enumerate(samples):
        for j, sample2 in enumerate(samples):

            if j <= i:
                continue

            distance = matrix[i][j]
            compared_sites = comparable_matrix[i][j]

            interpretation, cluster_class = interpret(distance)

            pairs.append({
                "sample1": sample1,
                "sample2": sample2,
                "snp_distance": distance,
                "comparable_sites": compared_sites,
                "interpretation": interpretation,
                "cluster_class": cluster_class
            })

    with pairs_file.open("w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "sample1",
                "sample2",
                "snp_distance",
                "comparable_sites",
                "interpretation",
                "cluster_class"
            ],
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(pairs)

    close_pairs = [
        p for p in pairs
        if p["cluster_class"] in ["close", "intermediate"]
    ]

    with clusters_file.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")

        writer.writerow([
            "cluster_id",
            "sample1",
            "sample2",
            "snp_distance",
            "interpretation"
        ])

        for idx, p in enumerate(close_pairs, start=1):
            writer.writerow([
                f"Cluster_{idx}",
                p["sample1"],
                p["sample2"],
                p["snp_distance"],
                p["interpretation"]
            ])

    write_heatmap_png(matrix, samples)

    css = """
body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#f8fafc;color:#0f172a}
.card{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:18px;margin-bottom:18px;box-shadow:0 1px 4px rgba(15,23,42,.08)}
table{border-collapse:collapse;width:100%;font-size:13px;background:white}
th,td{border:1px solid #e2e8f0;padding:8px;vertical-align:top}
th{background:#334155;color:white}
.badge{display:inline-block;border-radius:999px;padding:4px 8px;font-weight:700;font-size:12px}
.close{background:#dcfce7;color:#166534}
.intermediate{background:#fef9c3;color:#854d0e}
.distant{background:#e0f2fe;color:#075985}
.note{background:#eff6ff;border-left:5px solid #2563eb;padding:12px;border-radius:10px}
img{max-width:100%;height:auto;border-radius:14px}
"""

    out = [
        "<!doctype html>",
        "<html>",
        "<head>",
        "<meta charset='utf-8'>",
        "<title>Pairwise SNP Distance and Cluster Summary</title>",
        "<style>",
        css,
        "</style>",
        "</head>",
        "<body>",

        "<div class='card'>",
        "<h1>Pairwise SNP Distance and Cluster Summary</h1>",
        "<p>This section summarizes pairwise SNP distances calculated directly from the MTBC core genome alignment after excluding reference/non-sample sequences.</p>",
        f"<div class='note'>Interpretation thresholds used: &le;{likely_threshold} SNPs = genomically close and should be reviewed with epidemiological linkage data; {likely_threshold + 1}-{possible_threshold} SNPs = intermediate SNP distance requiring epidemiological metadata review; &gt;{possible_threshold} SNPs = not clustered by SNP threshold.</div>",
        "</div>",

        "<div class='card'>",
        "<h2>Pairwise SNP Distance Heatmap</h2>",
        "<p>Reference/non-sample sequences are excluded. Lower SNP distances indicate closer genomic relatedness.</p>",
        "<img src='pairwise_snp_heatmap.png' alt='Pairwise SNP Distance Heatmap'>",
        "</div>",

        "<div class='card'>",
        "<h2>Pairwise SNP distance interpretation</h2>",
        "<table>",
        "<thead>",
        "<tr>",
        "<th>Sample 1</th>",
        "<th>Sample 2</th>",
        "<th>SNP distance</th>",
        "<th>Comparable sites</th>",
        "<th>Interpretation</th>",
        "</tr>",
        "</thead>",
        "<tbody>"
    ]

    for p in pairs:
        cls = html.escape(p["cluster_class"])

        out.append("<tr>")
        out.append(f"<td>{html.escape(p['sample1'])}</td>")
        out.append(f"<td>{html.escape(p['sample2'])}</td>")
        out.append(f"<td>{html.escape(str(p['snp_distance']))}</td>")
        out.append(f"<td>{html.escape(str(p['comparable_sites']))}</td>")
        out.append(
            f"<td><span class='badge {cls}'>{html.escape(p['interpretation'])}</span></td>"
        )
        out.append("</tr>")

    out.extend([
        "</tbody>",
        "</table>",
        "</div>",

        "<div class='card'>",
        "<h2>Genomically close sample pairs requiring epidemiological review</h2>",
        "<table>",
        "<thead>",
        "<tr>",
        "<th>Cluster ID</th>",
        "<th>Sample 1</th>",
        "<th>Sample 2</th>",
        "<th>SNP distance</th>",
        "<th>Interpretation</th>",
        "</tr>",
        "</thead>",
        "<tbody>"
    ])

    if close_pairs:
        for idx, p in enumerate(close_pairs, start=1):
            cls = html.escape(p["cluster_class"])

            out.append("<tr>")
            out.append(f"<td>Cluster_{idx}</td>")
            out.append(f"<td>{html.escape(p['sample1'])}</td>")
            out.append(f"<td>{html.escape(p['sample2'])}</td>")
            out.append(f"<td>{html.escape(str(p['snp_distance']))}</td>")
            out.append(
                f"<td><span class='badge {cls}'>{html.escape(p['interpretation'])}</span></td>"
            )
            out.append("</tr>")
    else:
        out.append(
            "<tr><td colspan='5'>No genomically close or intermediate-distance sample pairs were detected using the configured SNP thresholds.</td></tr>"
        )

    out.extend([
        "</tbody>",
        "</table>",
        "</div>",

        "</body>",
        "</html>"
    ])

    html_file.write_text(
        "\n".join(out),
        encoding="utf-8"
    )

    status_file.write_text(
        f"SNP distance clustering completed successfully. {len(excluded)} reference/non-sample sequence(s) excluded. {len(samples)} sample sequence(s) analyzed.\n",
        encoding="utf-8"
    )

except Exception as exc:
    traceback.print_exc()

    write_fallback_outputs(
        "SNP distance clustering encountered an error but fallback outputs were generated: " + str(exc)
    )

    raise SystemExit(0)

PY
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 100 HDD"
    timeout: "24 hours"
    continueOnReturnCode: [0]
  }

  output {
    File pairwise_snp_distance_matrix = "snp_distance/pairwise_snp_distance_matrix.tsv"
    File pairwise_snp_distance_pairs = "snp_distance/pairwise_snp_distance_pairs.tsv"
    File snp_cluster_summary = "snp_distance/snp_cluster_summary.tsv"
    File snp_distance_cluster_html = "snp_distance/snp_distance_cluster_report.html"
    File pairwise_snp_heatmap = "snp_distance/pairwise_snp_heatmap.png"
    File snp_distance_command_log = "logs/snp_distance.command.log"
    File snp_distance_status = "snp_distance/snp_distance_status.txt"
    File excluded_reference_or_non_sample_sequences = "snp_distance/excluded_reference_or_non_sample_sequences.tsv"
  }
}
task GUBBINS_RECOMBINATION {
  input {
    String docker_image = "staphb/gubbins:3.4.1"
    File core_full_alignment
    Int cpu = 16
    Int memory_gb = 64
  }

  command <<<
    set -uo pipefail
    mkdir -p gubbins logs

    if command -v run_gubbins.py >/dev/null 2>&1; then
      GUBBINS_BIN="$(command -v run_gubbins.py)"
    elif [ -x /usr/local/bin/run_gubbins.py ]; then
      GUBBINS_BIN="/usr/local/bin/run_gubbins.py"
    elif [ -x /opt/conda/bin/run_gubbins.py ]; then
      GUBBINS_BIN="/opt/conda/bin/run_gubbins.py"
    else
      echo "ERROR: run_gubbins.py executable not found." >&2
      exit 127
    fi

    export TMPDIR=/tmp
    export TMP=/tmp
    export TEMP=/tmp

    workdir="/tmp/gubbins_work_${RANDOM}_${RANDOM}"
    mkdir -p "$workdir"

    cp "~{core_full_alignment}" "$workdir/core.full.aln"

    if [ ! -s "$workdir/core.full.aln" ]; then
      echo "ERROR: core.full.aln is missing or empty." >&2
      exit 1
    fi

    echo "Using Gubbins: ${GUBBINS_BIN}" > logs/gubbins.command.log
    echo "Input alignment: core.full.aln" >> logs/gubbins.command.log

    cd "$workdir"

    if "$GUBBINS_BIN" \
      --threads ~{cpu} \
      --prefix gubbins \
      core.full.aln \
      > "$OLDPWD/logs/gubbins.run.log" 2>&1; then

      status="success"

    else
      status="gubbins_failed"
      echo "WARNING: Gubbins failed. Passing original core alignment forward." >> "$OLDPWD/logs/gubbins.command.log"
    fi

    cd "$OLDPWD"

    if [ "$status" = "success" ] && [ -s "$workdir/gubbins.filtered_polymorphic_sites.fasta" ]; then
      cp "$workdir/gubbins.filtered_polymorphic_sites.fasta" gubbins/gubbins.filtered_polymorphic_sites.fasta
    else
      cp "~{core_full_alignment}" gubbins/gubbins.filtered_polymorphic_sites.fasta
    fi

    if [ "$status" = "success" ] && [ -s "$workdir/gubbins.final_tree.tre" ]; then
      cp "$workdir/gubbins.final_tree.tre" gubbins/gubbins.final_tree.tre
    else
      echo "(Gubbins_failed_or_skipped:0.0);" > gubbins/gubbins.final_tree.tre
    fi

    if [ "$status" = "success" ] && [ -s "$workdir/gubbins.recombination_predictions.gff" ]; then
      cp "$workdir/gubbins.recombination_predictions.gff" gubbins/gubbins.recombination_predictions.gff
    else
      cat > gubbins/gubbins.recombination_predictions.gff <<'EOF'
##gff-version 3
# Gubbins failed or was skipped; no recombination predictions available.
EOF
    fi

    echo "$status" > gubbins/gubbins_status.txt

    rm -rf "$workdir"
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 400 HDD"
    timeout: "120 hours"
  }

  output {
    File filtered_alignment = "gubbins/gubbins.filtered_polymorphic_sites.fasta"
    File gubbins_tree = "gubbins/gubbins.final_tree.tre"
    File recombination_predictions = "gubbins/gubbins.recombination_predictions.gff"
    File gubbins_status = "gubbins/gubbins_status.txt"
    File gubbins_command_log = "logs/gubbins.command.log"
    File gubbins_run_log = "logs/gubbins.run.log"
  }
}
task IQTREE2_PHYLOGENY {
  input {
    String docker_image = "gmboowa/iqtree2-python:2.3.4"
    File alignment
    String model = "GTR+G"
    Int bootstrap_replicates = 1000
    Int cpu = 16
    Int memory_gb = 64
    Boolean midpoint_root_tree = true

    # Samples with missing/gap/ambiguous content >= this threshold
    # will be excluded from IQ-TREE phylogenetic inference only.
    # They remain part of the wider workflow and can be reported downstream.
    Float max_missing_fraction_for_tree = 0.50

    # Minimum number of non-reference samples required to run IQ-TREE.
    Int min_non_reference_samples_for_tree = 3
  }

  command <<<
    set -uo pipefail
    mkdir -p iqtree logs

    ###########################################################################
    # Guaranteed placeholder outputs
    #
    # Cromwell fails output collection if a declared File output does not exist.
    # These placeholders are overwritten later when filtering and IQ-TREE run
    # successfully, but they protect the task from missing-output failures if
    # filtering exits early or IQ-TREE is skipped.
    ###########################################################################

    echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note" > iqtree/excluded_from_iqtree.tsv
    echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction" > iqtree/included_in_iqtree.tsv
    echo "Alignment filtering has not yet completed." > iqtree/alignment_filtering_summary.txt
    echo "IQ-TREE log not available." > iqtree/iqtree.log
    echo "IQ-TREE report not available." > iqtree/iqtree.report
    echo "" > iqtree/support_labels.txt
    echo "not_started" > iqtree/iqtree_status.txt
    echo "(IQTREE_not_started:0.0);" > iqtree/final.treefile
    echo "IQ-TREE was not executed or produced no run log." > logs/iqtree.run.log
    echo "IQ-TREE command log initialized." > logs/iqtree.command.log

    if command -v iqtree2 >/dev/null 2>&1; then
      IQTREE_BIN="$(command -v iqtree2)"
    elif command -v iqtree >/dev/null 2>&1; then
      IQTREE_BIN="$(command -v iqtree)"
    elif [ -x /usr/local/bin/iqtree2 ]; then
      IQTREE_BIN="/usr/local/bin/iqtree2"
    else
      echo "ERROR: IQ-TREE executable not found." >&2
      echo "iqtree_executable_not_found" > iqtree/iqtree_status.txt
      exit 127
    fi

    cp "~{alignment}" iqtree/mtbc_core_snp_alignment.original.fasta
    cp "~{alignment}" iqtree/mtbc_core_snp_alignment.fasta

    if [ ! -s iqtree/mtbc_core_snp_alignment.original.fasta ]; then
      echo "ERROR: alignment file is missing or empty." >&2
      echo "alignment_missing_or_empty" > iqtree/iqtree_status.txt
      echo "ERROR: alignment file is missing or empty." > iqtree/alignment_filtering_summary.txt
      exit 1
    fi

    echo "Using IQ-TREE: ${IQTREE_BIN}" > logs/iqtree.command.log
    echo "Model: ~{model}" >> logs/iqtree.command.log
    echo "Bootstrap replicates: ~{bootstrap_replicates}" >> logs/iqtree.command.log
    echo "CPU threads: ~{cpu}" >> logs/iqtree.command.log
    echo "Maximum missing fraction allowed for tree: ~{max_missing_fraction_for_tree}" >> logs/iqtree.command.log
    echo "Minimum non-reference samples required for tree: ~{min_non_reference_samples_for_tree}" >> logs/iqtree.command.log
    echo "" >> logs/iqtree.command.log

    ###########################################################################
    # Pre-IQ-TREE filtering
    #
    # This does NOT remove samples from the whole workflow.
    # It only removes samples from IQ-TREE phylogenetic inference when their
    # core-SNP alignment sequence is too poor for tree reconstruction.
    #
    # Outputs:
    #   iqtree/mtbc_core_snp_alignment.original.fasta = original alignment
    #   iqtree/mtbc_core_snp_alignment.fasta          = filtered alignment for IQ-TREE
    #   iqtree/excluded_from_iqtree.tsv               = samples excluded from tree only
    #   iqtree/included_in_iqtree.tsv                 = samples retained for tree
    #   iqtree/alignment_filtering_summary.txt        = filtering summary
    ###########################################################################

    python3 <<'PY'
from collections import OrderedDict
import sys

input_fasta = "iqtree/mtbc_core_snp_alignment.original.fasta"
clean_fasta = "iqtree/mtbc_core_snp_alignment.fasta"
excluded_tsv = "iqtree/excluded_from_iqtree.tsv"
included_tsv = "iqtree/included_in_iqtree.tsv"
summary_txt = "iqtree/alignment_filtering_summary.txt"

max_missing_fraction = float("~{max_missing_fraction_for_tree}")
min_non_reference_samples = int("~{min_non_reference_samples_for_tree}")

records = OrderedDict()
current_name = None
current_seq = []

with open(input_fasta, "r", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        line = line.rstrip("\n")
        if not line:
            continue

        if line.startswith(">"):
            if current_name is not None:
                records[current_name] = "".join(current_seq)

            current_name = line[1:].strip().split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())

if current_name is not None:
    records[current_name] = "".join(current_seq)

def write_empty_failure_outputs(message, exit_code):
    with open(excluded_tsv, "w") as out:
        out.write("sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note\n")

    with open(included_tsv, "w") as out:
        out.write("sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\n")

    with open(summary_txt, "w") as out:
        out.write(message.rstrip() + "\n")

    with open(clean_fasta, "w") as out:
        out.write("")

    sys.exit(exit_code)

if not records:
    write_empty_failure_outputs(
        "ERROR: No FASTA records were found in the alignment.",
        10
    )

lengths = {name: len(seq) for name, seq in records.items()}
unique_lengths = sorted(set(lengths.values()))

if len(unique_lengths) != 1:
    with open(excluded_tsv, "w") as out:
        out.write("sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note\n")

    with open(included_tsv, "w") as out:
        out.write("sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\n")

    with open(summary_txt, "w") as out:
        out.write("ERROR: Alignment records do not all have the same length.\n")
        out.write("sample\tsequence_length\n")
        for name, length in lengths.items():
            out.write(f"{name}\t{length}\n")

    with open(clean_fasta, "w") as out:
        out.write("")

    sys.exit(11)

alignment_length = unique_lengths[0]

included = OrderedDict()
excluded = []

threshold_percent_text = f"{max_missing_fraction * 100:.0f}%"

for name, seq in records.items():
    seq_upper = seq.upper()

    acgt_count = sum(1 for base in seq_upper if base in {"A", "C", "G", "T"})
    missing_count = sum(1 for base in seq_upper if base not in {"A", "C", "G", "T"})
    missing_fraction = missing_count / alignment_length if alignment_length > 0 else 1.0
    missing_percent_text = f"{missing_fraction * 100:.2f}%"

    reason = ""
    exclusion_note = ""

    if alignment_length == 0:
        reason = "empty_sequence"
        exclusion_note = (
            f"{name} was excluded from IQ-TREE because its core-SNP alignment "
            f"sequence was empty."
        )

    elif acgt_count == 0:
        reason = "no_usable_acgt_bases"
        exclusion_note = (
            f"{name} was excluded from IQ-TREE because it had "
            f"{missing_percent_text} missing/ambiguous/gap content in the "
            f"core-SNP alignment and no usable ACGT bases."
        )

    elif missing_fraction >= max_missing_fraction:
        reason = f"missing_fraction_ge_{max_missing_fraction}"
        exclusion_note = (
            f"{name} was excluded from IQ-TREE because "
            f"{missing_percent_text} of its core-SNP alignment was "
            f"missing/ambiguous/gap content, exceeding the maximum allowed "
            f"threshold of {threshold_percent_text}."
        )

    if reason:
        excluded.append({
            "sample": name,
            "alignment_length": alignment_length,
            "acgt_count": acgt_count,
            "missing_count": missing_count,
            "missing_fraction": missing_fraction,
            "threshold": max_missing_fraction,
            "reason": reason,
            "exclusion_note": exclusion_note
        })
    else:
        included[name] = seq

non_reference_included = [
    name for name in included
    if name.lower() not in {"reference", "ref", "h37rv", "h37rv_reference"}
]

with open(excluded_tsv, "w") as out:
    out.write("sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note\n")
    for row in excluded:
        out.write(
            f"{row['sample']}\t"
            f"{row['alignment_length']}\t"
            f"{row['acgt_count']}\t"
            f"{row['missing_count']}\t"
            f"{row['missing_fraction']:.6f}\t"
            f"{row['threshold']:.6f}\t"
            f"{row['reason']}\t"
            f"{row['exclusion_note']}\n"
        )

with open(included_tsv, "w") as out:
    out.write("sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\n")
    for name, seq in included.items():
        seq_upper = seq.upper()
        acgt_count = sum(1 for base in seq_upper if base in {"A", "C", "G", "T"})
        missing_count = sum(1 for base in seq_upper if base not in {"A", "C", "G", "T"})
        missing_fraction = missing_count / alignment_length if alignment_length > 0 else 1.0

        out.write(
            f"{name}\t"
            f"{alignment_length}\t"
            f"{acgt_count}\t"
            f"{missing_count}\t"
            f"{missing_fraction:.6f}\n"
        )

with open(clean_fasta, "w") as out:
    for name, seq in included.items():
        out.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

with open(summary_txt, "w") as out:
    out.write(f"Original sequences: {len(records)}\n")
    out.write(f"Included sequences: {len(included)}\n")
    out.write(f"Excluded sequences: {len(excluded)}\n")
    out.write(f"Included non-reference samples: {len(non_reference_included)}\n")
    out.write(f"Alignment length: {alignment_length}\n")
    out.write(f"Maximum allowed missing fraction: {max_missing_fraction}\n")
    out.write("\n")

    if excluded:
        out.write("Excluded samples from IQ-TREE only:\n")
        for row in excluded:
            out.write(
                f"- {row['sample']}: "
                f"{row['missing_fraction']:.2%} missing/ambiguous/gap content; "
                f"{row['reason']}; "
                f"{row['exclusion_note']}\n"
            )
    else:
        out.write("Excluded samples from IQ-TREE only: none\n")

    out.write("\n")
    out.write(
        "Note: Excluded samples are not removed from the wider workflow. "
        "They are excluded only from IQ-TREE phylogenetic inference because "
        "their core-SNP alignment sequence did not meet the minimum quality "
        "requirements for tree reconstruction.\n"
    )

if len(non_reference_included) < min_non_reference_samples:
    sys.exit(12)

sys.exit(0)
PY

    filter_exit_code=$?

    echo "Alignment filtering exit code: ${filter_exit_code}" >> logs/iqtree.command.log

    if [ -s iqtree/alignment_filtering_summary.txt ]; then
      cat iqtree/alignment_filtering_summary.txt >> logs/iqtree.command.log
    else
      echo "Alignment filtering summary was not generated." >> logs/iqtree.command.log
      echo "Alignment filtering summary was not generated." > iqtree/alignment_filtering_summary.txt
    fi

    echo "" >> logs/iqtree.command.log

    if [ ! -f iqtree/excluded_from_iqtree.tsv ]; then
      echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note" > iqtree/excluded_from_iqtree.tsv
    fi

    if [ ! -f iqtree/included_in_iqtree.tsv ]; then
      echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction" > iqtree/included_in_iqtree.tsv
    fi

    if [ ! -f iqtree/mtbc_core_snp_alignment.fasta ]; then
      cp iqtree/mtbc_core_snp_alignment.original.fasta iqtree/mtbc_core_snp_alignment.fasta || true
    fi

    if [ "${filter_exit_code}" -eq 10 ] || [ "${filter_exit_code}" -eq 11 ]; then
      status="alignment_filtering_failed"
      echo "ERROR: Alignment filtering failed before IQ-TREE." >> logs/iqtree.command.log

    elif [ "${filter_exit_code}" -eq 12 ]; then
      status="too_few_samples_after_filtering"
      echo "WARNING: Too few valid non-reference samples remained after filtering. Generating fallback tree." >> logs/iqtree.command.log

    elif [ "${filter_exit_code}" -eq 0 ]; then
      echo "Proceeding with filtered alignment:" >> logs/iqtree.command.log
      echo "  iqtree/mtbc_core_snp_alignment.fasta" >> logs/iqtree.command.log
      echo "" >> logs/iqtree.command.log

      if "$IQTREE_BIN" \
        -s iqtree/mtbc_core_snp_alignment.fasta \
        -m ~{model} \
        -B ~{bootstrap_replicates} \
        -alrt ~{bootstrap_replicates} \
        -bnni \
        -nt ~{cpu} \
        -pre iqtree/MTBC_core_SNP_phylogeny \
        >> logs/iqtree.run.log 2>&1; then

        status="success"

      else
        status="iqtree_failed_after_filtering"
        echo "WARNING: IQ-TREE failed even after problematic-sample filtering. Generating fallback tree." >> logs/iqtree.command.log
      fi

    else
      status="unknown_filtering_error"
      echo "ERROR: Unknown filtering error before IQ-TREE. Generating fallback tree." >> logs/iqtree.command.log
    fi

    if [ "$status" = "success" ] && [ -s iqtree/MTBC_core_SNP_phylogeny.treefile ]; then
      cp iqtree/MTBC_core_SNP_phylogeny.treefile iqtree/final.treefile
    else
      echo "(IQTREE_failed:0.0);" > iqtree/final.treefile
    fi

    if [ -s iqtree/MTBC_core_SNP_phylogeny.log ]; then
      cp iqtree/MTBC_core_SNP_phylogeny.log iqtree/iqtree.log
    else
      echo "IQ-TREE log not available." > iqtree/iqtree.log
    fi

    if [ -s iqtree/MTBC_core_SNP_phylogeny.iqtree ]; then
      cp iqtree/MTBC_core_SNP_phylogeny.iqtree iqtree/iqtree.report
    else
      echo "IQ-TREE report not available." > iqtree/iqtree.report
    fi

    if [ ! -s logs/iqtree.run.log ]; then
      echo "IQ-TREE was not executed or produced no run log." > logs/iqtree.run.log
    fi

    grep -oE '\)[0-9]+(\.[0-9]+)?(/[0-9]+(\.[0-9]+)?)?:' iqtree/final.treefile > iqtree/support_labels.txt || true

    if [ ! -f iqtree/support_labels.txt ]; then
      echo "" > iqtree/support_labels.txt
    fi

    echo "$status" > iqtree/iqtree_status.txt

    if [ ! -f iqtree/excluded_from_iqtree.tsv ]; then
      echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note" > iqtree/excluded_from_iqtree.tsv
    fi

    if [ ! -f iqtree/included_in_iqtree.tsv ]; then
      echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction" > iqtree/included_in_iqtree.tsv
    fi

    if [ ! -f iqtree/alignment_filtering_summary.txt ]; then
      echo "Alignment filtering summary not available." > iqtree/alignment_filtering_summary.txt
    fi

    if [ ! -f iqtree/mtbc_core_snp_alignment.original.fasta ]; then
      echo "" > iqtree/mtbc_core_snp_alignment.original.fasta
    fi

    if [ ! -f iqtree/mtbc_core_snp_alignment.fasta ]; then
      echo "" > iqtree/mtbc_core_snp_alignment.fasta
    fi

    if [ ! -f iqtree/final.treefile ]; then
      echo "(IQTREE_failed:0.0);" > iqtree/final.treefile
    fi

    if [ ! -f iqtree/iqtree.log ]; then
      echo "IQ-TREE log not available." > iqtree/iqtree.log
    fi

    if [ ! -f iqtree/iqtree.report ]; then
      echo "IQ-TREE report not available." > iqtree/iqtree.report
    fi

    if [ ! -f iqtree/iqtree_status.txt ]; then
      echo "unknown_status" > iqtree/iqtree_status.txt
    fi

    if [ ! -f logs/iqtree.command.log ]; then
      echo "IQ-TREE command log not available." > logs/iqtree.command.log
    fi

    if [ ! -f logs/iqtree.run.log ]; then
      echo "IQ-TREE run log not available." > logs/iqtree.run.log
    fi
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 400 HDD"
    timeout: "240 hours"
  }

  output {
    File final_tree = "iqtree/final.treefile"
    File iqtree_log = "iqtree/iqtree.log"
    File iqtree_report = "iqtree/iqtree.report"
    File support_labels = "iqtree/support_labels.txt"
    File iqtree_status = "iqtree/iqtree_status.txt"
    File iqtree_command_log = "logs/iqtree.command.log"
    File iqtree_run_log = "logs/iqtree.run.log"

    File original_alignment = "iqtree/mtbc_core_snp_alignment.original.fasta"
    File filtered_alignment = "iqtree/mtbc_core_snp_alignment.fasta"
    File excluded_from_iqtree = "iqtree/excluded_from_iqtree.tsv"
    File included_in_iqtree = "iqtree/included_in_iqtree.tsv"
    File alignment_filtering_summary = "iqtree/alignment_filtering_summary.txt"
  }
}
task TREE_VISUALIZATION {
  input {
    File input_tree
    File? tbprofiler_summary_tsv
    File? resistance_profile_summary_tsv
    File? iqtree_excluded_samples_tsv
    String docker_image = "gmboowa/ete3-render:1.18"
    Int width = 2400
    Int height = 1600
    String image_format = "png"
    String title = "MTBC Core-SNP Phylogenetic Tree"
  }

  command <<<
    set -uo pipefail

    mkdir -p tree_visualization
    export QT_QPA_PLATFORM=offscreen
    export MPLBACKEND=Agg

    python3 - <<'PY'
from pathlib import Path
import csv
import re
import base64
import traceback

outdir = Path("tree_visualization")
outdir.mkdir(exist_ok=True)

image_format = "~{image_format}".lower().strip()
if image_format not in {"png", "svg", "pdf"}:
    image_format = "png"

out_img = outdir / f"phylogenetic_tree.{image_format}"
cleaned_tree = outdir / "phylogenetic_tree.cleaned.nwk"
log = outdir / "render.log"
excluded_copy = outdir / "excluded_from_iqtree.tsv"

tree_input = "~{input_tree}"
tbprofiler_summary_path = "~{if defined(tbprofiler_summary_tsv) then tbprofiler_summary_tsv else ""}"
resistance_profile_summary_path = "~{if defined(resistance_profile_summary_tsv) then resistance_profile_summary_tsv else ""}"
iqtree_excluded_samples_path = "~{if defined(iqtree_excluded_samples_tsv) then iqtree_excluded_samples_tsv else ""}"
requested_width = int("~{width}")
requested_height = int("~{height}")
title = "~{title}"

log.write_text("TREE_VISUALIZATION started\n")

if iqtree_excluded_samples_path and Path(iqtree_excluded_samples_path).exists():
    try:
        excluded_copy.write_text(
            Path(iqtree_excluded_samples_path).read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8"
        )
        with open(log, "a") as fh:
            fh.write(f"Copied IQ-TREE exclusion file: {iqtree_excluded_samples_path}\n")
    except Exception as exc:
        with open(log, "a") as fh:
            fh.write(f"WARNING: Could not copy IQ-TREE exclusion file: {repr(exc)}\n")
else:
    excluded_copy.write_text(
        "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\treason\n",
        encoding="utf-8"
    )

RESISTANCE_COLORS = {
    # Deliberately unique colors for every final resistance category.
    # These colors are used consistently in the tree, legends, badges, and
    # surveillance tables so categories are not visually collapsed.
    "Sensitive": "#2A9D8F",                               # teal
    "Hr-TB": "#1D9BF0",                                  # blue
    "RR-TB": "#F59E0B",                                  # amber
    "MDR-TB": "#E11D48",                                 # rose
    "MDR/RR-TB": "#F97316",                              # orange
    "Pre-XDR-TB": "#EF4444",                             # red
    "XDR-TB": "#6D28D9",                                 # purple
    "Monoresistance": "#FACC15",                         # yellow
    "Polyresistance": "#06B6D4",                         # cyan
    "Other drug resistance": "#2563EB",                   # royal blue
    "Resistance not determined by TB-Profiler": "#94A3B8", # slate gray
    "Unknown": "#6B7280"                                  # dark gray
}

CANONICAL_RESISTANCE_ORDER = [
    "Sensitive",
    "Hr-TB",
    "RR-TB",
    "MDR-TB",
    "MDR/RR-TB",
    "Pre-XDR-TB",
    "XDR-TB",
    "Monoresistance",
    "Polyresistance",
    "Other drug resistance",
    "Resistance not determined by TB-Profiler",
    "Unknown"
]

def write_placeholder_image(message):
    msg = str(message).replace("<", "&lt;").replace(">", "&gt;")
    svg_text = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="500" viewBox="0 0 1200 500">
<rect width="100%" height="100%" fill="#ffffff"/>
<rect x="40" y="40" width="1120" height="420" rx="24" fill="#f8fafc" stroke="#cbd5e1" stroke-width="3"/>
<text x="600" y="170" text-anchor="middle" font-family="Arial" font-size="34" font-weight="700" fill="#0f172a">Tree visualization was not rendered</text>
<text x="600" y="225" text-anchor="middle" font-family="Arial" font-size="22" fill="#475569">The Newick tree was preserved for downstream reporting.</text>
<text x="600" y="285" text-anchor="middle" font-family="Arial" font-size="18" fill="#b91c1c">{msg[:180]}</text>
</svg>'''

    if image_format == "svg":
        out_img.write_text(svg_text, encoding="utf-8")
    elif image_format == "pdf":
        out_img.write_text("Tree visualization was not rendered. Newick tree preserved.\n" + str(message) + "\n", encoding="utf-8")
    else:
        png_b64 = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/p9sAAAAASUVORK5CYII="
        out_img.write_bytes(base64.b64decode(png_b64))

def normalize_tree_leaf_name(name):
    s = str(name or "").strip().strip("'").strip('"')
    s = s.split("/")[-1]
    s = s.split("\\")[-1]

    suffixes = [
        ".consensus.subs",
        ".consensus",
        ".aligned",
        ".snps",
        ".subs",
        ".tree",
        ".nwk",
        ".newick",
        ".fa",
        ".fasta",
        ".fna",
        ".bam",
        ".sam",
        ".vcf",
        ".gz"
    ]

    changed = True
    while changed:
        changed = False
        for suffix in suffixes:
            if s.lower().endswith(suffix.lower()):
                s = s[:-len(suffix)]
                changed = True

    s = re.sub(r"_R[12](_001)?$", "", s)
    s = re.sub(r"_[12]$", "", s)
    s = re.sub(r"\.R[12](_001)?$", "", s)
    s = re.sub(r"\.[12]$", "", s)
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"^\d{3,}[_-](?=(?:[DES]RR|SAM[END]?|ERS|SRS|DRS|SRX|ERX|DRX|[A-Za-z]))", "", s)
    return s.strip()

def normalize_metadata_sample_id(name):
    return normalize_tree_leaf_name(name)

def is_reference_tip(name):
    normalized = normalize_tree_leaf_name(name).lower()
    reference_names = {
        "reference",
        "ref",
        "h37rv",
        "h37rv_siena",
        "h37rvsiena",
        "nc_000962",
        "nc_000962.3",
        "mycobacterium_tuberculosis_h37rv",
        "mycobacterium_tuberculosis_h37rv_siena"
    }
    return normalized in reference_names

def normalize_profile_label(value):
    raw = str(value or "").strip()
    low = raw.lower().replace("_", "-")

    if low in {"", "unknown", "not reported", "na", "n/a", "none"}:
        return "Unknown"

    if low in {"sensitive", "susceptible", "no resistance detected by tb-profiler", "no resistance detected"}:
        return "Sensitive"

    if "resistance not determined" in low:
        return "Resistance not determined by TB-Profiler"

    if "xdr" in low and "pre" not in low:
        return "XDR-TB"

    if "pre-xdr" in low or "pre xdr" in low or "prexdr" in low:
        return "Pre-XDR-TB"

    # Preserve the umbrella only when no drug-level calls are available.
    # Drug-level classification elsewhere converts rif-only to RR-TB and rif+INH to MDR-TB.
    if re.search(r"\bmdr\s*/\s*rr[- ]?tb\b", low) or re.search(r"\bmdr[- ]?rr[- ]?tb\b", low):
        return "MDR/RR-TB"

    if "mdr" in low:
        return "MDR-TB"

    if re.search(r"\brr[- ]?tb\b", low):
        return "RR-TB"

    if re.search(r"\bhr[- ]?tb\b", low) or "isoniazid-resistant" in low:
        return "Hr-TB"

    if "mono" in low:
        return "Monoresistance"

    if "poly" in low:
        return "Polyresistance"

    if "other drug resistance" in low:
        return "Other drug resistance"

    return "Unknown"


def profile_color(profile):
    return RESISTANCE_COLORS.get(profile, RESISTANCE_COLORS["Unknown"])

def normalize_drug_name(x):
    s = str(x or "").strip().lower()
    s = re.sub(r"[_\-]+", " ", s)
    s = re.sub(r"\s+", " ", s)

    aliases = {
        "inh": "isoniazid",
        "isoniazid": "isoniazid",
        "h": "isoniazid",
        "rif": "rifampicin",
        "rmp": "rifampicin",
        "rifampin": "rifampicin",
        "rifampicin": "rifampicin",
        "r": "rifampicin",
        "pza": "pyrazinamide",
        "pyrazinamide": "pyrazinamide",
        "z": "pyrazinamide",
        "emb": "ethambutol",
        "ethambutol": "ethambutol",
        "e": "ethambutol",
        "sm": "streptomycin",
        "str": "streptomycin",
        "streptomycin": "streptomycin",
        "s": "streptomycin",
        "levo": "levofloxacin",
        "levofloxacin": "levofloxacin",
        "lfx": "levofloxacin",
        "moxi": "moxifloxacin",
        "moxifloxacin": "moxifloxacin",
        "mfx": "moxifloxacin",
        "ofx": "ofloxacin",
        "ofloxacin": "ofloxacin",
        "amikacin": "amikacin",
        "amk": "amikacin",
        "kanamycin": "kanamycin",
        "kan": "kanamycin",
        "capreomycin": "capreomycin",
        "cap": "capreomycin",
        "bedaquiline": "bedaquiline",
        "bdq": "bedaquiline",
        "linezolid": "linezolid",
        "lzd": "linezolid",
        "clofazimine": "clofazimine",
        "cfz": "clofazimine",
        "ethionamide": "ethionamide",
        "eto": "ethionamide",
        "prothionamide": "prothionamide",
        "pto": "prothionamide",
        "cycloserine": "cycloserine",
        "cs": "cycloserine",
        "para aminosalicylic acid": "para-aminosalicylic acid",
        "pas": "para-aminosalicylic acid",
        "delamanid": "delamanid",
        "dlm": "delamanid",
        "pretomanid": "pretomanid",
        "pa": "pretomanid"
    }

    return aliases.get(s, s)

def split_drug_tokens(value):
    if not value:
        return []

    txt = str(value or "")
    txt = txt.replace(";", ",").replace("|", ",").replace("/", ",")
    txt = re.sub(r"\band\b", ",", txt, flags=re.IGNORECASE)

    out = []

    for part in txt.split(","):
        p = normalize_drug_name(part)
        if p and p not in {"none", "none reported", "not reported", "susceptible", "sensitive", "unknown", "na", "n/a"}:
            out.append(p)

    return out

def classify_resistance_from_drugs(dr_type, resistant_drugs):
    """
    Final tree-label resistance classifier.

    Important fix:
    TB-Profiler may use the umbrella string MDR/RR-TB. That umbrella must not
    be displayed as the final isolate category when drug-level calls are
    available. Drug-level calls are authoritative for separating:
      - rifampicin only              -> RR-TB
      - rifampicin + isoniazid      -> MDR-TB
      - rifampicin + FQ             -> Pre-XDR-TB
      - rifampicin + FQ + BDQ/LZD   -> XDR-TB
    The raw dr_type string is used only when no parseable resistant_drugs field
    is available.
    """
    drugs = set(split_drug_tokens(resistant_drugs))

    if drugs:
        isoniazid = "isoniazid" in drugs
        rifampicin = "rifampicin" in drugs

        fluoroquinolones = {
            "levofloxacin",
            "moxifloxacin",
            "ofloxacin",
            "gatifloxacin",
            "ciprofloxacin"
        }

        group_a_additional = {
            "bedaquiline",
            "linezolid"
        }

        has_fq = bool(drugs.intersection(fluoroquinolones))
        has_group_a_additional = bool(drugs.intersection(group_a_additional))

        if rifampicin and has_fq and has_group_a_additional:
            return "XDR-TB", RESISTANCE_COLORS["XDR-TB"]

        if rifampicin and has_fq:
            return "Pre-XDR-TB", RESISTANCE_COLORS["Pre-XDR-TB"]

        if rifampicin and isoniazid:
            return "MDR-TB", RESISTANCE_COLORS["MDR-TB"]

        if rifampicin and not isoniazid:
            return "RR-TB", RESISTANCE_COLORS["RR-TB"]

        if isoniazid and not rifampicin:
            return "Hr-TB", RESISTANCE_COLORS["Hr-TB"]

        if len(drugs) == 1:
            return "Monoresistance", RESISTANCE_COLORS["Monoresistance"]

        return "Polyresistance", RESISTANCE_COLORS["Polyresistance"]

    profile = normalize_profile_label(dr_type)
    return profile, profile_color(profile)


def clean_single_bootstrap_value(raw_value):
    raw = str(raw_value or "").strip().strip("'\"")
    if not raw:
        return ""

    numeric_tokens = re.findall(r"\d+(?:\.\d+)?", raw)
    if not numeric_tokens:
        return ""

    support_value = numeric_tokens[-1] if len(numeric_tokens) >= 2 else numeric_tokens[0]

    if support_value.isdigit():
        if support_value == "100100":
            support_value = "100"
        elif len(support_value) >= 4 and support_value.endswith("100"):
            support_value = "100"
        elif len(support_value) == 4:
            support_value = support_value[2:]

    try:
        value = float(support_value)

        if 0 < value <= 1:
            value *= 100

        if value <= 0 or value > 100:
            return ""

        return str(int(round(value)))

    except Exception:
        return ""

try:
    if not tree_input:
        raise ValueError("input_tree not provided")

    tree_path = Path(tree_input)

    if not tree_path.exists() or tree_path.stat().st_size == 0:
        raise FileNotFoundError(f"Tree missing or empty: {tree_path}")

    raw_newick = tree_path.read_text(errors="replace").strip()

    if not raw_newick:
        raise ValueError("Tree file is empty")

    cleaned_tree.write_text(raw_newick + ("\n" if not raw_newick.endswith("\n") else ""), encoding="utf-8")

    if "IQTREE_failed" in raw_newick or "alignment_filtering_failed" in raw_newick or raw_newick.strip() in {"();", "(IQTREE_failed:0.0);"}:
        raise ValueError("IQ-TREE did not produce a valid phylogenetic tree. A fallback Newick tree was detected.")

    metadata = {}
    resistance_profile_map = {}
    resistance_profile_source = "none"

    if resistance_profile_summary_path and Path(resistance_profile_summary_path).exists():
        with open(resistance_profile_summary_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            for r in reader:
                sample_raw = (
                    r.get("sample_id") or
                    r.get("sample") or
                    r.get("Sample ID") or
                    ""
                ).strip()

                if not sample_raw:
                    continue

                sample = normalize_metadata_sample_id(sample_raw)

                if not sample:
                    continue

                profile_raw = (
                    r.get("resistance_profile") or
                    r.get("dr_type") or
                    r.get("Resistance profile") or
                    ""
                )

                resistant_drugs = (
                    r.get("resistant_drugs") or
                    r.get("Predicted resistant drugs") or
                    r.get("drug_resistance") or
                    ""
                )

                category, color = classify_resistance_from_drugs(profile_raw, resistant_drugs)

                resistance_profile_map[sample] = {
                    "category": category,
                    "color": color,
                    "resistant_drugs": resistant_drugs
                }

        resistance_profile_source = "resistance_profile_summary_tsv"

    if tbprofiler_summary_path and Path(tbprofiler_summary_path).exists():
        with open(tbprofiler_summary_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            for r in reader:
                sample_raw = (r.get("sample") or r.get("sample_id") or "").strip()

                if not sample_raw:
                    continue

                sample = normalize_metadata_sample_id(sample_raw)

                if not sample:
                    continue

                main_lineage = (r.get("main_lineage") or "").strip()
                sub_lineage = (r.get("sub_lineage") or "").strip()

                lineage = sub_lineage or main_lineage

                if lineage.lower() in {"not reported", "none", "na", "n/a", "unknown"}:
                    lineage = main_lineage

                if sample in resistance_profile_map:
                    category = resistance_profile_map[sample]["category"]
                    color = resistance_profile_map[sample]["color"]
                else:
                    category, color = classify_resistance_from_drugs(r.get("dr_type"), r.get("resistant_drugs"))

                metadata[sample] = {
                    "lineage": lineage,
                    "category": category,
                    "color": color
                }

    for sample, profile_info in resistance_profile_map.items():
        if sample not in metadata:
            metadata[sample] = {
                "lineage": "Lineage NA",
                "category": profile_info["category"],
                "color": profile_info["color"]
            }

    from ete3 import Tree, TreeStyle, TextFace, NodeStyle

    t = Tree(str(tree_path), format=1)

    original_to_normalized_leaf_names = {}
    removed_refs = []

    for leaf in list(t.get_leaves()):
        original_name = str(leaf.name or "").strip()
        normalized_name = normalize_tree_leaf_name(original_name)
        original_to_normalized_leaf_names[original_name] = normalized_name

        if is_reference_tip(original_name):
            removed_refs.append(original_name)
            leaf.detach()
            continue

        if normalized_name:
            leaf.name = normalized_name

    if len(t.get_leaves()) < 2:
        raise ValueError("Tree has fewer than two non-reference tips after normalizing labels and removing only the true reference tip.")

    try:
        midpoint = t.get_midpoint_outgroup()

        if midpoint:
            t.set_outgroup(midpoint)

    except Exception as e:
        with open(log, "a") as fh:
            fh.write(f"WARNING: midpoint rooting failed: {repr(e)}\n")

    n_leaves = len(t.get_leaves())

    if n_leaves <= 5:
        auto_width = max(requested_width, 3800)
        label_font = 14
        lineage_font = 11
        resistance_font = 11
        bootstrap_font = 10
        tip_node_size = 8
        branch_width = 2
        branch_vertical_margin = 18
        margin_right = 1500
    elif n_leaves <= 10:
        auto_width = max(requested_width, 4200)
        label_font = 13
        lineage_font = 10
        resistance_font = 10
        bootstrap_font = 9
        tip_node_size = 7
        branch_width = 2
        branch_vertical_margin = 14
        margin_right = 1600
    elif n_leaves <= 25:
        auto_width = max(requested_width, 5000)
        label_font = 11
        lineage_font = 9
        resistance_font = 9
        bootstrap_font = 8
        tip_node_size = 6
        branch_width = 2
        branch_vertical_margin = 8
        margin_right = 1700
    elif n_leaves <= 50:
        auto_width = max(requested_width, 5600)
        label_font = 9
        lineage_font = 8
        resistance_font = 8
        bootstrap_font = 7
        tip_node_size = 5
        branch_width = 1
        branch_vertical_margin = 5
        margin_right = 1800
    elif n_leaves <= 100:
        auto_width = max(requested_width, 6400)
        label_font = 8
        lineage_font = 7
        resistance_font = 7
        bootstrap_font = 6
        tip_node_size = 4
        branch_width = 1
        branch_vertical_margin = 3
        margin_right = 1900
    else:
        auto_width = max(requested_width, 7200)
        label_font = 7
        lineage_font = 6
        resistance_font = 6
        bootstrap_font = 5
        tip_node_size = 3
        branch_width = 1
        branch_vertical_margin = 2
        margin_right = 2000

    for node in t.traverse():
        ns = NodeStyle()
        ns["hz_line_width"] = branch_width
        ns["vt_line_width"] = branch_width
        ns["size"] = 0 if not node.is_leaf() else tip_node_size
        node.set_style(ns)

        if not node.is_leaf():
            raw_name = str(getattr(node, "name", "") or "").strip()
            clean_support = clean_single_bootstrap_value(raw_name)

            if not clean_support:
                raw_support = str(getattr(node, "support", "") or "").strip()
                clean_support = clean_single_bootstrap_value(raw_support)

            node.name = ""

            if clean_support:
                node.add_face(TextFace(clean_support, fsize=bootstrap_font, fgcolor="#b91c1c"), column=0, position="branch-top")

    metadata_matched = 0
    metadata_missing = []

    for leaf in t.get_leaves():
        sample = normalize_tree_leaf_name(leaf.name)
        meta = metadata.get(sample, {})

        if meta:
            metadata_matched += 1
        else:
            metadata_missing.append(sample)

        lineage = meta.get("lineage", "") or "Lineage NA"
        category = meta.get("category", "Unknown")
        category = normalize_profile_label(category)
        color = meta.get("color", profile_color(category))

        leaf.name = ""
        leaf.add_face(TextFace(sample, fsize=label_font, fgcolor="#111827"), column=0, position="branch-right")
        leaf.add_face(TextFace(f"  {lineage}", fsize=lineage_font, fgcolor="#2563eb"), column=1, position="branch-right")
        leaf.add_face(TextFace(f"  {category}", fsize=resistance_font, fgcolor=color), column=2, position="branch-right")

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.mode = "r"
    ts.scale = 120
    ts.branch_vertical_margin = branch_vertical_margin
    ts.margin_left = 20
    ts.margin_right = margin_right
    ts.margin_top = 20
    ts.margin_bottom = 20
    ts.title.add_face(TextFace(title, fsize=18, bold=True), column=0)

    t.render(str(out_img), w=auto_width, units="px", tree_style=ts)
    cleaned_tree.write_text(t.write(format=1) + "\n", encoding="utf-8")

    with open(log, "a") as fh:
        fh.write("TREE_VISUALIZATION completed successfully\n")
        fh.write(f"Input tree: {tree_path}\n")
        fh.write(f"Output image: {out_img}\n")
        fh.write(f"Leaves rendered: {n_leaves}\n")
        fh.write(f"Removed reference tips: {', '.join(removed_refs) if removed_refs else 'none'}\n")
        fh.write(f"Resistance profile source: {resistance_profile_source}\n")
        fh.write(f"Canonical resistance profiles loaded: {len(resistance_profile_map)}\n")
        fh.write(f"Metadata rows loaded after normalization: {len(metadata)}\n")
        fh.write(f"Tree leaves matched to metadata after normalization: {metadata_matched}\n")

        if excluded_copy.exists():
            excluded_rows = 0
            try:
                with open(excluded_copy, newline="", encoding="utf-8", errors="replace") as excluded_handle:
                    excluded_reader = csv.DictReader(excluded_handle, delimiter="\t")
                    for excluded_row in excluded_reader:
                        if excluded_row.get("sample", "").strip():
                            excluded_rows += 1
                fh.write(f"IQ-TREE-excluded samples available for reporting: {excluded_rows}\n")
            except Exception as exc:
                fh.write(f"WARNING: Could not summarize IQ-TREE exclusion file: {repr(exc)}\n")

        if metadata_missing:
            fh.write("WARNING: Tree leaves without matching metadata after normalization: " + ", ".join(metadata_missing[:100]) + "\n")
        else:
            fh.write("All rendered tree leaves matched metadata after normalization.\n")

except Exception as e:
    with open(log, "a") as fh:
        fh.write("TREE_VISUALIZATION failed, but fallback outputs were created.\n")
        fh.write(f"ERROR: {repr(e)}\n")
        fh.write(traceback.format_exc() + "\n")

    if tree_input and Path(tree_input).exists() and Path(tree_input).stat().st_size > 0:
        try:
            raw = Path(tree_input).read_text(errors="replace")
            cleaned_tree.write_text(raw if raw.endswith("\n") else raw + "\n", encoding="utf-8")
        except Exception:
            cleaned_tree.write_text("();\n", encoding="utf-8")
    else:
        cleaned_tree.write_text("();\n", encoding="utf-8")

    write_placeholder_image(repr(e))
PY
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 50 HDD"
    timeout: "24 hours"
  }

  output {
    File tree_image = "tree_visualization/phylogenetic_tree.~{image_format}"
    File cleaned_tree = "tree_visualization/phylogenetic_tree.cleaned.nwk"
    File render_log = "tree_visualization/render.log"
    File excluded_from_iqtree = "tree_visualization/excluded_from_iqtree.tsv"
  }
}
task TB_SURVEILLANCE_SUMMARY_VISUALS {
  input {
    String docker_image = "python:3.11-slim"
    File tbprofiler_summary_tsv
    File? species_typing_tsv
    File? pairwise_snp_distance_matrix
    File? mean_depth_tsv
    File? sample_metadata_tsv
    Int cpu = 4
    Int memory_gb = 64
  }

  command <<<
    set -euo pipefail
    mkdir -p surveillance_summary

    # Defensive placeholders: Cromwell requires these declared outputs even if
    # plotting/summary rendering encounters malformed or unusually large fields.
    printf "lineage	count
Not available	0
" > surveillance_summary/lineage_distribution.tsv
    cat > surveillance_summary/lineage_distribution.svg <<'SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240"><rect width="100%" height="100%" fill="#fff"/><text x="450" y="110" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">Lineage Distribution</text><text x="450" y="145" text-anchor="middle" font-family="Arial" font-size="14" fill="#475569">No lineage summary was generated.</text></svg>
SVG
    cat > surveillance_summary/snp_distance_heatmap.svg <<'SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240"><rect width="100%" height="100%" fill="#fff"/><text x="450" y="110" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">SNP Distance Heatmap</text><text x="450" y="145" text-anchor="middle" font-family="Arial" font-size="14" fill="#475569">No SNP distance heatmap was generated.</text></svg>
SVG
    printf "sample	integrated_mtbc_status	mtbc_support_source	kraken_mtbc_supported	kraken_species	kraken_mtbc_reads	mtbc_percent	tbprofiler_lineage_status	tbprofiler_main_lineage	tbprofiler_sub_lineage	lineage_group	resistance_profile	drug_resistance_detected	resistant_drugs	mean_depth	selected_for_phylogeny	included_in_tree	selection_basis	tbprofiler_status
" > surveillance_summary/tb_surveillance_metadata.tsv
    printf "sample	mean_depth	kraken_mtbc_supported	mtbc_percent	selected_for_phylogeny	included_in_tree	reason
" > surveillance_summary/qc_filtering_rationale.tsv
    cat > surveillance_summary/surveillance_summary.html <<'HTML'
<!doctype html><html><head><meta charset="utf-8"><title>Surveillance Summary</title></head><body><h1>Surveillance Summary</h1><p>Summary was not generated; fallback output created to allow workflow completion.</p></body></html>
HTML

    tb_tsv="~{tbprofiler_summary_tsv}"
    species_tsv="~{if defined(species_typing_tsv) then species_typing_tsv else ""}"
    snp_matrix="~{if defined(pairwise_snp_distance_matrix) then pairwise_snp_distance_matrix else ""}"
    depth_tsv="~{if defined(mean_depth_tsv) then mean_depth_tsv else ""}"
    user_metadata_tsv="~{if defined(sample_metadata_tsv) then select_first([sample_metadata_tsv]) else ""}"

    set +e
    python3 - "$tb_tsv" "$species_tsv" "$snp_matrix" "$depth_tsv" "$user_metadata_tsv" <<'PY'
import csv
import html
import re
import sys
from pathlib import Path
from collections import Counter

# Some TB-Profiler/Kraken/Snippy-derived TSV fields can exceed Python's
# conservative csv module default of 131072 bytes, especially when mutation
# evidence or long JSON-derived strings are carried forward. Raise the limit so
# DictReader can parse large but valid TSV rows.
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2147483647)

tb_tsv = Path(sys.argv[1])
species_tsv = Path(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] else None
snp_matrix = Path(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else None
depth_tsv = Path(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else None
user_metadata_path = Path(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] else None

outdir = Path("surveillance_summary")
outdir.mkdir(exist_ok=True)

REFERENCE_NAMES = {
    "reference",
    "ref",
    "h37rv",
    "h37rvsiena",
    "nc_000962",
    "nc_000962.3",
    "mycobacterium_tuberculosis_h37rv"
}

def safe(x):
    return html.escape(str(x if x is not None else ""))

def is_reference_name(name):
    n = str(name or "").strip().lower()
    return (
        n in REFERENCE_NAMES
        or n.startswith("reference")
        or n.startswith("ref_")
        or "h37rv" in n
    )

def read_tsv(path):
    if path and path.exists() and path.stat().st_size > 0:
        with path.open() as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    return []

def clean_text(x):
    return str(x or "").strip()

def normalize_missing(x, replacement="Not reported"):
    v = clean_text(x)
    if v.lower() in ["", "none", "none reported", "not reported", "unknown", "na", "n/a"]:
        return replacement
    return v

def normalize_sample_id(name):
    s = str(name or "").strip()
    s = s.split("/")[-1]
    s = s.split("\\")[-1]
    s = re.sub(r"(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq|\.gz)$", "", s, flags=re.IGNORECASE)
    s = re.sub(
        r"(_R?1_paired|_R?2_paired|_R?1|_R?2|_1_paired|_2_paired|_1|_2|\.R?1|\.R?2|\.1|\.2)$",
        "",
        s
    )
    s = re.sub(r"^\d{3,}[_-](?=(?:[DES]RR|SAM[END]?|ERS|SRS|DRS|SRX|ERX|DRX|[A-Za-z]))", "", s)
    return s.strip()

def extract_mtbc_percent(evidence):
    text = str(evidence or "")

    m = re.search(
        r"MTBC\s+support:\s*[^;]*?;\s*(\d+(?:\.\d+)?)\s*%",
        text,
        flags=re.IGNORECASE
    )
    if m:
        return m.group(1)

    m = re.search(
        r"MTBC\s+support[^\d]*(\d+(?:\.\d+)?)\s*%",
        text,
        flags=re.IGNORECASE
    )
    if m:
        return m.group(1)

    percents = re.findall(r"(\d+(?:\.\d+)?)\s*%", text)
    return percents[-1] if percents else "Not available"

def normalize_lineage(main_lineage, sub_lineage):
    text = " ".join([str(main_lineage or ""), str(sub_lineage or "")]).strip()
    low = text.lower()

    unresolved = {
        "",
        "not reported",
        "none",
        "unknown",
        "na",
        "n/a",
        "not reported not reported",
        "not resolved by tb-profiler",
        "not resolved by tb-profiler not resolved by tb-profiler",
        "lineage not resolved"
    }

    if low in unresolved:
        return "Not resolved by TB-Profiler"

    m = re.search(r"lineage\s*([1-9])", text, flags=re.IGNORECASE)
    if m:
        return "L" + m.group(1)

    m = re.search(r"\bL([1-9])\b", text, flags=re.IGNORECASE)
    if m:
        return "L" + m.group(1)

    m = re.search(r"(^|[^0-9])([1-9])\.[0-9]", text)
    if m:
        return "L" + m.group(2)

    return text

def read_snp_matrix(path):
    if not path or not path.exists() or path.stat().st_size == 0:
        return [], []

    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, [])

        raw_samples = header[1:]
        keep_idx = [
            i for i, s in enumerate(raw_samples)
            if s and not is_reference_name(s)
        ]

        samples = [normalize_sample_id(raw_samples[i]) for i in keep_idx]
        matrix = []

        for row in reader:
            if not row:
                continue

            row_sample = row[0]

            if is_reference_name(row_sample):
                continue

            vals = []

            for i in keep_idx:
                try:
                    vals.append(float(row[i + 1]))
                except Exception:
                    vals.append(None)

            matrix.append(vals)

    return samples, matrix

def resistance_detected(profile, drugs):
    p = str(profile or "").strip().lower()
    d = str(drugs or "").strip().lower()

    no_res_profiles = {
        "",
        "none",
        "none reported",
        "not reported",
        "unknown",
        "not determined",
        "resistance not determined by tb-profiler",
        "no resistance detected by tb-profiler"
    }

    no_res_drugs = {
        "",
        "none",
        "none reported",
        "not reported",
        "unknown",
        "na",
        "n/a"
    }

    if p in no_res_profiles and d in no_res_drugs:
        return "NO"

    if "resistance not determined" in p:
        return "Not determined"

    if "no resistance detected" in p:
        return "NO"

    return "YES"

def infer_kraken_mtbc_support(row):
    value = (
        row.get("MTBC_Supported")
        or row.get("mtbc_supported")
        or row.get("MTBC supported")
        or ""
    )

    value_upper = str(value or "").strip().upper()

    if value_upper in {"YES", "Y", "TRUE", "1"}:
        return "YES"

    if value_upper in {"NO", "N", "FALSE", "0"}:
        return "NO"

    species = (
        row.get("Species_Identified")
        or row.get("Species identified")
        or row.get("species")
        or row.get("Species")
        or ""
    )

    species_text = str(species or "").strip().lower()

    mtbc_terms = [
        "mycobacterium tuberculosis complex",
        "tuberculosis complex",
        "mtbc",
        "mycobacterium tuberculosis",
        "m. tuberculosis",
        "mycobacterium africanum",
        "mycobacterium bovis",
        "mycobacterium caprae",
        "mycobacterium canettii",
        "mycobacterium microti",
        "mycobacterium pinnipedii",
        "mycobacterium orygis",
        "mycobacterium mungi",
        "mycobacterium suricattae",
        "dassie bacillus"
    ]

    if any(t in species_text for t in mtbc_terms):
        return "YES"

    if "mycobacterium" in species_text:
        return "NO"

    try:
        mtbc_pct = float(row.get("MTBC_Percent") or row.get("mtbc_percent") or 0)
    except Exception:
        mtbc_pct = 0.0
    try:
        mtbc_reads = int(float(row.get("MTBC_Reads") or row.get("mtbc_reads") or 0))
    except Exception:
        mtbc_reads = 0

    if mtbc_pct >= 10.0 and mtbc_reads >= 1000:
        return "YES"

    return "NO"

def split_integrated_mtbc_status_from_kraken(mtbc_supported):
    if str(mtbc_supported or "").strip().upper() == "YES":
        return "MTBC supported", "Kraken2/Bracken"

    return "MTBC not supported", "Kraken2/Bracken"

tb_rows = read_tsv(tb_tsv)
species_rows = read_tsv(species_tsv)
depth_rows = read_tsv(depth_tsv)

snp_samples, snp_values = read_snp_matrix(snp_matrix)
tree_sample_set = set(normalize_sample_id(s) for s in snp_samples)

depth_by_sample = {}

for r in depth_rows:
    sid = (
        r.get("sample")
        or r.get("Sample")
        or r.get("Sample_ID")
        or r.get("sample_id")
        or ""
    )

    depth = (
        r.get("mean_depth")
        or r.get("MeanDepth")
        or r.get("Mean_Depth")
        or r.get("mean_depth_value")
        or "Not available"
    )

    if sid:
        depth_by_sample[normalize_sample_id(sid)] = depth

species_by_sample = {}

for r in species_rows:
    sample = (
        r.get("Sample_ID")
        or r.get("sample")
        or r.get("sample_id")
        or r.get("Sample ID")
        or ""
    )

    species = (
        r.get("Species_Identified")
        or r.get("Species identified")
        or r.get("species")
        or r.get("Species")
        or ""
    )

    evidence = r.get("Evidence") or r.get("evidence") or ""

    mtbc_supported = infer_kraken_mtbc_support(r)

    mtbc_reads = (
        r.get("MTBC_Reads")
        or r.get("mtbc_reads")
        or r.get("MTBC reads")
        or "Not available"
    )

    mtbc_percent = (
        r.get("MTBC_Percent")
        or r.get("mtbc_percent")
        or r.get("MTBC percent")
        or extract_mtbc_percent(evidence)
    )

    selection_basis = (
        r.get("Selection_Basis")
        or r.get("selection_basis")
        or r.get("Selection basis")
        or ""
    )

    if sample:
        species_by_sample[normalize_sample_id(sample)] = {
            "species": species,
            "evidence": evidence,
            "mtbc_supported": mtbc_supported,
            "mtbc_reads": mtbc_reads,
            "mtbc_percent": mtbc_percent,
            "selection_basis": selection_basis
        }

lineage_counts = Counter()
metadata_rows = []
qc_rows = []

for r in tb_rows:
    sample = r.get("sample", "")
    sample_norm = normalize_sample_id(sample)

    sp = species_by_sample.get(sample_norm, {})
    kraken_species = sp.get("species", "")
    kraken_evidence = sp.get("evidence", "")
    kraken_mtbc_supported = sp.get("mtbc_supported", "NO")
    kraken_mtbc_reads = sp.get("mtbc_reads", "Not available")
    kraken_selection_basis = sp.get("selection_basis", "")

    main_lineage = normalize_missing(
        r.get("main_lineage", ""),
        "Not resolved by TB-Profiler"
    )

    sub_lineage = normalize_missing(
        r.get("sub_lineage", ""),
        "Not resolved by TB-Profiler"
    )

    lineage_group = normalize_lineage(main_lineage, sub_lineage)
    lineage_counts[lineage_group] += 1

    ###########################################################################
    # Kraken2/Bracken-based phylogeny selection
    #
    # selected_for_phylogeny comes from mtbc_selected in TB_PROFILER_AND_MTBC_FILTER,
    # which should now be based only on Kraken2/Bracken MTBC support.
    #
    # TB-Profiler lineage, species, and resistance are annotations only.
    ###########################################################################

    selected_for_phylogeny = (
        "YES"
        if str(r.get("mtbc_selected", "")).strip().upper() == "YES"
        else "NO"
    )

    if kraken_mtbc_supported == "YES" and selected_for_phylogeny != "YES":
        selected_for_phylogeny = "YES"

    if tree_sample_set:
        included_in_tree = "YES" if sample_norm in tree_sample_set else "NO"
    else:
        included_in_tree = "Not determined"

    mean_depth = depth_by_sample.get(sample_norm, "Not available")

    mtbc_percent = sp.get("mtbc_percent") or extract_mtbc_percent(kraken_evidence)

    resistance_profile = normalize_missing(
        r.get("dr_type", ""),
        "Resistance not determined by TB-Profiler"
    )

    resistant_drugs = normalize_missing(
        r.get("resistant_drugs", ""),
        "None reported"
    )

    drug_resistance = resistance_detected(resistance_profile, resistant_drugs)

    mtbc_status, support_source = split_integrated_mtbc_status_from_kraken(kraken_mtbc_supported)

    if lineage_group != "Not resolved by TB-Profiler":
        lineage_status = "Resolved by TB-Profiler"
    else:
        lineage_status = "Not resolved by TB-Profiler"

    selection_reason = (
        r.get("mtbc_selection_reason", "")
        or kraken_selection_basis
        or "Selected for phylogeny if Kraken2/Bracken species typing supports MTBC"
    )

    if selected_for_phylogeny == "NO" and not selection_reason:
        selection_reason = "Not selected for phylogeny because Kraken2/Bracken did not support MTBC"

    metadata_rows.append({
        "sample": sample,
        "integrated_mtbc_status": mtbc_status,
        "mtbc_support_source": support_source,
        "kraken_mtbc_supported": kraken_mtbc_supported,
        "kraken_species": kraken_species,
        "kraken_mtbc_reads": kraken_mtbc_reads,
        "mtbc_percent": mtbc_percent,
        "tbprofiler_lineage_status": lineage_status,
        "tbprofiler_main_lineage": main_lineage,
        "tbprofiler_sub_lineage": sub_lineage,
        "lineage_group": lineage_group,
        "resistance_profile": resistance_profile,
        "drug_resistance_detected": drug_resistance,
        "resistant_drugs": resistant_drugs,
        "mean_depth": mean_depth,
        "selected_for_phylogeny": selected_for_phylogeny,
        "included_in_tree": included_in_tree,
        "selection_basis": selection_reason,
        "tbprofiler_status": r.get("status", "")
    })

    qc_rows.append({
        "sample": sample,
        "mean_depth": mean_depth,
        "kraken_mtbc_supported": kraken_mtbc_supported,
        "mtbc_percent": mtbc_percent,
        "selected_for_phylogeny": selected_for_phylogeny,
        "included_in_tree": included_in_tree,
        "reason": selection_reason
    })

with (outdir / "lineage_distribution.tsv").open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["lineage", "count"])

    items = sorted(lineage_counts.items()) if lineage_counts else [("Not resolved by TB-Profiler", 0)]

    for lineage, count in items:
        writer.writerow([lineage, count])

# Optional user-supplied sample metadata TSV.
# Expected format: tab-delimited file with one row per sample and a sample identifier column.
# Supported identifier column names: sample, sample_id, isolate, isolate_id, accession, run_accession, sra_accession.
# Extra metadata columns are carried forward into tb_surveillance_metadata.tsv with a user_ prefix.
user_metadata_extra_columns = []
user_metadata_by_sample = {}
if user_metadata_path and Path(user_metadata_path).exists():
    with open(user_metadata_path, newline="") as meta_handle:
        reader = csv.DictReader(meta_handle, delimiter="\t")
        if reader.fieldnames:
            sample_keys = [
                "sample",
                "sample_id",
                "isolate",
                "isolate_id",
                "accession",
                "run_accession",
                "sra_accession"
            ]
            sample_key = next((k for k in sample_keys if k in reader.fieldnames), None)
            if sample_key:
                user_metadata_extra_columns = [c for c in reader.fieldnames if c != sample_key]
                for meta_row in reader:
                    sample_id = (meta_row.get(sample_key) or "").strip()
                    if sample_id:
                        user_metadata_by_sample[sample_id] = {
                            f"user_{c}": (meta_row.get(c) or "") for c in user_metadata_extra_columns
                        }

if user_metadata_by_sample:
    for row in metadata_rows:
        row.update(user_metadata_by_sample.get(row.get("sample", ""), {}))

with (outdir / "tb_surveillance_metadata.tsv").open("w", newline="") as out:
    fieldnames = [
        "sample",
        "integrated_mtbc_status",
        "mtbc_support_source",
        "kraken_mtbc_supported",
        "kraken_species",
        "kraken_mtbc_reads",
        "mtbc_percent",
        "tbprofiler_lineage_status",
        "tbprofiler_main_lineage",
        "tbprofiler_sub_lineage",
        "lineage_group",
        "resistance_profile",
        "drug_resistance_detected",
        "resistant_drugs",
        "mean_depth",
        "selected_for_phylogeny",
        "included_in_tree",
        "selection_basis",
        "tbprofiler_status"
    ]
    fieldnames = fieldnames + [f"user_{c}" for c in user_metadata_extra_columns]

    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(metadata_rows)

with (outdir / "qc_filtering_rationale.tsv").open("w", newline="") as out:
    fieldnames = [
        "sample",
        "mean_depth",
        "kraken_mtbc_supported",
        "mtbc_percent",
        "selected_for_phylogeny",
        "included_in_tree",
        "reason"
    ]

    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(qc_rows)

def write_lineage_svg(counts, outfile):
    width, height = 1000, 520
    ml, mb, mt = 120, 80, 70
    pw, ph = width - ml - 60, height - mt - mb

    items = sorted(counts.items()) if counts else [("Not resolved by TB-Profiler", 0)]
    maxc = max([c for _, c in items] + [1])

    gap = 24
    bw = max(45, int((pw - gap * (len(items) + 1)) / max(len(items), 1)))

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width/2}" y="35" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700" fill="#0f172a">Lineage Distribution</text>',
        f'<text x="{width/2}" y="58" text-anchor="middle" font-family="Arial" font-size="12" fill="#475569">TB-Profiler lineage calls where available; unresolved Kraken-supported MTBC samples are shown separately.</text>',
        f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#334155" stroke-width="2"/>',
        f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#334155" stroke-width="2"/>'
    ]

    for i, (lab, c) in enumerate(items):
        x = ml + gap + i * (bw + gap)
        bh = int((c / maxc) * (ph - 20)) if maxc else 0
        y = mt + ph - bh

        svg += [
            f'<rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="8" fill="#2563eb"/>',
            f'<text x="{x+bw/2}" y="{y-8}" text-anchor="middle" font-family="Arial" font-size="14" font-weight="700" fill="#0f172a">{c}</text>',
            f'<text x="{x+bw/2}" y="{mt+ph+30}" text-anchor="middle" font-family="Arial" font-size="13" fill="#0f172a">{safe(lab)}</text>'
        ]

    svg.append(
        f'<text x="{width/2}" y="{height-24}" text-anchor="middle" font-family="Arial" font-size="12" fill="#475569">Lineage groups are summarized from TB-Profiler main-lineage and sub-lineage fields. Lineage does not determine phylogeny selection.</text></svg>'
    )

    outfile.write_text("\n".join(svg), encoding="utf-8")

def heat_color(v, maxv):
    if v is None:
        return "#e5e7eb"

    if maxv <= 0:
        return "#dcfce7"

    r = max(0, min(v / maxv, 1))

    if r <= 0.25:
        return "#dcfce7"
    if r <= 0.5:
        return "#fef9c3"
    if r <= 0.75:
        return "#fed7aa"

    return "#fecaca"

def write_heatmap(path, outfile):
    samples, matrix = read_snp_matrix(path)

    if not samples or not matrix:
        outfile.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240">'
            '<rect width="100%" height="100%" fill="#fff"/>'
            '<text x="450" y="100" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">SNP Distance Heatmap</text>'
            '<text x="450" y="140" text-anchor="middle" font-family="Arial" font-size="14" fill="#475569">No SNP distance matrix available after reference filtering.</text>'
            '</svg>',
            encoding="utf-8"
        )
        return

    n = len(samples)
    cell = 46 if n <= 12 else max(18, int(760 / n))

    left = 220
    top = 210

    width = left + n * cell + 110
    height = top + n * cell + 140

    vals = [
        v
        for row in matrix
        for v in row
        if v is not None
    ]

    maxv = max(vals + [1])

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width/2}" y="45" text-anchor="middle" font-family="Arial" font-size="28" font-weight="700" fill="#0f172a">Pairwise SNP Distance Heatmap</text>'
    ]

    for i, s in enumerate(samples):
        x = left + i * cell + cell / 2
        label_y = top - 32

        svg.append(
            f'<text x="{x}" y="{label_y}" transform="rotate(-55 {x},{label_y})" text-anchor="start" font-family="Arial" font-size="11" fill="#0f172a">{safe(s)}</text>'
        )

    for i, s in enumerate(samples):
        y = top + i * cell + cell / 2 + 4

        svg.append(
            f'<text x="{left - 10}" y="{y}" text-anchor="end" font-family="Arial" font-size="11" fill="#0f172a">{safe(s)}</text>'
        )

    for i, row in enumerate(matrix):
        for j, v in enumerate(row[:n]):
            x = left + j * cell
            y = top + i * cell

            if v is None:
                lab = ""
            else:
                lab = str(int(v)) if float(v).is_integer() else f"{v:.1f}"

            svg.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{heat_color(v, maxv)}" stroke="#ffffff" stroke-width="1"/>'
            )

            if n <= 18:
                svg.append(
                    f'<text x="{x + cell/2}" y="{y + cell/2 + 4}" text-anchor="middle" font-family="Arial" font-size="10" fill="#111827">{safe(lab)}</text>'
                )

    footer_y = height - 58

    svg.append(
        f'<text x="{width/2}" y="{footer_y}" text-anchor="middle" font-family="Arial" font-size="12" fill="#475569">Reference sequence excluded. Lower values indicate closer genomic relatedness.</text>'
    )

    svg.append(
        f'<text x="{width/2}" y="{footer_y + 22}" text-anchor="middle" font-family="Arial" font-size="12" fill="#475569">Color scale: green = closest, yellow/orange = intermediate, red = more distant.</text>'
    )

    svg.append("</svg>")

    outfile.write_text("\n".join(svg), encoding="utf-8")

write_lineage_svg(lineage_counts, outdir / "lineage_distribution.svg")
write_heatmap(snp_matrix, outdir / "snp_distance_heatmap.svg")

summary_rows = []
summary_rows.append('<!doctype html>')
summary_rows.append('<html>')
summary_rows.append('<head>')
summary_rows.append('<meta charset="utf-8">')
summary_rows.append('<title>Surveillance Summary</title>')
summary_rows.append('<style>')
summary_rows.append('body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#f8fafc;color:#0f172a}')
summary_rows.append('.card{background:#ffffff;border:1px solid #e2e8f0;border-radius:14px;padding:18px;margin-bottom:18px;box-shadow:0 1px 4px rgba(15,23,42,.08)}')
summary_rows.append('table{border-collapse:collapse;width:100%;font-size:13px;background:#ffffff}')
summary_rows.append('th,td{border:1px solid #e2e8f0;padding:8px;vertical-align:top}')
summary_rows.append('th{background:#1d4ed8;color:#ffffff}')
summary_rows.append('.muted{color:#475569}')
summary_rows.append('.yes{display:inline-block;border-radius:999px;padding:4px 8px;background:#dcfce7;color:#166534;font-weight:700}')
summary_rows.append('.no{display:inline-block;border-radius:999px;padding:4px 8px;background:#fee2e2;color:#991b1b;font-weight:700}')
summary_rows.append('.unknown{display:inline-block;border-radius:999px;padding:4px 8px;background:#e5e7eb;color:#374151;font-weight:700}')
summary_rows.append('</style>')
summary_rows.append('</head>')
summary_rows.append('<body>')
summary_rows.append('<div class="card">')
summary_rows.append('<h1>Surveillance Summary</h1>')
summary_rows.append('<p class="muted">This summary integrates Kraken2/Bracken species typing, TB-Profiler lineage and resistance annotations, mean-depth metadata, and SNP-distance outputs.</p>')
summary_rows.append('<p class="muted"><strong>Phylogeny selection rule:</strong> samples are selected for Snippy/core-SNP/IQ-TREE only when Kraken2/Bracken supports MTBC. TB-Profiler lineage, species, and resistance are reported as annotations and do not determine phylogeny selection.</p>')
summary_rows.append('<p class="muted"><strong>Tree inclusion rule:</strong> after Kraken2/Bracken-based selection, IQ-TREE-related outputs reflect only samples that passed downstream alignment and SNP-distance/tree generation steps.</p>')
summary_rows.append('</div>')

summary_rows.append('<div class="card">')
summary_rows.append('<h2>Kraken2/Bracken-based phylogeny selection and tree inclusion</h2>')
summary_rows.append('<table>')
summary_rows.append('<thead><tr>')
summary_rows.append('<th>Sample</th>')
summary_rows.append('<th>Kraken MTBC supported</th>')
summary_rows.append('<th>Selected for phylogeny</th>')
summary_rows.append('<th>Included in tree/SNP matrix</th>')
summary_rows.append('<th>MTBC percent</th>')
summary_rows.append('<th>Mean depth</th>')
summary_rows.append('<th>Selection basis</th>')
summary_rows.append('</tr></thead>')
summary_rows.append('<tbody>')

for row in metadata_rows:
    selected = row.get("selected_for_phylogeny", "")
    included = row.get("included_in_tree", "")
    kraken_supported = row.get("kraken_mtbc_supported", "")

    selected_class = "yes" if selected == "YES" else "no"
    kraken_class = "yes" if kraken_supported == "YES" else "no"

    if included == "YES":
        included_class = "yes"
    elif included == "NO":
        included_class = "no"
    else:
        included_class = "unknown"

    summary_rows.append(
        "<tr>"
        f"<td>{safe(row.get('sample', ''))}</td>"
        f"<td><span class='{kraken_class}'>{safe(kraken_supported)}</span></td>"
        f"<td><span class='{selected_class}'>{safe(selected)}</span></td>"
        f"<td><span class='{included_class}'>{safe(included)}</span></td>"
        f"<td>{safe(row.get('mtbc_percent', ''))}</td>"
        f"<td>{safe(row.get('mean_depth', ''))}</td>"
        f"<td>{safe(row.get('selection_basis', ''))}</td>"
        "</tr>"
    )

summary_rows.append('</tbody></table>')
summary_rows.append('</div>')
summary_rows.append('</body></html>')

(outdir / "surveillance_summary.html").write_text(
    "\n".join(summary_rows),
    encoding="utf-8"
)
PY
    py_rc=$?
    set -e

    if [ "$py_rc" -ne 0 ]; then
      echo "WARNING: TB_SURVEILLANCE_SUMMARY_VISUALS Python rendering failed with exit code ${py_rc}." >&2
      echo "WARNING: Keeping fallback surveillance summary outputs so downstream MERGE_TB_REPORTS can proceed." >&2
    fi

    # Final output guard: never let this optional visualization task fail during
    # output collection because of a missing declared file.
    [ -s surveillance_summary/lineage_distribution.tsv ] || printf "lineage	count
Not available	0
" > surveillance_summary/lineage_distribution.tsv
    [ -s surveillance_summary/tb_surveillance_metadata.tsv ] || printf "sample	integrated_mtbc_status	mtbc_support_source	kraken_mtbc_supported	kraken_species	kraken_mtbc_reads	mtbc_percent	tbprofiler_lineage_status	tbprofiler_main_lineage	tbprofiler_sub_lineage	lineage_group	resistance_profile	drug_resistance_detected	resistant_drugs	mean_depth	selected_for_phylogeny	included_in_tree	selection_basis	tbprofiler_status
" > surveillance_summary/tb_surveillance_metadata.tsv
    [ -s surveillance_summary/qc_filtering_rationale.tsv ] || printf "sample	mean_depth	kraken_mtbc_supported	mtbc_percent	selected_for_phylogeny	included_in_tree	reason
" > surveillance_summary/qc_filtering_rationale.tsv
    [ -s surveillance_summary/lineage_distribution.svg ] || cat > surveillance_summary/lineage_distribution.svg <<'SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240"><rect width="100%" height="100%" fill="#fff"/><text x="450" y="110" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">Lineage Distribution</text><text x="450" y="145" text-anchor="middle" font-family="Arial" font-size="14" fill="#475569">No lineage summary was generated.</text></svg>
SVG
    [ -s surveillance_summary/snp_distance_heatmap.svg ] || cat > surveillance_summary/snp_distance_heatmap.svg <<'SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240"><rect width="100%" height="100%" fill="#fff"/><text x="450" y="110" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">SNP Distance Heatmap</text><text x="450" y="145" text-anchor="middle" font-family="Arial" font-size="14" fill="#475569">No SNP distance heatmap was generated.</text></svg>
SVG
    [ -s surveillance_summary/surveillance_summary.html ] || cat > surveillance_summary/surveillance_summary.html <<'HTML'
<!doctype html><html><head><meta charset="utf-8"><title>Surveillance Summary</title></head><body><h1>Surveillance Summary</h1><p>Summary was not generated; fallback output created to allow workflow completion.</p></body></html>
HTML
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 100 HDD"
    timeout: "24 hours"
  }

  output {
    File lineage_distribution_tsv = "surveillance_summary/lineage_distribution.tsv"
    File lineage_distribution_svg = "surveillance_summary/lineage_distribution.svg"
    File snp_distance_heatmap_svg = "surveillance_summary/snp_distance_heatmap.svg"
    File surveillance_metadata_tsv = "surveillance_summary/tb_surveillance_metadata.tsv"
    File qc_filtering_rationale_tsv = "surveillance_summary/qc_filtering_rationale.tsv"
    File surveillance_summary_html = "surveillance_summary/surveillance_summary.html"
  }
}
task MERGE_TB_REPORTS {
  input {
    String docker_image = "python:3.11-slim"
    File? tbprofiler_html
    File? tbprofiler_summary_tsv
    File? resistance_profile_summary_tsv
    File? tbprofiler_mutation_evidence_tsv
    File? tbprofiler_mutation_evidence_html
    File? mtbc_samples_txt
    File? species_typing_html
    File? species_typing_tsv
    File? ntm_species_summary_html
    File? ntm_species_summary_tsv
    File? mycobacteria_routing_summary_tsv
    File? mycobacteria_routing_status
    File? qc_summary_html
    File? trimming_report_html
    File? variant_summary_html
    File? iqtree_report
    File? iqtree_excluded_samples_tsv
    File? iqtree_included_samples_tsv
    File? iqtree_filtering_summary_txt
    File? iqtree_status
    File? tree_image
    File? phylogenetic_tree_newick
    File? pairwise_tree_newick
    File? nonsynonymous_mutations_tsv
    File? nonsynonymous_mutations_html
    File? pairwise_snp_distance_matrix
    File? pairwise_snp_distance_pairs
    File? snp_cluster_summary
    File? snp_distance_cluster_html
    File? lineage_distribution_tsv
    File? lineage_distribution_svg
    File? snp_distance_heatmap_svg
    File? surveillance_metadata_tsv
    File? qc_filtering_rationale_tsv
    File? surveillance_summary_html
  }

  command <<<
    set -uo pipefail
    mkdir -p final_report

    tb_tsv="~{if defined(tbprofiler_summary_tsv) then tbprofiler_summary_tsv else ""}"
    resistance_tsv="~{if defined(resistance_profile_summary_tsv) then resistance_profile_summary_tsv else ""}"
    species_html="~{if defined(species_typing_html) then species_typing_html else ""}"
    species_tsv="~{if defined(species_typing_tsv) then species_typing_tsv else ""}"
    ntm_html="~{if defined(ntm_species_summary_html) then ntm_species_summary_html else ""}"
    ntm_tsv="~{if defined(ntm_species_summary_tsv) then ntm_species_summary_tsv else ""}"
    routing_tsv="~{if defined(mycobacteria_routing_summary_tsv) then mycobacteria_routing_summary_tsv else ""}"
    routing_status_txt="~{if defined(mycobacteria_routing_status) then mycobacteria_routing_status else ""}"
    nonsyn_tsv="~{if defined(nonsynonymous_mutations_tsv) then nonsynonymous_mutations_tsv else ""}"
    mutation_tsv="~{if defined(tbprofiler_mutation_evidence_tsv) then tbprofiler_mutation_evidence_tsv else ""}"
    snp_pairs_tsv="~{if defined(pairwise_snp_distance_pairs) then pairwise_snp_distance_pairs else ""}"
    snp_cluster_tsv="~{if defined(snp_cluster_summary) then snp_cluster_summary else ""}"
    lineage_tsv="~{if defined(lineage_distribution_tsv) then lineage_distribution_tsv else ""}"
    lineage_svg="~{if defined(lineage_distribution_svg) then lineage_distribution_svg else ""}"
    snp_heatmap_svg="~{if defined(snp_distance_heatmap_svg) then snp_distance_heatmap_svg else ""}"
    surveillance_metadata_tsv="~{if defined(surveillance_metadata_tsv) then surveillance_metadata_tsv else ""}"
    qc_rationale_tsv="~{if defined(qc_filtering_rationale_tsv) then qc_filtering_rationale_tsv else ""}"
    tree_png="~{if defined(tree_image) then tree_image else ""}"
    qc_html="~{if defined(qc_summary_html) then qc_summary_html else ""}"
    trim_html="~{if defined(trimming_report_html) then trimming_report_html else ""}"
    variant_html="~{if defined(variant_summary_html) then variant_summary_html else ""}"
    iqtree_txt="~{if defined(iqtree_report) then iqtree_report else ""}"
    iqtree_excluded_tsv="~{if defined(iqtree_excluded_samples_tsv) then iqtree_excluded_samples_tsv else ""}"
    iqtree_included_tsv="~{if defined(iqtree_included_samples_tsv) then iqtree_included_samples_tsv else ""}"
    iqtree_filtering_summary_txt="~{if defined(iqtree_filtering_summary_txt) then iqtree_filtering_summary_txt else ""}"
    iqtree_status_txt="~{if defined(iqtree_status) then iqtree_status else ""}"

    if [ -n "$tree_png" ] && [ -f "$tree_png" ]; then
      cp "$tree_png" final_report/mtbc_tree.png || true
    fi

    if [ -n "$iqtree_excluded_tsv" ] && [ -f "$iqtree_excluded_tsv" ]; then
      cp "$iqtree_excluded_tsv" final_report/excluded_from_iqtree.tsv || true
    fi

    if [ -n "$iqtree_included_tsv" ] && [ -f "$iqtree_included_tsv" ]; then
      cp "$iqtree_included_tsv" final_report/included_in_iqtree.tsv || true
    fi

    if [ -n "$iqtree_filtering_summary_txt" ] && [ -f "$iqtree_filtering_summary_txt" ]; then
      cp "$iqtree_filtering_summary_txt" final_report/alignment_filtering_summary.txt || true
    fi

    if [ -n "$iqtree_status_txt" ] && [ -f "$iqtree_status_txt" ]; then
      cp "$iqtree_status_txt" final_report/iqtree_status.txt || true
    fi

    if [ -n "$lineage_svg" ] && [ -f "$lineage_svg" ]; then
      cp "$lineage_svg" final_report/lineage_distribution.svg || true
    fi

    if [ -n "$snp_heatmap_svg" ] && [ -f "$snp_heatmap_svg" ]; then
      cp "$snp_heatmap_svg" final_report/snp_distance_heatmap.svg || true
    fi

    if [ -n "$surveillance_metadata_tsv" ] && [ -f "$surveillance_metadata_tsv" ]; then
      cp "$surveillance_metadata_tsv" final_report/tb_surveillance_metadata.tsv || true
    fi

    if [ -n "$qc_rationale_tsv" ] && [ -f "$qc_rationale_tsv" ]; then
      cp "$qc_rationale_tsv" final_report/qc_filtering_rationale.tsv || true
    fi

    if [ -n "$resistance_tsv" ] && [ -f "$resistance_tsv" ]; then
      cp "$resistance_tsv" final_report/resistance_profile_summary.tsv || true
    fi

    if [ -z "$tb_tsv" ] || [ ! -f "$tb_tsv" ]; then
      echo -e "sample\tspecies\tmain_lineage\tsub_lineage\tdr_type\tresistant_drugs\tresistance_mutations\tkey_mutations\tjson_file\tmtbc_selected\tmtbc_selection_reason\tstatus" > final_report/empty.tsv
      tb_tsv="final_report/empty.tsv"
    fi

    if [ -z "$resistance_tsv" ] || [ ! -f "$resistance_tsv" ]; then
      echo -e "sample_id\tresistance_profile\tresistant_drugs\tresistance_mutations\tkey_mutations\tstatus" > final_report/empty_resistance_profile_summary.tsv
      resistance_tsv="final_report/empty_resistance_profile_summary.tsv"
    fi

    if [ -z "$species_tsv" ] || [ ! -f "$species_tsv" ]; then
      echo -e "Sample_ID\tSpecies_Identified\tEvidence" > final_report/empty_species_typing.tsv
      species_tsv="final_report/empty_species_typing.tsv"
    fi

    if [ -z "$nonsyn_tsv" ] || [ ! -f "$nonsyn_tsv" ]; then
      echo -e "sample\tgene\teffect\taa_change\tnt_change\tproduct" > final_report/empty_nonsyn.tsv
      nonsyn_tsv="final_report/empty_nonsyn.tsv"
    fi

    if [ -z "$ntm_tsv" ] || [ ! -f "$ntm_tsv" ]; then
      echo -e "Sample_ID\tMost_Probable_Mycobacterium_Species\tMycobacteria_Category\tMTBC_Supported\tMTBC_Reads\tMTBC_Percent\tEvidence\tWorkflow_Decision\tSelection_Basis" > final_report/empty_ntm_species_summary.tsv
      ntm_tsv="final_report/empty_ntm_species_summary.tsv"
    fi

    if [ -z "$routing_tsv" ] || [ ! -f "$routing_tsv" ]; then
      echo -e "Sample_ID\tMost_Probable_Mycobacterium_Species\tMycobacteria_Category\tMTBC_Supported\tMTBC_Reads\tMTBC_Percent\tEvidence\tWorkflow_Decision\tSelection_Basis" > final_report/empty_mycobacteria_routing_summary.tsv
      routing_tsv="final_report/empty_mycobacteria_routing_summary.tsv"
    fi

    if [ -z "$routing_status_txt" ] || [ ! -f "$routing_status_txt" ]; then
      echo "Mycobacteria routing summary was not provided." > final_report/empty_mycobacteria_routing_status.txt
      routing_status_txt="final_report/empty_mycobacteria_routing_status.txt"
    fi

    if [ -z "$mutation_tsv" ] || [ ! -f "$mutation_tsv" ]; then
      echo -e "sample\tdrug\tgene\tmutation\tchange\tconfidence\tevidence\tsource_json" > final_report/empty_mutation_evidence.tsv
      mutation_tsv="final_report/empty_mutation_evidence.tsv"
    fi

    if [ -z "$snp_pairs_tsv" ] || [ ! -f "$snp_pairs_tsv" ]; then
      echo -e "sample1\tsample2\tsnp_distance\tcomparable_sites\tinterpretation\tcluster_class" > final_report/empty_snp_pairs.tsv
      snp_pairs_tsv="final_report/empty_snp_pairs.tsv"
    fi

    if [ -z "$snp_cluster_tsv" ] || [ ! -f "$snp_cluster_tsv" ]; then
      echo -e "cluster_id\tsample1\tsample2\tsnp_distance\tinterpretation" > final_report/empty_snp_clusters.tsv
      snp_cluster_tsv="final_report/empty_snp_clusters.tsv"
    fi

    if [ -z "$iqtree_excluded_tsv" ] || [ ! -f "$iqtree_excluded_tsv" ]; then
      echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note" > final_report/empty_excluded_from_iqtree.tsv
      iqtree_excluded_tsv="final_report/empty_excluded_from_iqtree.tsv"
    fi

    if [ -z "$iqtree_included_tsv" ] || [ ! -f "$iqtree_included_tsv" ]; then
      echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction" > final_report/empty_included_in_iqtree.tsv
      iqtree_included_tsv="final_report/empty_included_in_iqtree.tsv"
    fi

    if [ -z "$iqtree_filtering_summary_txt" ] || [ ! -f "$iqtree_filtering_summary_txt" ]; then
      echo "IQ-TREE alignment filtering summary was not provided." > final_report/empty_alignment_filtering_summary.txt
      iqtree_filtering_summary_txt="final_report/empty_alignment_filtering_summary.txt"
    fi

    if [ -z "$iqtree_status_txt" ] || [ ! -f "$iqtree_status_txt" ]; then
      echo "unknown_status" > final_report/empty_iqtree_status.txt
      iqtree_status_txt="final_report/empty_iqtree_status.txt"
    fi

    python3 - "$tb_tsv" "$resistance_tsv" "$species_tsv" "$ntm_tsv" "$routing_tsv" "$routing_status_txt" "$nonsyn_tsv" "$mutation_tsv" "$snp_pairs_tsv" "$snp_cluster_tsv" "$lineage_tsv" "$surveillance_metadata_tsv" "$qc_rationale_tsv" "$qc_html" "$trim_html" "$variant_html" "$iqtree_txt" "$species_html" "$ntm_html" "$iqtree_excluded_tsv" "$iqtree_included_tsv" "$iqtree_filtering_summary_txt" "$iqtree_status_txt" <<'PY'
import base64
import csv
import html
import os
import re
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone

summary_tsv, resistance_tsv, species_tsv, ntm_tsv, routing_tsv, routing_status_txt, nonsyn_tsv, mutation_tsv, snp_pairs_tsv, snp_cluster_tsv, lineage_tsv, surveillance_metadata_tsv, qc_rationale_tsv, qc_html, trim_html, variant_html, iqtree_txt, species_html, ntm_html, iqtree_excluded_tsv, iqtree_included_tsv, iqtree_filtering_summary_txt, iqtree_status_txt = sys.argv[1:24]

outdir = Path("final_report")
outdir.mkdir(exist_ok=True)

run_started_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
run_stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_UTC")

(outdir / "run_metadata.txt").write_text(
    f"Workflow report generation timestamp: {run_started_utc}\n"
    f"Run stamp: {run_stamp}\n"
)

RESISTANCE_COLORS = {
    # Deliberately unique colors for every final resistance category.
    # These colors are used consistently in the tree, legends, badges, and
    # surveillance tables so categories are not visually collapsed.
    "Sensitive": "#2A9D8F",                               # teal
    "Hr-TB": "#1D9BF0",                                  # blue
    "RR-TB": "#F59E0B",                                  # amber
    "MDR-TB": "#E11D48",                                 # rose
    "MDR/RR-TB": "#F97316",                              # orange
    "Pre-XDR-TB": "#EF4444",                             # red
    "XDR-TB": "#6D28D9",                                 # purple
    "Monoresistance": "#FACC15",                         # yellow
    "Polyresistance": "#06B6D4",                         # cyan
    "Other drug resistance": "#2563EB",                   # royal blue
    "Resistance not determined by TB-Profiler": "#94A3B8", # slate gray
    "Unknown": "#6B7280"                                  # dark gray
}

RESISTANCE_ORDER = [
    "Sensitive",
    "Hr-TB",
    "RR-TB",
    "MDR-TB",
    "MDR/RR-TB",
    "Pre-XDR-TB",
    "XDR-TB",
    "Monoresistance",
    "Polyresistance",
    "Other drug resistance",
    "Resistance not determined by TB-Profiler",
    "Unknown"
]

REFERENCE_NAMES = {
    "reference",
    "ref",
    "h37rv",
    "h37rvsiena",
    "nc_000962",
    "nc_000962.3",
    "mycobacterium_tuberculosis_h37rv"
}

def safe(v):
    return html.escape(str(v if v is not None else ""))

def clean(v):
    return str(v if v is not None else "").strip()

def lower_clean(v):
    return clean(v).lower()

def normalize_sample_id(name):
    s = clean(name)
    s = s.split("/")[-1]
    s = s.split("\\")[-1]
    s = re.sub(r"(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq|\.gz)$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"(_R?1_paired|_R?2_paired|_R?1|_R?2|_1_paired|_2_paired|_1|_2|\.R?1|\.R?2|\.1|\.2)$", "", s)
    # Remove synthetic leading numeric prefixes introduced by Terra/WDL set-mode,
    # e.g. 0000_ERR17070185 -> ERR17070185. Analysis files can keep the raw ID;
    # this function is used for report display, matching, and downloadable tables.
    s = re.sub(r"^\d{3,}[_-](?=(?:[DES]RR|SAM[END]?|ERS|SRS|DRS|SRX|ERX|DRX|[A-Za-z]))", "", s)
    return s.strip()

def is_reference_name(name):
    n = lower_clean(name)
    return n in REFERENCE_NAMES or n.startswith("reference") or n.startswith("ref_") or "h37rv" in n

def read_optional_file(path):
    p = Path(path)
    if path and p.exists() and p.stat().st_size > 0:
        try:
            return p.read_text(errors="replace")
        except Exception:
            return ""
    return ""

def read_tsv_rows(path):
    try:
        p = Path(path)
        if path and p.exists() and p.is_file() and p.stat().st_size > 0:
            with p.open() as fh:
                return list(csv.DictReader(fh, delimiter="\t"))
    except Exception:
        return []
    return []

def normalize_missing(v, replacement="Not reported"):
    x = clean(v)
    if lower_clean(x) in ["", "none", "none reported", "not reported", "unknown", "na", "n/a"]:
        return replacement
    return x

rows = read_tsv_rows(summary_tsv)
resistance_rows = read_tsv_rows(resistance_tsv)
species_rows = read_tsv_rows(species_tsv)
ntm_rows = read_tsv_rows(ntm_tsv)
routing_rows = read_tsv_rows(routing_tsv)
routing_status_text = read_optional_file(routing_status_txt).strip() or "Mycobacteria routing summary was not provided."

# -------------------------------------------------------------------------
# Single source of truth for MTBC/NTM routing in the final report.
# -------------------------------------------------------------------------
# The final HTML must not infer MTBC retention from TB-Profiler, Snippy, or
# IQ-TREE outputs. Those are downstream results. Routing comes first and is
# determined by MYCOBACTERIA_SAMPLE_ROUTER using Kraken2/Bracken support.
# These maps are used across dashboard counters, NTM tables, TB-Profiler rows,
# mutation evidence, lineage summaries, QC metadata, and tree notes.
# -------------------------------------------------------------------------

routing_by_sample = {}
mtbc_routed_samples = set()
non_mtbc_routed_samples = set()
ntm_display_rows_from_routing = []

for rr in routing_rows:
    sid = normalize_sample_id(rr.get("Sample_ID") or rr.get("sample") or rr.get("sample_id") or "")
    if not sid:
        continue
    category = clean(rr.get("Mycobacteria_Category") or rr.get("Category") or "")
    routing_by_sample[sid] = rr
    if category == "MTBC":
        mtbc_routed_samples.add(sid)
    else:
        non_mtbc_routed_samples.add(sid)
        ntm_display_rows_from_routing.append(rr)

def is_sample_routed_mtbc(sample):
    sid = normalize_sample_id(sample)
    if sid in mtbc_routed_samples:
        return True
    if sid in non_mtbc_routed_samples:
        return False
    # Fallback for older runs without router TSV: trust explicit TB-Profiler
    # mtbc_selected only when router information is unavailable.
    return False

def is_sample_routed_non_mtbc(sample):
    sid = normalize_sample_id(sample)
    return sid in non_mtbc_routed_samples

nonsyn_rows = read_tsv_rows(nonsyn_tsv)
mutation_rows = read_tsv_rows(mutation_tsv)
snp_pair_rows = read_tsv_rows(snp_pairs_tsv)
snp_cluster_rows = read_tsv_rows(snp_cluster_tsv)
lineage_rows = read_tsv_rows(lineage_tsv)
surveillance_metadata_rows = read_tsv_rows(surveillance_metadata_tsv)
qc_rationale_rows = read_tsv_rows(qc_rationale_tsv)
iqtree_excluded_rows = read_tsv_rows(iqtree_excluded_tsv)
iqtree_included_rows = read_tsv_rows(iqtree_included_tsv)
iqtree_filtering_summary_text = read_optional_file(iqtree_filtering_summary_txt)
iqtree_status_text = read_optional_file(iqtree_status_txt).strip() or "unknown_status"

# -------------------------------------------------------------------------
# IQ-TREE inclusion/exclusion reporting maps
# -------------------------------------------------------------------------
# Important distinction:
#   Selected for MTBC workflow = sample retained by MTBC support logic
#   Included in IQ-TREE        = sample retained after core-SNP alignment
#                                quality filtering for IQ-TREE inference
# -------------------------------------------------------------------------

iqtree_included_map = {}
iqtree_excluded_map = {}

for r in iqtree_included_rows:
    sample = normalize_sample_id(r.get("sample") or r.get("Sample") or "")
    if sample:
        iqtree_included_map[sample] = r

for r in iqtree_excluded_rows:
    sample = normalize_sample_id(r.get("sample") or r.get("Sample") or "")
    if sample:
        iqtree_excluded_map[sample] = r

def get_iqtree_inclusion_info(sample):
    sample_id = normalize_sample_id(sample)

    if sample_id in iqtree_excluded_map:
        r = iqtree_excluded_map[sample_id]
        note = clean(r.get("exclusion_note"))
        reason = clean(r.get("reason"))
        missing_fraction = clean(r.get("missing_fraction"))
        threshold = clean(r.get("threshold"))

        if not note:
            try:
                missing_percent = f"{float(missing_fraction) * 100:.2f}%"
            except Exception:
                missing_percent = missing_fraction or "not reported"

            try:
                threshold_percent = f"{float(threshold) * 100:.0f}%" if threshold else "50%"
            except Exception:
                threshold_percent = threshold or "50%"

            reason_low = reason.lower()

            if reason_low == "no_usable_acgt_bases":
                note = (
                    f"Excluded from IQ-TREE because the core-SNP alignment had "
                    f"{missing_percent} missing/ambiguous/gap content and no usable ACGT bases."
                )
            elif reason_low.startswith("missing_fraction_ge_"):
                note = (
                    f"Excluded from IQ-TREE because {missing_percent} of the core-SNP alignment "
                    f"was missing/ambiguous/gap content, exceeding the allowed threshold of {threshold_percent}."
                )
            elif reason_low == "empty_sequence":
                note = "Excluded from IQ-TREE because the core-SNP alignment sequence was empty."
            else:
                note = "Excluded from IQ-TREE because it did not meet tree-building alignment-quality requirements."

        return {
            "status": "NO",
            "note": note,
            "reason": reason,
            "row": r
        }

    if sample_id in iqtree_included_map:
        return {
            "status": "YES",
            "note": "Retained after IQ-TREE core-SNP alignment quality filtering.",
            "reason": "included_after_alignment_filtering",
            "row": iqtree_included_map[sample_id]
        }

    return {
        "status": "Not reported",
        "note": "No IQ-TREE inclusion/exclusion record was available for this sample.",
        "reason": "not_reported",
        "row": {}
    }

def iqtree_inclusion_html(sample):
    info = get_iqtree_inclusion_info(sample)
    status = info["status"]
    note = info["note"]

    if status == "YES":
        return green_badge("YES") + f"<br><small>{safe(note)}</small>"

    if status == "NO":
        return red_badge("NO") + f"<br><small>{safe(note)}</small>"

    return neutral_badge("Not reported") + f"<br><small>{safe(note)}</small>"

def selected_for_mtbc_workflow_html(mtbc_selected, mtbc_reason):
    selected = clean(mtbc_selected).upper()

    if selected == "YES":
        return green_badge("YES") + f"<br><small>{safe(mtbc_reason)}</small>"

    if selected == "NO":
        return red_badge("NO") + f"<br><small>{safe(mtbc_reason)}</small>"

    return neutral_badge("Not reported") + f"<br><small>{safe(mtbc_reason)}</small>"

iqtree_included_count = len([
    s for s in iqtree_included_map
    if s and not is_reference_name(s) and (not routing_by_sample or s in mtbc_routed_samples)
])

iqtree_excluded_count = len([
    s for s in iqtree_excluded_map
    if s and not is_reference_name(s) and (not routing_by_sample or s in mtbc_routed_samples)
])

def normalize_drug_name(x):
    s = lower_clean(x)
    s = re.sub(r"[_\-]+", " ", s)
    s = re.sub(r"\s+", " ", s)

    aliases = {
        "inh": "isoniazid",
        "isoniazid": "isoniazid",
        "h": "isoniazid",
        "rif": "rifampicin",
        "rmp": "rifampicin",
        "rifampin": "rifampicin",
        "rifampicin": "rifampicin",
        "r": "rifampicin",
        "pza": "pyrazinamide",
        "pyrazinamide": "pyrazinamide",
        "z": "pyrazinamide",
        "emb": "ethambutol",
        "ethambutol": "ethambutol",
        "e": "ethambutol",
        "sm": "streptomycin",
        "str": "streptomycin",
        "streptomycin": "streptomycin",
        "s": "streptomycin",
        "levo": "levofloxacin",
        "levofloxacin": "levofloxacin",
        "lfx": "levofloxacin",
        "moxi": "moxifloxacin",
        "moxifloxacin": "moxifloxacin",
        "mfx": "moxifloxacin",
        "ofx": "ofloxacin",
        "ofloxacin": "ofloxacin",
        "amikacin": "amikacin",
        "amk": "amikacin",
        "kanamycin": "kanamycin",
        "kan": "kanamycin",
        "capreomycin": "capreomycin",
        "cap": "capreomycin",
        "bedaquiline": "bedaquiline",
        "bdq": "bedaquiline",
        "linezolid": "linezolid",
        "lzd": "linezolid",
        "clofazimine": "clofazimine",
        "cfz": "clofazimine",
        "ethionamide": "ethionamide",
        "eto": "ethionamide",
        "prothionamide": "prothionamide",
        "pto": "prothionamide",
        "cycloserine": "cycloserine",
        "cs": "cycloserine",
        "para aminosalicylic acid": "para-aminosalicylic acid",
        "para-aminosalicylic-acid": "para-aminosalicylic acid",
        "pas": "para-aminosalicylic acid",
        "delamanid": "delamanid",
        "dlm": "delamanid",
        "pretomanid": "pretomanid",
        "pa": "pretomanid",
    }

    return aliases.get(s, s)

def split_drug_tokens(value):
    if not value:
        return []

    txt = clean(value)
    txt = txt.replace(";", ",")
    txt = txt.replace("|", ",")
    txt = txt.replace("/", ",")
    txt = re.sub(r"\band\b", ",", txt, flags=re.IGNORECASE)

    out = []

    for part in txt.split(","):
        p = normalize_drug_name(part)
        if p and p not in {
            "none",
            "none reported",
            "not reported",
            "susceptible",
            "sensitive",
            "unknown",
            "na",
            "n/a",
            "not determined",
            "resistance not determined by tb-profiler",
            "no resistance detected by tb-profiler"
        }:
            out.append(p)

    return out

def normalize_resistance_profile(value):
    raw = clean(value)
    low = raw.lower().replace("_", "-")

    if low in {"", "none", "none reported", "not reported", "unknown", "na", "n/a", "not determined"}:
        return "Unknown"

    if "resistance not determined" in low:
        return "Resistance not determined by TB-Profiler"

    if low in {"sensitive", "susceptible"} or "no resistance detected" in low:
        return "Sensitive"

    if "xdr" in low and "pre" not in low:
        return "XDR-TB"

    if "pre-xdr" in low or "pre xdr" in low or "prexdr" in low:
        return "Pre-XDR-TB"

    # Preserve the umbrella only when no drug-level calls are available.
    # Drug-level classification elsewhere converts rif-only to RR-TB and rif+INH to MDR-TB.
    if re.search(r"\bmdr\s*/\s*rr[- ]?tb\b", low) or re.search(r"\bmdr[- ]?rr[- ]?tb\b", low):
        return "MDR/RR-TB"

    if "mdr" in low:
        return "MDR-TB"

    if re.search(r"\brr[- ]?tb\b", low):
        return "RR-TB"

    if re.search(r"\bhr[- ]?tb\b", low) or "isoniazid-resistant" in low:
        return "Hr-TB"

    if "mono" in low:
        return "Monoresistance"

    if "poly" in low:
        return "Polyresistance"

    if "other drug resistance" in low:
        return "Other drug resistance"

    return "Unknown"


def classify_resistance_fallback(dr_type, resistant_drugs):
    """
    Final report resistance classifier.

    Critical fix for RR-TB versus MDR/RR-TB:
    The raw TB-Profiler dr_type can contain the umbrella label MDR/RR-TB.
    That umbrella should not override the actual resistant_drugs list. When
    drug-level calls are available, classify strictly from those drugs.
    """
    drugs = set(split_drug_tokens(resistant_drugs))

    if drugs:
        isoniazid = "isoniazid" in drugs
        rifampicin = "rifampicin" in drugs

        fluoroquinolones = {
            "levofloxacin",
            "moxifloxacin",
            "ofloxacin",
            "gatifloxacin",
            "ciprofloxacin"
        }

        group_a_additional = {
            "bedaquiline",
            "linezolid"
        }

        has_fq = bool(drugs.intersection(fluoroquinolones))
        has_group_a_additional = bool(drugs.intersection(group_a_additional))

        if rifampicin and has_fq and has_group_a_additional:
            return "XDR-TB"

        if rifampicin and has_fq:
            return "Pre-XDR-TB"

        if rifampicin and isoniazid:
            return "MDR-TB"

        if rifampicin and not isoniazid:
            return "RR-TB"

        if isoniazid and not rifampicin:
            return "Hr-TB"

        if len(drugs) == 1:
            return "Monoresistance"

        return "Polyresistance"

    # Use dr_type only if no resistant_drugs could be parsed.
    return normalize_resistance_profile(dr_type)


def build_resistance_profile_map():
    profile_map = {}

    for r in resistance_rows:
        sample = (
            r.get("sample_id") or
            r.get("sample") or
            r.get("Sample ID") or
            ""
        )

        sample_id = normalize_sample_id(sample)

        if not sample_id:
            continue

        profile_raw = (
            r.get("resistance_profile") or
            r.get("dr_type") or
            r.get("Resistance profile") or
            ""
        )

        resistant_drugs = (
            r.get("resistant_drugs") or
            r.get("Predicted resistant drugs") or
            r.get("drug_resistance") or
            ""
        )

        resistance_mutations = (
            r.get("resistance_mutations") or
            r.get("Resistance-associated mutations") or
            ""
        )

        key_mutations = (
            r.get("key_mutations") or
            r.get("All key mutations") or
            ""
        )

        status = r.get("status") or ""

        # Always classify using resistant_drugs first. This prevents the TB-Profiler
        # umbrella label MDR/RR-TB from masking true RR-TB when resistant_drugs
        # contains rifampicin but not isoniazid.
        profile = classify_resistance_fallback(profile_raw, resistant_drugs)

        profile_map[sample_id] = {
            "resistance_profile": profile,
            "resistant_drugs": normalize_missing(resistant_drugs, "None reported"),
            "resistance_mutations": normalize_missing(resistance_mutations, "None reported"),
            "key_mutations": normalize_missing(key_mutations, "None reported"),
            "status": status
        }

    return profile_map

resistance_profile_map = build_resistance_profile_map()

def get_resistance_info(sample, fallback_row=None):
    sample_id = normalize_sample_id(sample)

    if sample_id in resistance_profile_map:
        return resistance_profile_map[sample_id]

    fallback_row = fallback_row or {}
    fallback_profile = classify_resistance_fallback(
        fallback_row.get("dr_type") or fallback_row.get("resistance_profile") or "",
        fallback_row.get("resistant_drugs") or ""
    )

    return {
        "resistance_profile": fallback_profile,
        "resistant_drugs": normalize_missing(fallback_row.get("resistant_drugs"), "None reported"),
        "resistance_mutations": normalize_missing(fallback_row.get("resistance_mutations"), "None reported"),
        "key_mutations": normalize_missing(fallback_row.get("key_mutations"), "None reported"),
        "status": fallback_row.get("status") or ""
    }

def is_drug_resistant_category(category):
    return category not in {
        "Sensitive",
        "Resistance not determined by TB-Profiler",
        "Unknown",
        ""
    }

def badge(label):
    label = normalize_resistance_profile(label)
    color = RESISTANCE_COLORS.get(label, RESISTANCE_COLORS["Unknown"])
    text_color = "#ffffff"
    return f'<span class="res-badge" style="background:{color};color:{text_color} !important;">{safe(label)}</span>'

def green_badge(label):
    return f'<span class="badge badge-green" style="background:#28A745 !important;color:white !important;font-weight:700;">{safe(label)}</span>'

def red_badge(label):
    return f'<span class="badge badge-red" style="background:#b91c1c !important;color:white !important;font-weight:700;">{safe(label)}</span>'

def neutral_badge(label):
    return f'<span class="badge" style="background:#9ca3af !important;color:white !important;font-weight:700;">{safe(label)}</span>'

def species_badge(label):
    text = lower_clean(label)

    if text == "mycobacterium tuberculosis":
        return (
            '<span class="badge" '
            'style="background:#c2e3f4 !important;color:#111111 !important;font-weight:700;">'
            f'{safe(label)}</span>'
        )

    return green_badge(label)

def cluster_badge(label):
    text = lower_clean(label)

    if "genomically close" in text or text == "close" or "likely recent" in text or "likely" in text:
        return '<span class="badge" style="background:#22c55e !important;color:white !important;font-weight:700;">Genomically close; review epidemiological linkage</span>'

    if "intermediate" in text or "possible" in text:
        return '<span class="badge" style="background:#facc15 !important;color:#111111 !important;font-weight:700;">Intermediate SNP distance; review metadata</span>'

    if "not clustered" in text or "distant" in text:
        return '<span class="badge" style="background:#38bdf8 !important;color:#111111 !important;font-weight:700;">Not clustered by SNP threshold</span>'

    return '<span class="badge" style="background:#9ca3af !important;color:white !important;font-weight:700;">Not interpreted</span>'

total_samples = max(len(species_rows), len(routing_by_sample), len(rows))

# Dashboard counters are driven by routing, not by downstream task artifacts.
mtbc_retained = len([s for s in mtbc_routed_samples if s and not is_reference_name(s)])
non_mtbc = len([s for s in non_mtbc_routed_samples if s and not is_reference_name(s)])

if not routing_by_sample:
    # Backward-compatible fallback for older reports without router TSV.
    mtbc_retained = sum(1 for r in rows if clean(r.get("mtbc_selected")).upper() == "YES")
    non_mtbc = max(total_samples - mtbc_retained, 0)

drug_resistant = 0

for r in rows:
    sample = r.get("sample")
    if routing_by_sample and not is_sample_routed_mtbc(sample):
        continue
    info = get_resistance_info(sample, r)
    category = info["resistance_profile"]

    if is_drug_resistant_category(category):
        drug_resistant += 1

def build_qc_section():
    qc_content = read_optional_file(qc_html)
    trim_content = read_optional_file(trim_html)
    variant_content = read_optional_file(variant_html)

    sample_ids = []

    for r in species_rows:
        sid = r.get("Sample_ID") or r.get("sample") or ""
        if sid:
            sample_ids.append(normalize_sample_id(sid))

    if not sample_ids:
        sample_ids = [normalize_sample_id(r.get("sample", "")) for r in rows if r.get("sample")]

    if not sample_ids:
        body = '<tr><td colspan="5">No sample-level QC records available.</td></tr>'
    else:
        body = []

        for s in sample_ids:
            body.append(
                "<tr>"
                f"<td>{display_sample_id(s)}</td>"
                "<td>Reported in MultiQC</td>"
                "<td>See trimming report</td>"
                '<td><span class="badge badge-green" style="background:#28A745 !important;color:white !important;font-weight:700;">PASS</span></td>'
                '<td><span class="badge badge-green" style="background:#28A745 !important;color:white !important;font-weight:700;">Proceed</span></td>'
                "</tr>"
            )

        body = "".join(body)

    embedded = ""

    if qc_content:
        embedded += '<details><summary>Embedded QC summary report</summary><div class="embedded-report">' + qc_content + "</div></details>"

    if trim_content:
        embedded += '<details><summary>Embedded trimming report</summary><div class="embedded-report">' + trim_content + "</div></details>"

    if variant_content:
        embedded += '<details><summary>Embedded variant summary report</summary><div class="embedded-report">' + variant_content + "</div></details>"

    if not embedded:
        embedded = '<div class="note">QC, trimming, and variant reports were not provided as separate HTML inputs, but sample-level workflow decisions are summarized below.</div>'

    return f"""
<div class="section">
<h2>1. Sample QC and Trimming Summary</h2>
<div class="controls">
<input id="qcSearch" onkeyup="filterTable('qcSearch','qcTable')" placeholder="Search QC table...">
<button onclick="downloadCSV('qcTable','qc_summary.csv')">Download QC CSV</button>
</div>
<table id="qcTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('qcTable',0)">Sample ID</th>
<th class="status" onclick="sortTable('qcTable',1)">Raw reads</th>
<th class="status" onclick="sortTable('qcTable',2)">Trimmed reads</th>
<th class="status" onclick="sortTable('qcTable',3)">FastQC status</th>
<th class="status" onclick="sortTable('qcTable',4)">Workflow decision</th>
</tr>
</thead>
<tbody>
{body}
</tbody>
</table>
{embedded}
</div>
"""

def build_species_section():
    if not species_rows:
        body = '<tr><td colspan="3">No species typing results were generated.</td></tr>'
    else:
        body = []

        for r in species_rows:
            sample = r.get("Sample_ID") or r.get("sample") or ""
            species = r.get("Species_Identified") or r.get("species") or "No species-level Mycobacterium call"
            evidence = r.get("Evidence") or r.get("evidence") or "No supporting evidence available"

            body.append(
                "<tr>"
                f"<td>{safe(display_sample_id(sample))}</td>"
                f"<td>{species_badge(species)}</td>"
                f"<td>{safe(evidence)}</td>"
                "</tr>"
            )

        body = "".join(body)

    return f"""
<div class="section">
<h2>2. Species Typing using Kraken2 + Bracken</h2>
<div class="note">
Species typing was performed using Kraken2 against a custom Mycobacterium-only database embedded in the Docker image
<code>gmboowa/mycobacterium-kraken2-bracken:2026.05</code>. The table reports one most probable species-level call per sample based on the highest species-level Kraken2 assignment and supporting taxonomic evidence.
</div>

<div class="controls">
<input id="speciesSearch" onkeyup="filterTable('speciesSearch','speciesTable')" placeholder="Search species typing results...">
<button onclick="downloadCSV('speciesTable','species_typing_summary.csv')">Download Species Typing CSV</button>
</div>

<table id="speciesTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('speciesTable',0)">Sample ID</th>
<th class="species" onclick="sortTable('speciesTable',1)">Species Identified</th>
<th class="status" onclick="sortTable('speciesTable',2)">Evidence Supporting Call</th>
</tr>
</thead>
<tbody>
{body}
</tbody>
</table>
</div>
"""

def build_ntm_species_section():
    display_rows = ntm_display_rows_from_routing if ntm_display_rows_from_routing else [r for r in ntm_rows if clean(r.get("Mycobacteria_Category")) != "MTBC"]

    if not display_rows:
        body = '<tr><td colspan="7">No non-MTBC Mycobacterium or NTM samples were detected by Kraken2/Bracken.</td></tr>'
    else:
        body_parts = []
        for r in display_rows:
            sample = r.get("Sample_ID") or r.get("sample") or ""
            species = r.get("Most_Probable_Mycobacterium_Species") or r.get("Species_Identified") or r.get("species") or "Not resolved"
            category = r.get("Mycobacteria_Category") or "Non-MTBC Mycobacterium"
            mtbc_reads = r.get("MTBC_Reads") or "0"
            mtbc_percent = r.get("MTBC_Percent") or "0.00"
            evidence = r.get("Evidence") or "No supporting evidence reported"
            decision = r.get("Workflow_Decision") or "Reported and excluded from MTBC-specific downstream analyses"

            body_parts.append(
                "<tr>"
                f"<td>{safe(display_sample_id(sample))}</td>"
                f"<td><em>{safe(species)}</em></td>"
                f"<td>{safe(category)}</td>"
                f"<td>{safe(mtbc_reads)}</td>"
                f"<td>{safe(mtbc_percent)}</td>"
                f"<td>{safe(evidence)}</td>"
                f"<td>{safe(decision)}</td>"
                "</tr>"
            )
        body = "".join(body_parts)

    return f"""
<div class="section">
<h2>2b. Non-MTBC Mycobacteria Species Summary</h2>
<div class="note">
<strong>Routing status:</strong> {safe(routing_status_text)}<br><br>
Kraken2/Bracken species typing is used as the first routing layer. MTBC-supported samples are retained for TB-Profiler, drug-resistance interpretation, and MTBC-specific phylogenomics. Non-MTBC Mycobacterium samples are reported in this section and excluded from TB-Profiler and MTBC-specific downstream analyses.
</div>
<div class="controls">
<input id="ntmSearch" onkeyup="filterTable('ntmSearch','ntmTable')" placeholder="Search non-MTBC Mycobacteria results...">
<button onclick="downloadCSV('ntmTable','non_mtbc_mycobacteria_summary.csv')">Download Non-MTBC Mycobacteria CSV</button>
</div>
<table id="ntmTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('ntmTable',0)">Sample ID</th>
<th class="species" onclick="sortTable('ntmTable',1)">Most probable Mycobacterium species</th>
<th class="status" onclick="sortTable('ntmTable',2)">Category</th>
<th class="status" onclick="sortTable('ntmTable',3)">MTBC reads</th>
<th class="status" onclick="sortTable('ntmTable',4)">MTBC percent</th>
<th class="status" onclick="sortTable('ntmTable',5)">Evidence</th>
<th class="status" onclick="sortTable('ntmTable',6)">Workflow decision</th>
</tr>
</thead>
<tbody>
{body}
</tbody>
</table>
</div>
"""

def build_tb_rows():
    out = []

    for r in rows:
        sample = r.get("sample")
        if routing_by_sample and not is_sample_routed_mtbc(sample):
            continue
        species = r.get("species") or "Mycobacterium tuberculosis complex (inferred from available MTBC evidence)"
        main_lineage = normalize_missing(r.get("main_lineage"), "Not resolved by TB-Profiler")
        sub_lineage = normalize_missing(r.get("sub_lineage"), "Not resolved by TB-Profiler")
        lineage = f"{safe(main_lineage)} / {safe(sub_lineage)}"

        info = get_resistance_info(sample, r)
        category = info["resistance_profile"]
        resistant_drugs = info["resistant_drugs"]
        key_mutations = info["key_mutations"]

        if key_mutations and key_mutations != "None reported":
            resistance_detail = f"{resistant_drugs}<br><small><strong>Key mutations:</strong> {safe(key_mutations)}</small>"
        else:
            resistance_detail = safe(resistant_drugs)

        mtbc_selected = clean(r.get("mtbc_selected")).upper()
        mtbc_reason = r.get("mtbc_selection_reason") or ""

        selected_html = selected_for_mtbc_workflow_html(mtbc_selected, mtbc_reason)
        iqtree_html = iqtree_inclusion_html(sample)

        out.append(
            "<tr>"
            f"<td>{safe(display_sample_id(sample))}</td>"
            f"<td>{safe(species)}</td>"
            f"<td>{lineage}</td>"
            f"<td>{badge(category)}</td>"
            f"<td>{resistance_detail}</td>"
            f"<td>{selected_html}</td>"
            f"<td>{iqtree_html}</td>"
            "</tr>"
        )

    if not out:
        return '<tr><td colspan="7">No MTBC-routed TB-Profiler rows were available. Non-MTBC/NTM samples are intentionally excluded from this section.</td></tr>'

    return "".join(out)

def build_mutation_evidence_section():
    cleaned = [
        r for r in mutation_rows
        if any(clean(v) for v in r.values())
        and clean(r.get("sample"))
        and (not routing_by_sample or is_sample_routed_mtbc(r.get("sample")))
    ]

    if not cleaned:
        return """
<div class="section">
<h2>4. Resistance Mutation Evidence Summary</h2>
<div class="note">
No mutation-level TB-Profiler resistance evidence was generated or the task was skipped.
</div>
</div>
"""

    grouped = defaultdict(list)

    for r in cleaned:
        grouped[r.get("sample", "Unknown")].append(r)

    out = ["""
<div class="section">
<h2>4. Resistance Mutation Evidence Summary</h2>
<div class="note">
This section reports mutation-level drug-resistance evidence extracted from TB-Profiler JSON outputs. Results are grouped per sample to match the display style of the non-synonymous mutation summary. It complements the resistance profile by showing the underlying drug or evidence source, gene, mutation/change, confidence, and associated evidence fields where available.
</div>
"""]

    for s in sorted(grouped):
        out.append(f'<details><summary>Sample: {display_sample_id(s)} — {len(grouped[s])} mutation(s)</summary><table>')
        out.append("<thead><tr><th class='resistance'>Drug / Evidence source</th><th class='lineage'>Gene</th><th class='mutations'>Mutation</th><th class='mutations'>Change</th><th class='status'>Confidence</th><th class='status'>Evidence / associated drug(s)</th></tr></thead><tbody>")

        for r in grouped[s]:
            drug = normalize_missing(r.get("drug", ""), "Not reported")
            gene = normalize_missing(r.get("gene", ""), "Not reported")
            mutation = normalize_missing(r.get("mutation", ""), "Not reported")
            change = normalize_missing(r.get("change", ""), "Not reported")
            confidence = normalize_missing(r.get("confidence", ""), "Not reported")
            evidence = normalize_missing(r.get("evidence", ""), "Not reported")

            out.append(
                "<tr>"
                f"<td><span class='badge badge-red'>{safe(drug)}</span></td>"
                f"<td><span class='badge badge-blue'>{safe(gene)}</span></td>"
                f"<td>{safe(mutation)}</td>"
                f"<td>{safe(change)}</td>"
                f"<td><span class='badge badge-green'>{safe(confidence)}</span></td>"
                f"<td>{safe(evidence)}</td>"
                "</tr>"
            )

        out.append("</tbody></table></details>")

    out.append("</div>")
    return "".join(out)

def build_nonsyn_section():
    cleaned = [
        r for r in nonsyn_rows
        if any(clean(v) for v in r.values())
        and (not routing_by_sample or is_sample_routed_mtbc(r.get("sample")))
    ]

    if not cleaned:
        return """
<div class="section">
<h2>5. Non-synonymous Mutation Summary</h2>
<div class="note">
No mutations were detected or mutation analysis was skipped.
</div>
</div>
"""

    grouped = defaultdict(list)

    for r in cleaned:
        grouped[r.get("sample", "Unknown")].append(r)

    out = ["""
<div class="section">
<h2>5. Non-synonymous Mutation Summary</h2>
<div class="note">
<strong>Mechanism:</strong> mutations are grouped per sample from per-sample Snippy annotation outputs and filtered to configured TB drug-resistance-associated genes. This complements TB-Profiler and should not replace catalogue-based resistance interpretation.
</div>
"""]

    for s in sorted(grouped):
        out.append(f'<details><summary>Sample: {display_sample_id(s)} — {len(grouped[s])} mutation(s)</summary><table>')
        out.append("<thead><tr><th class='lineage'>Gene</th><th class='mutations'>Effect</th><th class='mutations'>AA change</th><th class='mutations'>NT change</th><th class='status'>Product</th></tr></thead><tbody>")

        for r in grouped[s]:
            out.append(
                "<tr>"
                f"<td><strong>{safe(r.get('gene'))}</strong></td>"
                f"<td>{safe(r.get('effect'))}</td>"
                f"<td>{safe(r.get('aa_change'))}</td>"
                f"<td>{safe(r.get('nt_change'))}</td>"
                f"<td>{safe(r.get('product'))}</td>"
                "</tr>"
            )

        out.append("</tbody></table></details>")

    out.append("</div>")
    return "".join(out)

def build_snp_distance_section():
    cleaned_pairs = [
        r for r in snp_pair_rows
        if clean(r.get("sample1"))
        and clean(r.get("sample2"))
        and not is_reference_name(r.get("sample1"))
        and not is_reference_name(r.get("sample2"))
    ]

    cleaned_clusters = [
        r for r in snp_cluster_rows
        if clean(r.get("cluster_id"))
        and clean(r.get("sample1"))
        and clean(r.get("sample2"))
        and not is_reference_name(r.get("sample1"))
        and not is_reference_name(r.get("sample2"))
    ]

    if not cleaned_pairs:
        return """
<div class="section">
<h2>7. Pairwise SNP Distance and Cluster Summary</h2>
<div class="note">
Pairwise SNP distance and cluster reporting was not generated or was skipped. This usually occurs when fewer than the required number of MTBC samples were available for phylogenomic analysis.
</div>
</div>
"""

    pair_body = []

    for r in cleaned_pairs:
        sample1 = r.get("sample1", "")
        sample2 = r.get("sample2", "")
        distance = r.get("snp_distance", "")
        comparable_sites = r.get("comparable_sites", "Not reported")
        interpretation = r.get("interpretation", "")
        cluster_class = r.get("cluster_class", "")

        pair_body.append(
            "<tr>"
            f"<td>{display_sample_id(sample1)}</td>"
            f"<td>{display_sample_id(sample2)}</td>"
            f"<td>{safe(distance)}</td>"
            f"<td>{safe(comparable_sites)}</td>"
            f"<td>{cluster_badge(interpretation or cluster_class)}</td>"
            f"<td>{safe(cluster_class)}</td>"
            "</tr>"
        )

    if cleaned_clusters:
        cluster_body = []

        for r in cleaned_clusters:
            cluster_body.append(
                "<tr>"
                f"<td>{safe(r.get('cluster_id'))}</td>"
                f"<td>{display_sample_id(r.get('sample1'))}</td>"
                f"<td>{display_sample_id(r.get('sample2'))}</td>"
                f"<td>{safe(r.get('snp_distance'))}</td>"
                f"<td>{cluster_badge(r.get('interpretation'))}</td>"
                "</tr>"
            )

        cluster_table = f"""
<h3>Genomically close sample pairs requiring epidemiological review</h3>
<table id="snpClusterTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('snpClusterTable',0)">Cluster ID</th>
<th class="sample" onclick="sortTable('snpClusterTable',1)">Sample 1</th>
<th class="sample" onclick="sortTable('snpClusterTable',2)">Sample 2</th>
<th class="mutations" onclick="sortTable('snpClusterTable',3)">SNP distance</th>
<th class="status" onclick="sortTable('snpClusterTable',4)">Interpretation</th>
</tr>
</thead>
<tbody>
{"".join(cluster_body)}
</tbody>
</table>
"""
    else:
        cluster_table = """
<div class="note">
No sample pairs met the configured SNP thresholds for genomically close or intermediate-distance review.
</div>
"""

    return f"""
<div class="section">
<h2>7. Pairwise SNP Distance and Cluster Summary</h2>
<div class="note">
Pairwise SNP distances were calculated from the MTBC core genome alignment after excluding reference/non-sample sequences. Interpretations are threshold-based and should be considered alongside epidemiological metadata, lineage, resistance profile, sequencing quality, and phylogenetic support.
</div>

<div class="controls">
<input id="snpSearch" onkeyup="filterTable('snpSearch','snpPairsTable')" placeholder="Search SNP distance pairs...">
<button onclick="downloadCSV('snpPairsTable','pairwise_snp_distance_pairs.csv')">Download SNP Distance CSV</button>
</div>

<table id="snpPairsTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('snpPairsTable',0)">Sample 1</th>
<th class="sample" onclick="sortTable('snpPairsTable',1)">Sample 2</th>
<th class="mutations" onclick="sortTable('snpPairsTable',2)">SNP distance</th>
<th class="status" onclick="sortTable('snpPairsTable',3)">Comparable sites</th>
<th class="status" onclick="sortTable('snpPairsTable',4)">Interpretation</th>
<th class="status" onclick="sortTable('snpPairsTable',5)">Cluster class</th>
</tr>
</thead>
<tbody>
{"".join(pair_body)}
</tbody>
</table>

{cluster_table}
</div>
"""

def display_sample_id(value):
    """Remove internal Terra shard/index prefixes from sample IDs shown in final report tables."""
    s = clean(value)
    return re.sub(r"^\d{3,6}[_-](?=ERR|SRR|DRR|SAM|ERS|SRS|[A-Za-z])", "", s)

def svg_file_to_inline_html(path, title="embedded SVG"):
    """Inline SVG content so Terra-downloaded final_report.html is self-contained."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    try:
        svg = Path(path).read_text(encoding="utf-8", errors="replace").strip()
        if "<svg" in svg.lower():
            return f'<div class="embedded-svg" role="img" aria-label="{safe(title)}">{svg}</div>'
    except Exception:
        return ""
    return ""

def lineage_color(index):
    palette = ["#2563eb", "#7c3aed", "#059669", "#dc2626", "#ea580c", "#0891b2", "#be185d", "#65a30d", "#9333ea", "#0f766e", "#b45309", "#4f46e5", "#c026d3", "#16a34a", "#e11d48"]
    return palette[index % len(palette)]

def build_lineage_bar_chart(cleaned_rows):
    valid = []
    for r in cleaned_rows:
        lineage = clean(r.get("lineage") or r.get("Lineage"))
        try:
            count = int(float(clean(r.get("count") or r.get("Count") or "0")))
        except Exception:
            continue
        if lineage and count > 0:
            valid.append((lineage, count))
    if not valid:
        return ""
    max_count = max(c for _, c in valid) or 1
    bars = []
    for idx, (lineage, count) in enumerate(valid):
        pct = max(4, int(round((count / max_count) * 100)))
        color = lineage_color(idx)
        bars.append(
            '<div class="lineage-bar-row">'
            f'<div class="lineage-label">{safe(lineage)}</div>'
            '<div class="lineage-bar-track">'
            f'<div class="lineage-bar" style="width:{pct}%; background:{color};"><span>{safe(count)}</span></div>'
            '</div></div>'
        )
    return '<div class="lineage-bar-chart">' + ''.join(bars) + '</div>'

def build_lineage_distribution_section():
    lineage_svg_exists = Path("final_report/lineage_distribution.svg").exists()

    # Recompute lineage counts from MTBC-routed TB-Profiler rows. This prevents
    # non-MTBC/NTM samples from appearing as L6 or Not resolved by TB-Profiler.
    lineage_counts = defaultdict(int)
    for rr in rows:
        sample = rr.get("sample")
        if routing_by_sample and not is_sample_routed_mtbc(sample):
            continue
        lineage = clean(rr.get("main_lineage")) or "Not resolved by TB-Profiler"
        if lower_clean(lineage) in ["", "none", "not reported", "unknown", "na", "n/a"]:
            lineage = "Not resolved by TB-Profiler"
        lineage_counts[lineage] += 1

    if lineage_counts:
        cleaned = [{"lineage": k, "count": str(v)} for k, v in sorted(lineage_counts.items())]
    else:
        cleaned = [
            r for r in lineage_rows
            if any(clean(v) for v in r.values())
        ]

    if not cleaned:
        return """
<div class="section">
<h2>6. Lineage Distribution Summary</h2>
<div class="note">Lineage distribution summary was not generated or was skipped.</div>
</div>
"""

    body = []

    for r in cleaned:
        body.append(
            "<tr>"
            f"<td>{safe(r.get('lineage'))}</td>"
            f"<td>{safe(r.get('count'))}</td>"
            "</tr>"
        )

    plot_html = build_lineage_bar_chart(cleaned)
    if not plot_html and lineage_svg_exists:
        plot_html = svg_file_to_inline_html("final_report/lineage_distribution.svg", "Lineage distribution")
    if not plot_html:
        plot_html = '<div class="note">Lineage distribution plot was not available.</div>'

    return f"""
<div class="section">
<h2>6. Lineage Distribution Summary</h2>
<div class="note">
This section summarizes TB-Profiler lineage calls where available. Samples supported as MTBC by Kraken2/Bracken but without TB-Profiler lineage resolution are shown as “Not resolved by TB-Profiler”.
</div>
{plot_html}
<div class="controls">
<input id="lineageSearch" onkeyup="filterTable('lineageSearch','lineageTable')" placeholder="Search lineage summary...">
<button onclick="downloadCSV('lineageTable','lineage_distribution.csv')">Download Lineage CSV</button>
</div>
<table id="lineageTable">
<thead>
<tr>
<th class="lineage" onclick="sortTable('lineageTable',0)">Lineage</th>
<th class="status" onclick="sortTable('lineageTable',1)">Count</th>
</tr>
</thead>
<tbody>
{''.join(body)}
</tbody>
</table>
</div>
"""

def build_snp_heatmap_section():
    heatmap_exists = Path("final_report/snp_distance_heatmap.svg").exists()

    heatmap_html = svg_file_to_inline_html("final_report/snp_distance_heatmap.svg", "SNP distance heatmap") if heatmap_exists else ""
    if not heatmap_html:
        heatmap_html = '<div class="note">SNP distance heatmap was not generated, was skipped, or the SVG was empty.</div>'

    return f"""
<div class="section">
<h2>8. SNP Distance Heatmap</h2>
<div class="note">
This heatmap visualizes pairwise SNP distances among MTBC isolates after excluding reference/non-sample sequences. Lower SNP distances indicate closer genomic relatedness and should be interpreted together with epidemiological metadata, lineage, resistance profile, and tree topology.
</div>
{heatmap_html}
</div>
"""

def build_qc_rationale_surveillance_metadata_section():
    qc_cleaned = [
        r for r in qc_rationale_rows
        if any(clean(v) for v in r.values())
    ]

    metadata_cleaned = [
        r for r in surveillance_metadata_rows
        if any(clean(v) for v in r.values())
    ]

    qc_body = []

    for r in qc_cleaned:
        sample = r.get("sample") or r.get("sample_id") or ""
        selected = r.get("selected_for_phylogeny", r.get("included", ""))
        selected_reason = r.get("reason") or ""

        selected_html = selected_for_mtbc_workflow_html(selected, selected_reason)
        iqtree_html = iqtree_inclusion_html(sample)

        qc_body.append(
            "<tr>"
            f"<td>{safe(display_sample_id(sample))}</td>"
            f"<td>{safe(r.get('mean_depth'))}</td>"
            f"<td>{safe(r.get('mtbc_percent'))}</td>"
            f"<td>{selected_html}</td>"
            f"<td>{iqtree_html}</td>"
            f"<td>{safe(selected_reason)}</td>"
            "</tr>"
        )

    metadata_body = []

    for r in metadata_cleaned:
        sample = r.get("sample") or r.get("sample_id") or ""
        info = get_resistance_info(sample, r)
        resistance_profile = info["resistance_profile"]
        resistant_drugs = info["resistant_drugs"]

        iqtree_html = iqtree_inclusion_html(sample)

        if is_drug_resistant_category(resistance_profile):
            drug_res_html = red_badge("YES")
        elif resistance_profile == "Sensitive":
            drug_res_html = green_badge("NO")
        else:
            drug_res_html = neutral_badge("Not determined")

        metadata_body.append(
            "<tr>"
            f"<td>{safe(display_sample_id(sample))}</td>"
            f"<td>{safe(r.get('integrated_mtbc_status') or r.get('tbprofiler_species'))}</td>"
            f"<td>{safe(r.get('mtbc_support_source'))}</td>"
            f"<td>{safe(r.get('tbprofiler_lineage_status'))}</td>"
            f"<td>{safe(r.get('kraken_species'))}</td>"
            f"<td>{safe(r.get('tbprofiler_main_lineage') or r.get('main_lineage'))}</td>"
            f"<td>{safe(r.get('tbprofiler_sub_lineage') or r.get('sub_lineage'))}</td>"
            f"<td>{safe(r.get('lineage_group'))}</td>"
            f"<td>{badge(resistance_profile)}</td>"
            f"<td>{drug_res_html}</td>"
            f"<td>{safe(resistant_drugs)}</td>"
            f"<td>{safe(r.get('mean_depth'))}</td>"
            f"<td>{iqtree_html}</td>"
            "</tr>"
        )

    if not qc_body:
        qc_body = ['<tr><td colspan="6">No QC filtering rationale records available.</td></tr>']

    if not metadata_body:
        metadata_body = ['<tr><td colspan="13">No surveillance metadata records available.</td></tr>']

    return f"""
<div class="section">
<h2>10. QC Filtering Rationale and Surveillance Metadata</h2>
<div class="note">
This section provides a transparent rationale for sample inclusion/exclusion and a surveillance-ready metadata table. Resistance profile, drug-resistance detected status, and resistant drugs are populated from the canonical <code>resistance_profile_summary.tsv</code> generated by TB-Profiler parsing, ensuring this table matches Section 3 and the phylogenetic tree labels.
</div>

<div class="controls">
<button onclick="downloadCSV('qcRationaleTable','qc_filtering_rationale.csv')">Download QC Rationale CSV</button>
<button onclick="downloadCSV('surveillanceMetadataTable','tb_surveillance_metadata.csv')">Download Surveillance Metadata CSV</button>
</div>

<h3>QC Filtering Rationale</h3>
<table id="qcRationaleTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('qcRationaleTable',0)">Sample</th>
<th class="status" onclick="sortTable('qcRationaleTable',1)">Mean depth</th>
<th class="status" onclick="sortTable('qcRationaleTable',2)">MTBC %</th>
<th class="status" onclick="sortTable('qcRationaleTable',3)">Selected for MTBC workflow</th>
<th class="status" onclick="sortTable('qcRationaleTable',4)">Included in IQ-TREE</th>
<th class="status" onclick="sortTable('qcRationaleTable',5)">Reason</th>
</tr>
</thead>
<tbody>
{''.join(qc_body)}
</tbody>
</table>

<h3>Surveillance Metadata</h3>
<table id="surveillanceMetadataTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('surveillanceMetadataTable',0)">Sample</th>
<th class="species" onclick="sortTable('surveillanceMetadataTable',1)">Integrated MTBC Status</th>
<th class="species" onclick="sortTable('surveillanceMetadataTable',2)">MTBC Support Source</th>
<th class="species" onclick="sortTable('surveillanceMetadataTable',3)">TB-Profiler Lineage Status</th>
<th class="species" onclick="sortTable('surveillanceMetadataTable',4)">Kraken Species</th>
<th class="lineage" onclick="sortTable('surveillanceMetadataTable',5)">TB-Profiler Main Lineage</th>
<th class="lineage" onclick="sortTable('surveillanceMetadataTable',6)">TB-Profiler Sub-lineage</th>
<th class="lineage" onclick="sortTable('surveillanceMetadataTable',7)">Lineage Group</th>
<th class="resistance" onclick="sortTable('surveillanceMetadataTable',8)">Resistance Profile</th>
<th class="resistance" onclick="sortTable('surveillanceMetadataTable',9)">Drug Resistance Detected</th>
<th class="mutations" onclick="sortTable('surveillanceMetadataTable',10)">Resistant Drugs</th>
<th class="status" onclick="sortTable('surveillanceMetadataTable',11)">Mean Depth</th>
<th class="status" onclick="sortTable('surveillanceMetadataTable',12)">Included in IQ-TREE</th>
</tr>
</thead>
<tbody>
{''.join(metadata_body)}
</tbody>
</table>
</div>
"""

def build_tree_legend():
    legend = ['<div class="legend">']

    for label in RESISTANCE_ORDER:
        color = RESISTANCE_COLORS[label]
        legend.append(f'<span><strong style="color:{color};">●</strong> {safe(label)}</span>')

    legend.append('<span><strong style="color:#b91c1c;">Bootstrap support (%)</strong></span>')
    legend.append('<span><strong>Scale bar</strong> substitutions/site</span>')
    legend.append("</div>")

    return "".join(legend)

def build_iqtree_section():
    txt = read_optional_file(iqtree_txt)

    if not txt:
        return ""

    return f"""
<details>
<summary>IQ-TREE report</summary>
<pre>{safe(txt[:30000])}</pre>
</details>
"""

def build_iqtree_exclusion_footnote():
    cleaned = [
        r for r in iqtree_excluded_rows
        if clean(r.get("sample"))
        and (not routing_by_sample or normalize_sample_id(r.get("sample")) in mtbc_routed_samples)
    ]

    included_count = len([
        r for r in iqtree_included_rows
        if clean(r.get("sample"))
        and (not routing_by_sample or normalize_sample_id(r.get("sample")) in mtbc_routed_samples)
    ])

    status_line = f"<p><strong>IQ-TREE status:</strong> {safe(iqtree_status_text)}"

    if included_count:
        status_line += f" &nbsp; | &nbsp; <strong>Samples retained for tree:</strong> {included_count}"

    status_line += "</p>"

    if not cleaned:
        return f"""
<div class="tree-notes tree-notes-pass">
<h4>Phylogenetic tree notes</h4>
{status_line}
<p>No samples were excluded by the IQ-TREE pre-filtering step. All samples that passed MTBC selection and core-alignment generation were eligible for phylogenetic inference.</p>
</div>
"""

    notes = []

    for r in cleaned:
        sample = clean(r.get("sample"))
        missing_fraction = clean(r.get("missing_fraction"))
        reason = clean(r.get("reason"))
        exclusion_note = clean(r.get("exclusion_note"))
        threshold = clean(r.get("threshold"))

        try:
            missing_percent = f"{float(missing_fraction) * 100:.2f}%"
        except Exception:
            missing_percent = missing_fraction or "not reported"

        try:
            threshold_percent = f"{float(threshold) * 100:.0f}%" if threshold else "50%"
        except Exception:
            threshold_percent = threshold or "50%"

        reason_low = reason.lower()

        if exclusion_note:
            note = exclusion_note
        elif reason_low == "no_usable_acgt_bases":
            note = (
                f"{sample} was excluded from IQ-TREE because it had "
                f"{missing_percent} missing/ambiguous/gap content in the "
                f"core-SNP alignment and no usable ACGT bases."
            )
        elif reason_low.startswith("missing_fraction_ge_"):
            if not threshold or threshold_percent == "50%":
                threshold_raw = reason_low.replace("missing_fraction_ge_", "")
                try:
                    threshold_percent = f"{float(threshold_raw) * 100:.0f}%"
                except Exception:
                    threshold_percent = threshold_percent

            note = (
                f"{sample} was excluded from IQ-TREE because "
                f"{missing_percent} of its core-SNP alignment was "
                f"missing/ambiguous/gap content, exceeding the maximum allowed "
                f"threshold of {threshold_percent}."
            )
        elif reason_low == "empty_sequence":
            note = (
                f"{sample} was excluded from IQ-TREE because its core-SNP "
                f"alignment sequence was empty."
            )
        else:
            note = (
                f"{sample} was excluded from IQ-TREE because it did not meet "
                f"the minimum alignment-quality requirements for phylogenetic inference."
            )

        if note.startswith(sample):
            note_remainder = note[len(sample):].lstrip()
            notes.append(f"<li><strong>{display_sample_id(sample)}</strong> {safe(note_remainder)}</li>")
        else:
            notes.append(f"<li>{safe(note)}</li>")

    summary_details = ""

    if iqtree_filtering_summary_text:
        summary_details = f"""
<details>
<summary>View IQ-TREE alignment filtering summary</summary>
<pre>{safe(iqtree_filtering_summary_text[:12000])}</pre>
</details>
"""

    return f"""
<div class="tree-notes">
<h4>Phylogenetic tree notes</h4>
{status_line}
<p><strong>{len(cleaned)} MTBC-routed sample(s)</strong> were excluded from IQ-TREE before phylogenetic inference because their core-SNP alignment sequence did not meet tree-building quality requirements. Non-MTBC/NTM samples are counted separately as excluded from the MTBC workflow and are not described as IQ-TREE-only exclusions.</p>
<ul>
{''.join(notes)}
</ul>
{summary_details}
</div>
"""

tree_image_path = Path("final_report/mtbc_tree.png")
if tree_image_path.exists() and tree_image_path.stat().st_size > 0:
    try:
        tree_b64 = base64.b64encode(tree_image_path.read_bytes()).decode("ascii")
        tree_html = (
            '<img class="mtbc-tree-img" '
            f'src="data:image/png;base64,{tree_b64}" '
            'alt="ETE3-rendered MTBC core-SNP phylogenetic tree">'
        )
    except Exception as exc:
        tree_html = (
            '<div class="note">Tree image was generated, but could not be embedded '
            f'into the final HTML report: {safe(str(exc))}</div>'
        )
else:
    tree_html = '<div class="note">Tree image was not provided or could not be copied into the final report folder.</div>'

html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>rMAP-TB Interactive Report</title>
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<style>
:root{{
  --bg:#f5f7fb;
  --card:#ffffff;
  --text:#1f2937;
  --muted:#6b7280;
  --border:#e5e7eb;
  --blue:#2563eb;
  --teal:#0f766e;
  --purple:#7c3aed;
  --red:#b91c1c;
  --orange:#d97706;
  --green:#087f5b;
  --dark:#12355b;
}}
body{{
  margin:0;
  font-family:Arial,Helvetica,sans-serif;
  background:var(--bg);
  color:var(--text);
}}
.header{{
  background:linear-gradient(135deg,#12355b,#0f766e);
  color:white;
  padding:28px 42px;
}}
.header h1{{
  margin:0;
  font-size:30px;
}}
.header p{{
  margin:8px 0 0;
  font-size:15px;
  opacity:.95;
}}
.container{{
  padding:28px 42px;
}}
.cards{{
  display:grid;
  grid-template-columns:repeat(4,minmax(150px,1fr));
  gap:16px;
  margin-bottom:24px;
}}
.card{{
  background:var(--card);
  border-radius:16px;
  padding:18px;
  box-shadow:0 2px 12px rgba(0,0,0,.08);
}}
.card h3{{
  margin:0;
  color:var(--muted);
  font-size:14px;
}}
.card .num{{
  font-size:30px;
  font-weight:bold;
  margin-top:8px;
}}
.blue{{color:var(--blue)}}
.green{{color:var(--green)}}
.orange{{color:var(--orange)}}
.red{{color:var(--red)}}
.section{{
  background:var(--card);
  border-radius:16px;
  padding:20px;
  margin-bottom:24px;
  box-shadow:0 2px 12px rgba(0,0,0,.08);
}}
.section h2{{
  margin-top:0;
  padding-bottom:10px;
  border-bottom:2px solid var(--border);
}}
.controls{{
  display:flex;
  flex-wrap:wrap;
  gap:10px;
  margin-bottom:14px;
}}
input,select,button{{
  border:1px solid var(--border);
  border-radius:10px;
  padding:9px 11px;
  font-size:13px;
  background:white;
}}
button{{
  cursor:pointer;
  background:#eef6ff;
  color:#12355b;
  font-weight:bold;
}}
button:hover{{
  background:#dbeafe;
}}
table{{
  width:100%;
  border-collapse:collapse;
  font-size:13px;
  overflow:hidden;
  border-radius:12px;
}}
th{{
  color:white;
  padding:10px;
  text-align:left;
  cursor:pointer;
  user-select:none;
}}
td{{
  padding:9px;
  border-bottom:1px solid var(--border);
  vertical-align:top;
}}
tr:hover td{{
  background:#f9fafb;
}}
th.sample{{background:#0f766e;}}
th.species{{background:#2563eb;}}
th.lineage{{background:#7c3aed;}}
th.resistance{{background:#b91c1c;}}
th.mutations{{background:#d97706;}}
th.status{{background:#087f5b;}}
.badge{{
  padding:4px 8px;
  border-radius:999px;
  color:#ffffff !important;
  font-size:12px;
  display:inline-block;
  font-weight:700 !important;
}}
.badge-green{{
  background:#28A745 !important;
  color:#ffffff !important;
  font-weight:700 !important;
}}
.badge-red{{
  background:#b91c1c !important;
  color:#ffffff !important;
}}
.badge-blue{{
  background:#2563eb !important;
  color:#ffffff !important;
}}
.badge-orange{{
  background:#d97706 !important;
  color:#ffffff !important;
}}
.res-badge{{
  padding:4px 8px;
  border-radius:999px;
  color:#ffffff !important;
  font-size:12px;
  display:inline-block;
  font-weight:bold;
}}
.note{{
  background:#eef6ff;
  border-left:5px solid #2563eb;
  padding:12px;
  border-radius:10px;
  margin:12px 0;
  line-height:1.45;
}}
.tree-notes{{
  margin-top:16px;
  padding:14px 16px;
  border-left:5px solid #f59e0b;
  background:#fffbeb;
  border-radius:12px;
  line-height:1.5;
}}
.tree-notes h4{{
  margin:0 0 8px 0;
  color:#92400e;
  font-size:16px;
}}
.tree-notes p{{
  margin:8px 0;
}}
.tree-notes ul{{
  margin:8px 0 0 20px;
  padding:0;
}}
.tree-notes li{{
  margin-bottom:8px;
}}
.tree-notes details{{
  margin-top:10px;
  background:#ffffff;
}}
.tree-notes-pass{{
  border-left-color:#28A745;
  background:#f0fdf4;
}}
.tree-notes-pass h4{{
  color:#166534;
}}
.grid2{{
  display:grid;
  grid-template-columns:minmax(0,1fr);
  gap:20px;
  align-items:start;
}}
.tree-panel{{
  background:#fbfdff;
  border:1px solid var(--border);
  border-radius:16px;
  padding:10px;
  overflow:auto;
  width:100%;
  box-sizing:border-box;
}}

.embedded-svg {{
  max-width:100%;
  overflow-x:auto;
  margin:1rem 0;
  border:1px solid var(--border);
  border-radius:14px;
  background:#fff;
  padding:1rem;
}}
.embedded-svg svg {{
  max-width:100%;
  height:auto;
}}
.lineage-bar-chart {{
  margin:1rem 0 1.25rem 0;
  padding:1rem;
  border:1px solid var(--border);
  border-radius:14px;
  background:#fff;
}}
.lineage-bar-row {{
  display:grid;
  grid-template-columns:minmax(120px,220px) 1fr;
  gap:.75rem;
  align-items:center;
  margin:.55rem 0;
}}
.lineage-label {{
  font-weight:700;
  color:#111827;
}}
.lineage-bar-track {{
  width:100%;
  background:#f3f4f6;
  border-radius:999px;
  overflow:hidden;
  min-height:1.8rem;
}}
.lineage-bar {{
  color:#fff;
  font-weight:700;
  min-height:1.8rem;
  display:flex;
  align-items:center;
  justify-content:flex-end;
  padding-right:.7rem;
  border-radius:999px;
}}
.mtbc-tree-img{{
  display:block;
  width:auto;
  max-width:none;
  height:auto;
  margin:0;
}}
.legend{{
  display:flex;
  flex-wrap:wrap;
  gap:10px;
  margin-top:10px;
  font-size:12px;
}}
.legend span{{
  border-radius:999px;
  padding:5px 9px;
  background:#f3f4f6;
}}
.details{{
  background:#fafafa;
  border:1px solid var(--border);
  border-radius:12px;
  padding:12px;
}}
details{{
  margin-top:12px;
  margin-bottom:10px;
  border:1px solid var(--border);
  border-radius:12px;
  padding:10px;
  background:#fff;
}}
summary{{
  cursor:pointer;
  font-weight:bold;
}}
pre{{
  white-space:pre-wrap;
  overflow:auto;
  background:#111827;
  color:#f9fafb;
  padding:14px;
  border-radius:12px;
  font-size:12px;
}}
.embedded-report{{
  overflow:auto;
  max-height:650px;
  border:1px solid var(--border);
  border-radius:12px;
  padding:10px;
  background:white;
}}
.footer{{
  font-size:12px;
  color:var(--muted);
  margin-top:20px;
}}
@media(max-width:900px){{
  .cards{{grid-template-columns:repeat(2,1fr);}}
  .container{{padding:18px;}}
  .header{{padding:22px;}}
}}
</style>
</head>

<body>
<div class="header">
<h1>rMAP-TB Interactive Report</h1>
<p>Trimming → QC → Species typing → MTBC/NTM routing → TB-Profiler for MTBC-supported samples → mutation evidence → lineage and surveillance summaries → SNP distance clustering and heatmap → core-SNP phylogenomics → final merged report</p>
<p><strong>Run generated:</strong> {safe(run_started_utc)} &nbsp; | &nbsp; <strong>Run stamp:</strong> {safe(run_stamp)}</p>
</div>

<div class="container">

<div class="cards">
<div class="card"><h3>Total paired samples</h3><div class="num blue">{total_samples}</div></div>
<div class="card"><h3>MTBC isolates retained</h3><div class="num green">{mtbc_retained}</div></div>
<div class="card"><h3>Non-MTBC excluded</h3><div class="num orange">{non_mtbc}</div></div>
<div class="card"><h3>Drug-resistant isolates</h3><div class="num red">{drug_resistant}</div></div>
</div>

{build_qc_section()}

{build_species_section()}

{build_ntm_species_section()}

<div class="section">
<h2>3. TB-Profiler Resistance, Species, and Lineage Report</h2>

<div class="note">
<strong>Interpretation note:</strong> Resistance profile classifications in this section are populated from the canonical <code>resistance_profile_summary.tsv</code> generated during TB-Profiler parsing. The same source is also used for Surveillance Metadata, tree labels, resistance badges, resistant drug lists, and drug-resistant isolate counts.
<br><br>
<strong>Reporting distinction:</strong> <em>Selected for MTBC workflow</em> means the sample was retained by the MTBC support logic and remains part of the wider workflow. <em>Included in IQ-TREE</em> means the sample also passed core-SNP alignment quality filtering and was actually eligible for IQ-TREE phylogenetic inference. Samples can therefore be selected for the MTBC workflow but excluded from IQ-TREE only.
<br><br>
<strong>WHO 2021+ resistance definitions:</strong>
<strong>Hr-TB:</strong> resistant to isoniazid and not resistant to rifampicin.
<strong>RR-TB:</strong> resistant to rifampicin without detected isoniazid resistance.
<strong>MDR-TB:</strong> resistant to at least both isoniazid and rifampicin.
<strong>MDR/RR-TB:</strong> umbrella term for MDR-TB or RR-TB; the final report separates this into <strong>RR-TB</strong> when rifampicin is detected without isoniazid, and <strong>MDR-TB</strong> when both rifampicin and isoniazid are detected.
<strong>Pre-XDR-TB:</strong> MDR/RR-TB that is also resistant to any fluoroquinolone.
<strong>XDR-TB:</strong> MDR/RR-TB that is resistant to any fluoroquinolone and at least one additional Group A drug, bedaquiline or linezolid.
</div>

<div class="controls">
<input id="tbSearch" onkeyup="filterTable('tbSearch','tbTable')" placeholder="Search TB-Profiler results...">
<select onchange="filterResistance(this.value)">
<option value="">All resistance profiles</option>
<option value="Sensitive">Sensitive only</option>
<option value="Hr-TB">Hr-TB only</option>
<option value="RR-TB">RR-TB only</option>
<option value="MDR-TB">MDR-TB only</option>
<option value="MDR/RR-TB">MDR/RR-TB umbrella only</option>
<option value="Pre-XDR-TB">Pre-XDR-TB only</option>
<option value="XDR-TB">XDR-TB only</option>
<option value="Other drug resistance">Other drug resistance only</option>
<option value="Resistance not determined by TB-Profiler">Resistance not determined only</option>
<option value="Unknown">Unknown only</option>
</select>
<button onclick="downloadCSV('tbTable','tbprofiler_summary.csv')">Download TB-Profiler CSV</button>
</div>

<table id="tbTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('tbTable',0)">Sample ID</th>
<th class="species" onclick="sortTable('tbTable',1)">Species</th>
<th class="lineage" onclick="sortTable('tbTable',2)">Lineage</th>
<th class="resistance" onclick="sortTable('tbTable',3)">Resistance profile</th>
<th class="mutations" onclick="sortTable('tbTable',4)">Resistant drugs / key mutations</th>
<th class="status" onclick="sortTable('tbTable',5)">Selected for MTBC workflow</th>
<th class="status" onclick="sortTable('tbTable',6)">Included in IQ-TREE</th>
</tr>
</thead>
<tbody>
{build_tb_rows()}
</tbody>
</table>
</div>

{build_mutation_evidence_section()}

{build_nonsyn_section()}

{build_lineage_distribution_section()}

{build_snp_distance_section()}

{build_snp_heatmap_section()}

<div class="section">
<h2>9. MTBC-only Core-SNP Phylogenetic Tree</h2>
<div class="grid2">
<div class="tree-panel">
{tree_html}
{build_tree_legend()}
{build_iqtree_exclusion_footnote()}
</div>
<div class="details">
<h3>Tree construction summary</h3>
<p><strong>Selected for MTBC workflow:</strong> {mtbc_retained} MTBC isolate(s) retained in the wider workflow.</p>
<p><strong>Included in IQ-TREE:</strong> {iqtree_included_count} non-reference MTBC isolate(s) retained after core-SNP alignment quality filtering.</p>
<p><strong>Excluded from IQ-TREE only:</strong> {iqtree_excluded_count} sample(s) excluded from phylogenetic inference because of alignment-quality issues.</p>
<p><strong>Excluded from MTBC workflow:</strong> {non_mtbc} non-MTBC or low-confidence isolate(s).</p>
<p><strong>Core alignment:</strong> Snippy-core alignment.</p>
<p><strong>Recombination:</strong> Optional Gubbins-filtered alignment when enabled.</p>
<p><strong>Tree:</strong> IQ-TREE2 maximum-likelihood phylogeny.</p>
<p><strong>Display:</strong> ETE3-rendered static tree image, shown inside an auto-scaling scrollable report panel.</p>
{build_iqtree_section()}
</div>
</div>
</div>

{build_qc_rationale_surveillance_metadata_section()}

<div class="section">
<h2>11. Pipeline Provenance and Software Versions</h2>
<p>The report documents all samples through QC, species typing, TB-Profiler analysis, mutation-level resistance evidence, lineage distribution, SNP distance clustering, SNP heatmap visualization, surveillance metadata, and MTBC-only phylogenomic reconstruction. Samples not classified as MTBC are excluded from the tree but retained in the workflow record for transparency. Samples with excessive missing, ambiguous, or gap-only content in the core-SNP alignment may also be excluded from IQ-TREE phylogenetic inference and are listed under the tree footnotes.</p>
<div class="note">
<strong>Interpretation:</strong> use close clustering together with bootstrap support, lineage, drug-resistance profile, mutation-level resistance evidence, lineage distribution, surveillance metadata, sample-exclusion notes, and SNP distances before making transmission inferences.
</div>
<table>
<thead>
<tr>
<th class="sample">Workflow component</th>
<th class="status">Description</th>
</tr>
</thead>
<tbody>
<tr><td>Species typing</td><td>Kraken2 + Bracken using <code>gmboowa/mycobacterium-kraken2-bracken:2026.05</code></td></tr>
<tr><td>TB resistance and lineage</td><td>TB-Profiler Docker image provided by workflow input</td></tr>
<tr><td>Canonical resistance profile</td><td><code>resistance_profile_summary.tsv</code> used for Section 3, Surveillance Metadata, tree labels, badges, resistant drugs, and drug-resistant isolate count</td></tr>
<tr><td>Mutation-level resistance evidence</td><td>Parsed from TB-Profiler JSON outputs and summarized by sample, drug or evidence source, gene, mutation/change, confidence, and evidence</td></tr>
<tr><td>Lineage distribution</td><td>Lineage counts and barplot generated from TB-Profiler lineage fields where resolved</td></tr>
<tr><td>Pairwise SNP distance and clustering</td><td>Pairwise SNP distances calculated from the MTBC core genome alignment after excluding reference/non-sample sequences and interpreted using configured SNP thresholds</td></tr>
<tr><td>SNP distance heatmap</td><td>SVG heatmap generated from the pairwise SNP distance matrix after reference filtering</td></tr>
<tr><td>Surveillance metadata</td><td>Downloadable metadata and QC filtering rationale TSV files generated for transparent surveillance reporting</td></tr>
<tr><td>Core-SNP phylogenomics</td><td>Snippy-core, optional Gubbins filtering, IQ-TREE2, and ETE3 tree rendering</td></tr>
<tr><td>IQ-TREE problematic-sample filtering</td><td>Samples with excessive missing, ambiguous, or gap-only sequence content in the core-SNP alignment are excluded from IQ-TREE and reported in <code>excluded_from_iqtree.tsv</code></td></tr>
<tr><td>Report generated</td><td>{safe(run_started_utc)}</td></tr>
<tr><td>Run stamp</td><td>{safe(run_stamp)}</td></tr>
</tbody>
</table>
</div>

<div class="footer">
Generated by rMAP-TB WDL workflow. Run generated: {safe(run_started_utc)}. Run stamp: {safe(run_stamp)}.
</div>

</div>

<script>
function filterTable(inputId, tableId){{
  const filter=document.getElementById(inputId).value.toLowerCase();
  const rows=document.getElementById(tableId).getElementsByTagName("tbody")[0].rows;
  for(let i=0;i<rows.length;i++){{
    rows[i].style.display=rows[i].innerText.toLowerCase().includes(filter)?"":"none";
  }}
}}

function filterResistance(value){{
  const rows=document.getElementById("tbTable").getElementsByTagName("tbody")[0].rows;
  for(let i=0;i<rows.length;i++){{
    rows[i].style.display=value===""||rows[i].cells[3].innerText.includes(value)?"":"none";
  }}
}}

function sortTable(tableId,col){{
  const table=document.getElementById(tableId);
  const tbody=table.tBodies[0];
  const rows=Array.from(tbody.rows);
  const asc=table.getAttribute("data-sort-col")!=col||table.getAttribute("data-sort-dir")!=="asc";
  rows.sort((a,b)=>{{
    const A=a.cells[col].innerText.replace(/,/g,'');
    const B=b.cells[col].innerText.replace(/,/g,'');
    const nA=parseFloat(A),nB=parseFloat(B);
    if(!isNaN(nA)&&!isNaN(nB))return asc?nA-nB:nB-nA;
    return asc?A.localeCompare(B):B.localeCompare(A);
  }});
  rows.forEach(r=>tbody.appendChild(r));
  table.setAttribute("data-sort-col",col);
  table.setAttribute("data-sort-dir",asc?"asc":"desc");
}}

function downloadCSV(tableId, filename){{
  const table=document.getElementById(tableId);
  let csv=[];
  for(const row of table.rows){{
    let cols=[];
    for(const cell of row.cells){{
      cols.push('"' + cell.innerText.replace(/"/g,'""') + '"');
    }}
    csv.push(cols.join(","));
  }}
  const blob=new Blob([csv.join("\\n")],{{type:"text/csv"}});
  const link=document.createElement("a");
  link.href=URL.createObjectURL(blob);
  link.download=filename;
  link.click();
}}
</script>

</body>
</html>
"""

(outdir / "integrated_tb_amr_mtbc_phylogenomics_report.html").write_text(html_out, encoding="utf-8")
PY

    py_rc=$?

    if [ "$py_rc" -ne 0 ]; then
      echo "WARNING: MERGE_TB_REPORTS Python report builder exited with code ${py_rc}; creating fallback report." >&2
      cat > final_report/integrated_tb_amr_mtbc_phylogenomics_report.html <<'EOF_REPORT'
<!doctype html>
<html>
<head><meta charset="utf-8"><title>rMAP-TB Report</title></head>
<body style="font-family:Arial,Helvetica,sans-serif;margin:24px;">
<h1>rMAP-TB Integrated Report</h1>
<p>The final merge report builder encountered an error, but upstream workflow outputs were generated. Please inspect the task logs for details.</p>
</body>
</html>
EOF_REPORT
      echo "Fallback report generated after merge error code ${py_rc}" > final_report/run_metadata.txt
    fi

    if [ ! -f final_report/lineage_distribution.svg ]; then
      cat > final_report/lineage_distribution.svg <<'EOF_SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240">
<rect width="100%" height="100%" fill="#ffffff"/>
<text x="450" y="110" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">Lineage Distribution</text>
<text x="450" y="145" text-anchor="middle" font-family="Arial" font-size="14" fill="#475569">Lineage distribution SVG was not generated.</text>
</svg>
EOF_SVG
    fi

    if [ ! -f final_report/snp_distance_heatmap.svg ]; then
      cat > final_report/snp_distance_heatmap.svg <<'EOF_SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240">
<rect width="100%" height="100%" fill="#ffffff"/>
<text x="450" y="110" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">SNP Distance Heatmap</text>
<text x="450" y="145" text-anchor="middle" font-family="Arial" font-size="14" fill="#475569">SNP distance heatmap SVG was not generated.</text>
</svg>
EOF_SVG
    fi

    [ -f final_report/tb_surveillance_metadata.tsv ] || echo -e "sample\tintegrated_mtbc_status\tmtbc_support_source\ttbprofiler_lineage_status\tkraken_species\ttbprofiler_main_lineage\ttbprofiler_sub_lineage\tlineage_group\tresistance_profile\tdrug_resistance_detected\tresistant_drugs\tmean_depth\tincluded_in_tree" > final_report/tb_surveillance_metadata.tsv
    [ -f final_report/qc_filtering_rationale.tsv ] || echo -e "sample\tmean_depth\tmtbc_percent\tselected_for_phylogeny\tincluded_in_tree\treason" > final_report/qc_filtering_rationale.tsv
    [ -f final_report/resistance_profile_summary.tsv ] || echo -e "sample_id\tresistance_profile\tresistant_drugs\tresistance_mutations\tkey_mutations\tstatus" > final_report/resistance_profile_summary.tsv
    [ -f final_report/excluded_from_iqtree.tsv ] || echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note" > final_report/excluded_from_iqtree.tsv
    [ -f final_report/included_in_iqtree.tsv ] || echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction" > final_report/included_in_iqtree.tsv
    [ -f final_report/alignment_filtering_summary.txt ] || echo "IQ-TREE alignment filtering summary was not provided." > final_report/alignment_filtering_summary.txt
    [ -f final_report/iqtree_status.txt ] || echo "unknown_status" > final_report/iqtree_status.txt
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 50 HDD"
    timeout: "12 hours"
  }

  output {
    File run_metadata = "final_report/run_metadata.txt"
    File final_report_html = "final_report/integrated_tb_amr_mtbc_phylogenomics_report.html"
    File? merged_lineage_distribution_svg = "final_report/lineage_distribution.svg"
    File? merged_snp_distance_heatmap_svg = "final_report/snp_distance_heatmap.svg"
    File? merged_surveillance_metadata_tsv = "final_report/tb_surveillance_metadata.tsv"
    File? merged_qc_filtering_rationale_tsv = "final_report/qc_filtering_rationale.tsv"
    File? merged_resistance_profile_summary_tsv = "final_report/resistance_profile_summary.tsv"
    File? merged_iqtree_excluded_samples_tsv = "final_report/excluded_from_iqtree.tsv"
    File? merged_iqtree_included_samples_tsv = "final_report/included_in_iqtree.tsv"
    File? merged_iqtree_filtering_summary_txt = "final_report/alignment_filtering_summary.txt"
    File? merged_iqtree_status_txt = "final_report/iqtree_status.txt"
  }
}
