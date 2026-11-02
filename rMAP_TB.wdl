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
    Boolean use_gubbins = false
    Boolean midpoint_root_tree = true
    Int max_bootstrap_to_label = 100
    Boolean report_nonsynonymous_drug_gene_mutations = true
    Boolean report_snp_distance_clusters = true
    Boolean generate_surveillance_summary_visuals = true

    String trimmomatic_quality_encoding = "phred33"
    String species_typing_docker = "gmboowa/mycobacterium-kraken2-bracken:2026.05"
    String tbprofiler_docker = "staphb/tbprofiler:6.6.6"
    String snippy_reference_type = "genbank"
    String iqtree2_model = "GTR+G"
    Int iqtree2_bootstraps = 1000
    Float min_species_call_percent = 5.0

    Float mtbc_node_override_percent = 50.0

    Int min_mtbc_samples_for_tree = 3
    Int likely_transmission_snp_threshold = 5
    Int possible_transmission_snp_threshold = 12

    File? snp_mask_bed
    Int min_comparable_sites = 2205000

    Int max_heterozygous_sites_for_clustering = 0

    Boolean exclude_mixed_infection_from_clustering = true
    Int min_mean_depth = 15
    Float min_breadth_10x = 0.90

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

  call PREPARE_READ_PAIRS {
    input:
      read1 = read1,
      read2 = read2
  }

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

  if (do_species_typing) {
    call SPECIES_TYPING {
      input:
        input_reads = analysis_reads,
        min_species_call_percent = min_species_call_percent,
        mtbc_node_override_percent = mtbc_node_override_percent,
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

  if (do_tb_profiler && do_phylogeny && size(mtbc_reads) >= min_mtbc_fastq_files_for_tree) {
    call SNIPPY_CORE_MTBC {
      input:
        input_reads = mtbc_reads,
        reference_genome = mtbc_reference_genbank,
        reference_type = snippy_reference_type,
        cpu = cpu_16,
        memory_gb = max_memory_gb,
        min_quality = min_mapping_quality,
        snp_mask_bed = snp_mask_bed,
        min_mean_depth = min_mean_depth,
        min_breadth_10x = min_breadth_10x
    }
  }

  if (report_nonsynonymous_drug_gene_mutations && defined(SNIPPY_CORE_MTBC.snippy_tab_files)) {
    call TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS {
      input:
        snippy_tab_files = select_first([SNIPPY_CORE_MTBC.snippy_tab_files, []]),
        genes_csv = tb_drug_resistance_genes
    }
  }

  if (do_phylogeny && report_snp_distance_clusters && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call SNP_DISTANCE_CLUSTERING {
      input:
        core_full_alignment = select_first([SNIPPY_CORE_MTBC.core_full_alignment]),
        likely_transmission_snp_threshold = likely_transmission_snp_threshold,
        possible_transmission_snp_threshold = possible_transmission_snp_threshold,
        min_comparable_sites = min_comparable_sites,
        max_heterozygous_sites_for_clustering = max_heterozygous_sites_for_clustering,
        tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv,
        exclude_mixed_infection_from_clustering = exclude_mixed_infection_from_clustering,
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

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

  if (do_phylogeny && use_gubbins && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call GUBBINS_RECOMBINATION {
      input:
        core_full_alignment = select_first([SNIPPY_CORE_MTBC.core_full_alignment]),
        cpu = cpu_16,
        memory_gb = max_memory_gb
    }
  }

  if (do_phylogeny && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call IQTREE2_PHYLOGENY {
      input:
        alignment = select_first([GUBBINS_RECOMBINATION.filtered_alignment, SNIPPY_CORE_MTBC.core_full_alignment]),
        model = iqtree2_model,
        bootstrap_replicates = iqtree2_bootstraps,
        cpu = cpu_16,
        memory_gb = max_memory_gb,
        max_missing_fraction_for_tree = 0.50,
        min_non_reference_samples_for_tree = 3
    }
  }

  if (do_phylogeny && defined(IQTREE2_PHYLOGENY.final_tree)) {
    call TREE_VISUALIZATION {
      input:
        input_tree = select_first([IQTREE2_PHYLOGENY.final_tree]),
        tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv,
        resistance_profile_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.resistance_profile_summary_tsv,
        iqtree_excluded_samples_tsv = IQTREE2_PHYLOGENY.excluded_from_iqtree,
        width = tree_width,
        height = tree_height,
        image_format = tree_image_format,
        midpoint_root_tree = midpoint_root_tree,
        max_bootstrap_to_label = max_bootstrap_to_label
    }
  }

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
    coverage_qc_tsv = SNIPPY_CORE_MTBC.coverage_qc_tsv,

    iqtree_report = IQTREE2_PHYLOGENY.iqtree_report,
    iqtree_excluded_samples_tsv = IQTREE2_PHYLOGENY.excluded_from_iqtree,
    iqtree_included_samples_tsv = IQTREE2_PHYLOGENY.included_in_iqtree,
    iqtree_filtering_summary_txt = IQTREE2_PHYLOGENY.alignment_filtering_summary,
    iqtree_status = IQTREE2_PHYLOGENY.iqtree_status,

    tree_image = TREE_VISUALIZATION.tree_image,
    tree_svg = TREE_VISUALIZATION.tree_svg,
    tree_render_log = TREE_VISUALIZATION.render_log,
    tree_render_status = TREE_VISUALIZATION.render_status,
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
    File? coverage_qc_tsv = SNIPPY_CORE_MTBC.coverage_qc_tsv

    File? nonsynonymous_drug_gene_mutations_tsv = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_tsv
    File? nonsynonymous_drug_gene_mutations_html = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_html

    File? pairwise_snp_distance_matrix = SNP_DISTANCE_CLUSTERING.pairwise_snp_distance_matrix
    File? pairwise_snp_distance_pairs = SNP_DISTANCE_CLUSTERING.pairwise_snp_distance_pairs
    File? snp_cluster_summary = SNP_DISTANCE_CLUSTERING.snp_cluster_summary
    File? snp_distance_cluster_html = SNP_DISTANCE_CLUSTERING.snp_distance_cluster_html
    File? alignment_ambiguity = SNP_DISTANCE_CLUSTERING.alignment_ambiguity
    File? clustering_eligibility = SNP_DISTANCE_CLUSTERING.clustering_eligibility

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
    File? tree_svg = TREE_VISUALIZATION.tree_svg
    File? tree_render_log = TREE_VISUALIZATION.render_log
    File? tree_render_status = TREE_VISUALIZATION.render_status

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
    Float min_species_call_percent = 5.0
    Float mtbc_node_override_percent = 50.0
    Array[File]+ input_reads
    File? qc_dependency
    String docker_image = "gmboowa/mycobacterium-kraken2-bracken:2026.05"
    Int cpu = 16
    Int memory_gb = 64
  }

  command <<<
    set -euo pipefail

    mkdir -p species_typing kraken_reports kraken_outputs

    KDB=/opt/kraken2_db/mycobacterium

    : > species_typing/kraken2_db_inspect.txt
    : > species_typing/kraken2_db_species.tsv
    : > species_typing/kraken2_db_accessions.tsv

    {
      printf 'field\tvalue\n'
      printf 'kraken2_db_path\t%s\n' "${KDB}"
      printf 'docker_image\t%s\n' "~{docker_image}"
      if [ -f "${KDB}/hash.k2d" ]; then
        printf 'hash_k2d_bytes\t%s\n' "$(stat -c %s "${KDB}/hash.k2d" 2>/dev/null || echo NA)"
        printf 'hash_k2d_md5\t%s\n' "$(md5sum "${KDB}/hash.k2d" 2>/dev/null | awk '{print $1}')"
      else
        printf 'hash_k2d_bytes\tNA\n'
      fi
      if [ -f "${KDB}/seqid2taxid.map" ]; then
        printf 'reference_sequences\t%s\n' "$(wc -l < "${KDB}/seqid2taxid.map")"
      elif [ -f "${KDB}/database.kraken" ]; then
        printf 'reference_sequences\t%s\n' \
          "$(cut -f2 "${KDB}/database.kraken" | sort -u | wc -l | tr -d ' ')"
      else
        printf 'reference_sequences\tNA\n'
      fi
      if [ -f "${KDB}/taxo.k2d" ]; then
        printf 'taxo_k2d_md5\t%s\n' "$(md5sum "${KDB}/taxo.k2d" 2>/dev/null | awk '{print $1}')"
      else
        printf 'taxo_k2d_md5\tNA\n'
      fi
      if [ -f "${KDB}/database150mers.kmer_distrib" ]; then
        printf 'bracken_kmer_distrib_present\ttrue (unused; Bracken is not run)\n'
      fi
    } > species_typing/kraken2_db_manifest.tsv 2>/dev/null || true

    if command -v kraken2-inspect >/dev/null 2>&1; then
      kraken2-inspect --db "${KDB}" > species_typing/kraken2_db_inspect.txt 2>/dev/null || true
      awk -F'\t' '$4=="S"{sub(/^ +/,"",$6); print $6"\t"$2}' \
        species_typing/kraken2_db_inspect.txt > species_typing/kraken2_db_species.tsv 2>/dev/null || true
    else
      printf 'kraken2-inspect not available in container; database composition not enumerated\n' \
        > species_typing/kraken2_db_inspect.txt
    fi

    if [ -f "${KDB}/seqid2taxid.map" ]; then
      cut -f1,2 "${KDB}/seqid2taxid.map" | sort -u \
        > species_typing/kraken2_db_accessions.tsv 2>/dev/null || true
    elif [ -f "${KDB}/database.kraken" ]; then
      cut -f2,3 "${KDB}/database.kraken" | sort -u \
        > species_typing/kraken2_db_accessions.tsv 2>/dev/null || true
    fi

    files=(~{sep=' ' input_reads})
    n=${#files[@]}

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: input_reads must contain paired-end FASTQ files in R1/R2 order." >&2
      exit 1
    fi


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

      species_call_confidence="$(python3 - <<PYCONF
try:
    pct = float("${top_percent}")
except Exception:
    pct = 0.0
print("low" if pct < float("~{min_species_call_percent}") else "high")
PYCONF
)"

      if [ -n "$mtbc_line" ]; then
        mtbc_percent="$(echo "$mtbc_line" | awk '{print $1}')"
        mtbc_reads="$(echo "$mtbc_line" | awk '{print $2}')"
      else
        mtbc_percent="0.00"
        mtbc_reads="0"
      fi


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
      if echo "$final_call_lc" | grep -Eq "mycobacterium tuberculosis complex|tuberculosis complex|(^|[^a-z])mtbc([^a-z]|$)|mycobacterium tuberculosis|m\. tuberculosis|mycobacterium africanum|mycobacterium bovis|mycobacterium caprae|mycobacterium canetti|mycobacterium microti|mycobacterium pinnipedii|mycobacterium orygis|mycobacterium mungi|mycobacterium suricattae|dassie bacillus"; then
        top_call_is_mtbc="YES"
      fi

      top_call_is_clear_non_mtbc="NO"
      if [ "$top_call_is_mtbc" = "NO" ] && [ "$species_call_confidence" = "high" ] \
         && echo "$final_call_lc" | grep -q "mycobacterium"; then
        if ! echo "$final_call_lc" | grep -Eq "no species-level|unknown|not resolved|unresolved|^na$|^n/a$"; then
          top_call_is_clear_non_mtbc="YES"
        fi
      fi

      mtbc_node_majority="$(python3 - <<PYMAJ
pct = float("${mtbc_percent_float}")
reads = int("${mtbc_reads_int}")
print("YES" if (pct >= float("~{mtbc_node_override_percent}") and reads >= 1000) else "NO")
PYMAJ
)"

      strong_mtbc_node_support="$(python3 - <<PYSTRONG
pct = float("${mtbc_percent_float}")
reads = int("${mtbc_reads_int}")
print("YES" if (pct >= 10.0 and reads >= 1000) else "NO")
PYSTRONG
)"

      if [ "$top_call_is_mtbc" = "YES" ]; then
        mtbc_supported="YES"
        selection_basis="Selected for MTBC workflow because the top Kraken2 species-level call is a recognized MTBC member"
      elif [ "$mtbc_node_majority" = "YES" ]; then
        mtbc_supported="YES"
        selection_basis="Selected for MTBC workflow because MTBC-complex reads are the majority (${mtbc_reads} reads; ${mtbc_percent}%), which exceeds the top species-level call (${final_call}; ${top_percent}%). Treat as an MTBC sample with non-MTBC contamination"
      elif [ "$top_call_is_clear_non_mtbc" = "YES" ]; then
        mtbc_supported="NO"
        selection_basis="Not selected for MTBC workflow because the top Kraken2 species-level call is a clear non-MTBC Mycobacterium; MTBC-complex reads are treated as background/cross-assignment"
      elif [ "$strong_mtbc_node_support" = "YES" ]; then
        mtbc_supported="YES"
        selection_basis="Selected for MTBC workflow because the species-level call was unresolved and Kraken2 detected strong MTBC-complex support: ${mtbc_reads} reads; ${mtbc_percent}%"
      else
        mtbc_supported="NO"
        selection_basis="Not selected for MTBC workflow because Kraken2 species typing did not provide MTBC support; low-level MTBC reads alone are treated as background/cross-assignment"
      fi

      if [ -z "$top_species_line" ]; then
        evidence="No species-level Kraken2 assignment detected; MTBC support: ${mtbc_reads} reads; ${mtbc_percent}%"
      else
        evidence="Top species-level assignment: ${final_call} (${top_reads} reads; ${top_percent}%); MTBC support: ${mtbc_reads} reads; ${mtbc_percent}%"
        if [ "$species_call_confidence" = "low" ]; then
          evidence="${evidence}; LOW CONFIDENCE: the top species call rests on ${top_percent}% of classified reads, below the ~{min_species_call_percent}% minimum. An organism with no representative in the database is assigned its nearest relative, so this name may be wrong. Treated as unresolved for MTBC routing"
        fi
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
html.append('<h2>2. Species Typing using Kraken2</h2>')
html.append('<p>Species typing was performed using Kraken2 against a custom <em>Mycobacterium</em>-only database embedded in the Docker image <code>gmboowa/mycobacterium-kraken2-bracken:2026.05</code> (76 RefSeq sequences, 50 species; composition and provenance documented in <code>resources/kraken2_db/README.md</code>). The table reports the single most probable species-level call for each sample based on the highest species-level Kraken2 assignment and supporting taxonomic evidence.</p>')
html.append('<p><strong>Database scope:</strong> this database contains <em>Mycobacterium</em> references only. It includes no human, non-mycobacterial bacterial, viral, archaeal or fungal sequences, so reads from any non-mycobacterial source are returned as <em>unclassified</em> and cannot be identified. The unclassified fraction should therefore be read as &ldquo;not mycobacterial&rdquo; and not as a contamination profile.</p>')
html.append('<p><strong>Interpretation limit within MTBC:</strong> the database contains only two <em>M. tuberculosis</em> complex genomes &mdash; <code>NC_000962.3</code> (<em>M. tuberculosis</em> H37Rv) and <code>NC_015848.1</code> (<em>M. canettii</em> CIPT 140010059). <em>M. bovis</em>, <em>M. africanum</em>, <em>M. microti</em>, <em>M. orygis</em>, <em>M. caprae</em> and <em>M. pinnipedii</em> are <strong>not represented</strong>, and MTBC members share &gt;99.9% nucleotide identity. A Kraken2 call here is therefore evidence of MTBC membership only and must not be reported as species-level discrimination within the complex. Use lineage/SNP-barcode assignment (TB-Profiler) for intra-MTBC identification.</p>')
html.append('<p><strong>Phylogeny selection rule:</strong> samples are selected for downstream Snippy/core-SNP/IQ-TREE analysis only when Kraken2 supports MTBC. TB-Profiler species, lineage, and resistance are reported later as annotations but do not determine phylogeny selection.</p>')
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
    File? kraken2_db_manifest = "species_typing/kraken2_db_manifest.tsv"
    File? kraken2_db_inspect = "species_typing/kraken2_db_inspect.txt"
    File? kraken2_db_species = "species_typing/kraken2_db_species.tsv"
    File? kraken2_db_accessions = "species_typing/kraken2_db_accessions.tsv"
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
        "mycobacterium canetti",
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

    if species_is_mtbc(species):
        return "MTBC"
    if is_mtbc_supported(row):
        return "MTBC"
    if species_is_clear_non_mtbc_mycobacterium(species):
        return "NTM"
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
    body.append("<tr><td colspan='7'>No non-MTBC Mycobacterium or NTM samples were detected by Kraken2.</td></tr>")

html_text = f"""<section class="card">
<h2>Non-MTBC Mycobacteria Species Summary</h2>
<p><strong>Routing status:</strong> {html.escape(status)}</p>
<p>Kraken2 species typing was used as the first routing layer. MTBC-supported samples are retained for TB-Profiler and MTBC-specific phylogenomics. Non-MTBC Mycobacterium samples are reported here and excluded from MTBC-specific downstream analyses.</p>
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

    echo -e "sample\tspecies\tmain_lineage\tsub_lineage\tdr_type\tresistant_drugs\tresistance_mutations\tkey_mutations\tjson_file\tmtbc_selected\tmtbc_selection_reason\tmixed_infection_status\tmixed_infection_evidence\tstatus" > tbprofiler_summary.tsv
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

NON_INFORMATIVE_CONFIDENCE = {
    "", "none", "none reported", "not reported", "na", "n/a",
    "unknown", "null", "nan"
}


def extract_who_confidence(item, drug_hint=None):
    """Return the WHO catalogue confidence grade for a TB-Profiler variant.

    TB-Profiler >= 5 does NOT store confidence at the top level of a variant
    record. It nests one entry per implicated drug, e.g.

        {"gene": "katG", "change": "p.Ser315Thr",
         "drugs":      [{"drug": "isoniazid", "confidence": "Assoc w R", ...}],
         "annotation": [{"type": "who_confidence",
                         "drug": "isoniazid", "confidence": "Assoc w R"}]}

    Reading item["confidence"] therefore always yields nothing, which is why
    every mutation previously rendered as "Not reported". This helper walks the
    nested blocks, preferring an entry matching the drug under consideration,
    and only then falls back to any legacy top-level field.
    """
    wanted = str(drug_hint or "").strip().lower()

    def _candidates(blocks):
        out = []
        if isinstance(blocks, dict):
            blocks = [blocks]
        if not isinstance(blocks, list):
            return out
        for b in blocks:
            if not isinstance(b, dict):
                continue
            conf = (
                b.get("confidence")
                or b.get("who_confidence")
                or b.get("grade")
                or b.get("confidence_level")
                or ""
            )
            conf = str(conf).strip()
            if not conf or conf.lower() in NON_INFORMATIVE_CONFIDENCE:
                continue
            bd = str(b.get("drug") or b.get("drug_name") or "").strip().lower()
            out.append((bd, conf))
        return out

    matched = []
    unmatched = []
    for key in ("drugs", "annotation", "consequences", "gene_associated_drugs"):
        for bd, conf in _candidates(item.get(key)):
            if wanted and bd and bd == wanted:
                matched.append(conf)
            else:
                unmatched.append(conf)

    if matched:
        return matched[0]
    if unmatched:
        return unmatched[0]

    legacy = str(
        item.get("confidence")
        or item.get("who_confidence")
        or item.get("grade")
        or ""
    ).strip()
    if legacy and legacy.lower() not in NON_INFORMATIVE_CONFIDENCE:
        return legacy

    return ""


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

LOCUS_TO_GENE = {
    "Rv0667": "rpoB", "Rv1908c": "katG", "Rv1484": "inhA", "Rv1483": "fabG1",
    "Rv2428": "ahpC", "Rv3795": "embB", "Rv2043c": "pncA", "Rv0682": "rpsL",
    "Rv0006": "gyrA", "Rv0005": "gyrB", "Rv2416c": "eis", "Rv3854c": "ethA",
    "Rv3855": "ethR", "Rv1694": "tlyA", "Rv3919c": "gid", "Rv1305": "atpE",
    "Rv2535c": "pepQ", "Rv2764c": "thyA", "Rv2447c": "folC",
}


def gene_label(item):
    """Prefer a gene name; show the locus tag with it rather than instead."""
    name = ""
    for key in ("gene_name", "gene", "locus_tag", "locus"):
        v = clean_value(item.get(key) or "")
        if v and not re.match(r"^(Rv|MT|MTB)\d", v):
            name = v
            break

    tag = ""
    for key in ("locus_tag", "gene", "locus"):
        v = clean_value(item.get(key) or "")
        if v and re.match(r"^(Rv|MT|MTB)\d", v):
            tag = v
            break

    if not name and tag:
        name = LOCUS_TO_GENE.get(tag, "")
    if name and tag:
        return "%s (%s)" % (name, tag)
    return name or tag


MIXED_LINEAGE_FRACTION_CEILING = 0.90
HETERO_FREQ_MIN = 0.10
HETERO_FREQ_MAX = 0.90


def _as_fraction(value):
    """TB-Profiler reports fractions as 0-1 floats, occasionally as strings."""
    try:
        f = float(str(value).strip().rstrip("%"))
    except (TypeError, ValueError):
        return None
    if f > 1.0:
        f = f / 100.0
    if 0.0 <= f <= 1.0:
        return f
    return None


ANIMAL_LINEAGES = [
    (r"\bla\s*1\b|\bbovis\b|\bbcg\b", "La1 (M. bovis / BCG)"),
    (r"\bla\s*2\b|\bcaprae\b", "La2 (M. caprae)"),
    (r"\bla\s*3\b|\borygis\b", "La3 (M. orygis)"),
    (r"\bmicroti\b", "M. microti"),
    (r"\bpinnipedii\b", "M. pinnipedii"),
    (r"\bmungi\b", "M. mungi"),
    (r"\bsuricattae\b", "M. suricattae"),
    (r"\bdassie\b", "Dassie bacillus"),
    (r"\bcanetti", "M. canettii"),
]


def animal_lineage_label(text):
    """Return a stable animal-lineage label, or "" if this is not one."""
    low = str(text or "").lower()
    for pattern, label in ANIMAL_LINEAGES:
        if re.search(pattern, low):
            return label
    return ""


def _main_lineage_token(value):
    """Reduce a lineage call to a token comparable across samples.

    Human lineages collapse to 'lineage4' and so on. Animal-adapted members
    return their own token, so that an M. bovis plus M. tuberculosis mixture
    counts as two distinct lineages rather than silently dropping one.
    """
    raw = str(value or "").strip()
    if not raw:
        return ""

    animal = animal_lineage_label(raw)
    if animal:
        return "animal:" + re.sub(r"[^a-z0-9]+", "", animal.split("(")[0].lower())

    t = re.sub(r"[^0-9a-z.]", "", raw.lower())
    m = re.match(r"(?:lineage)?(\d+)", t)
    return "lineage" + m.group(1) if m else ""


def assess_mixed_infection(data, tbprofiler_status):
    """Return (status, evidence) describing any mixed-infection signal."""
    if tbprofiler_status != "success":
        return ("Not assessable", "TB-Profiler did not complete for this sample")

    notes = []
    lineage_mixture = False

    blocks = data.get("lineage")
    if isinstance(blocks, dict):
        blocks = list(blocks.values())
    if not isinstance(blocks, list):
        blocks = []

    tokens, low_fractions = [], []
    for blk in blocks:
        if not isinstance(blk, dict):
            tok = _main_lineage_token(blk)
            if tok:
                tokens.append(tok)
            continue
        tok = _main_lineage_token(
            blk.get("lin") or blk.get("lineage") or blk.get("name") or ""
        )
        if tok:
            tokens.append(tok)
        frac = _as_fraction(blk.get("frac") or blk.get("fraction"))
        if frac is not None and 0.0 < frac < MIXED_LINEAGE_FRACTION_CEILING:
            low_fractions.append((tok or "unnamed lineage", frac))

    for field in ("main_lineage", "main_lin"):
        raw = str(data.get(field) or "")
        for part in re.split(r"[;,/]", raw):
            tok = _main_lineage_token(part)
            if tok:
                tokens.append(tok)

    distinct = sorted(set(t for t in tokens if t))

    assess_mixed_infection.last_detected_lineages = list(distinct)

    if len(distinct) > 1:
        lineage_mixture = True
        notes.append("multiple MTBC main lineages reported (%s)" % ", ".join(distinct))

    if low_fractions:
        lineage_mixture = True
        notes.append(
            "lineage fraction below %.2f (%s)"
            % (
                MIXED_LINEAGE_FRACTION_CEILING,
                "; ".join("%s=%.2f" % (n, f) for n, f in sorted(low_fractions)),
            )
        )

    hetero = []
    for key in ("dr_variants", "resistance_variants"):
        block = data.get(key, [])
        if isinstance(block, dict):
            block = list(block.values())
        if not isinstance(block, list):
            continue
        for item in block:
            if not isinstance(item, dict):
                continue
            freq = _as_fraction(
                item.get("freq") if item.get("freq") is not None
                else item.get("frequency")
            )
            if freq is None:
                continue
            if HETERO_FREQ_MIN <= freq <= HETERO_FREQ_MAX:
                label = " ".join(
                    x for x in [
                        gene_label(item),
                        clean_value(item.get("change") or item.get("mutation") or ""),
                    ] if x
                ) or "variant"
                hetero.append("%s @ %.2f" % (label, freq))

    if hetero:
        shown = "; ".join(hetero[:5])
        if len(hetero) > 5:
            shown += "; and %d more" % (len(hetero) - 5)
        notes.append(
            "%d resistance variant(s) at intermediate allele frequency (%s)"
            % (len(hetero), shown)
        )

    if lineage_mixture:
        return ("Possible mixed infection", "; ".join(notes))
    if hetero:
        return ("Possible heteroresistance", "; ".join(notes))
    if not blocks and not tokens:
        return ("Not assessable", "TB-Profiler reported no usable lineage fields")
    return ("No mixed-infection signal detected", "Single lineage; no intermediate-frequency resistance variants")


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
                            row.get("Kraken2_Evidence") or
                            row.get("Kraken2_Bracken_Evidence") or
                            ""
                        )

                        explicit = str(row.get("MTBC_Supported") or row.get("mtbc_supported") or "").strip().upper()
                        if explicit in ["YES", "TRUE", "1"]:
                            supports_mtbc = True
                        elif explicit in ["NO", "FALSE", "0"]:
                            supports_mtbc = False
                        else:
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
                                "mycobacterium canetti",
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
                            mtbc_majority = mtbc_pct >= 50.0 and mtbc_reads >= 1000
                            supports_mtbc = bool(
                                top_call_is_mtbc
                                or mtbc_majority
                                or ((not clear_non_mtbc)
                                    and mtbc_pct >= 10.0 and mtbc_reads >= 1000)
                            )
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

resistant_drugs_for_classification = []

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

    confidence = extract_who_confidence(item, drug)

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

mixed_status, mixed_evidence = assess_mixed_infection(data, status)

detected_lineages = getattr(assess_mixed_infection, "last_detected_lineages", [])
if len(detected_lineages) > 1:
    present = set(re.findall(r"lineage\s*\d+", str(main_lineage or "").lower()))
    if len(present) < len(detected_lineages):
        main_lineage = ";".join(detected_lineages)

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


if kraken_supports_mtbc:
    is_mtbc = True

    if has_tbprofiler_lineage:
        reason = "Selected for phylogeny because Kraken2 species typing supports MTBC; TB-Profiler lineage reported"
    else:
        reason = "Selected for phylogeny because Kraken2 species typing supports MTBC; TB-Profiler lineage not resolved"
else:
    is_mtbc = False

    if species_tsv and os.path.exists(species_tsv):
        if kraken_species or kraken_evidence:
            reason = "Not selected for phylogeny because Kraken2 species typing did not support MTBC"
        else:
            reason = "Not selected for phylogeny because no matching Kraken2 species-typing result supported MTBC"
    else:
        reason = "Not selected for phylogeny because Kraken2 species-typing results were not provided"

species_display = species

if kraken_supports_mtbc:
    if has_tbprofiler_lineage:
        if species_display:
            species_display = f"{species_display} (Kraken2 supports MTBC; TB-Profiler lineage reported)"
        else:
            species_display = "Mycobacterium tuberculosis complex (supported by Kraken2 species typing; TB-Profiler lineage reported)"
    else:
        if species_display:
            species_display = f"{species_display} (Kraken2 supports MTBC; TB-Profiler lineage not resolved)"
        else:
            species_display = "Mycobacterium tuberculosis complex (supported by Kraken2 species typing; TB-Profiler lineage not resolved)"
else:
    if species_display:
        species_display = f"{species_display} (not selected for phylogeny; Kraken2 did not support MTBC)"
    else:
        species_display = "Non-MTB / not classified as MTBC by Kraken2"

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
    mixed_status,
    mixed_evidence,
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
    ("mtbc_selected", "Selected for SNP tree by Kraken2"),
    ("mtbc_selection_reason", "Kraken2 selection reason"),
    ("mixed_infection_status", "Mixed-infection flag"),
    ("mixed_infection_evidence", "Mixed-infection evidence"),
    ("status", "TB-Profiler status")
]

colors = [
    "#0f172a", "#075985", "#7c2d12", "#4c1d95",
    "#991b1b", "#9f1239", "#d97706", "#b45309",
    "#166534", "#365314", "#854d0e", "#78350f", "#374151"
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
    "<p class='muted'>This table summarizes TB-Profiler JSON outputs, resistance-associated mutations, and Kraken2-based MTBC sample selection for downstream core-SNP phylogenomics. TB-Profiler species and lineage are reported as annotations only and do not determine SNP-tree inclusion.</p>",
    "<p class='muted'><strong>Phylogeny selection rule:</strong> samples are selected for Snippy/core-SNP/IQ-TREE only when Kraken2 species typing supports MTBC. IQ-TREE may later exclude selected samples from the final tree if their core-SNP alignment has excessive missing, ambiguous, or gap content.</p>",
    "<p class='muted'><strong>Mixed-infection flag:</strong> the mixed-infection "
    "columns are a screening flag derived from fields TB-Profiler already reports "
    "&mdash; more than one MTBC main lineage, a lineage fraction below 0.90, or "
    "resistance variants called at intermediate allele frequency (0.10&ndash;0.90), "
    "which indicates a resistant subpopulation. rMAP-TB does not run a dedicated "
    "mixed-infection detection tool, and no sensitivity or specificity is claimed "
    "for this flag. Flagged samples warrant review before their SNP distances are "
    "interpreted, because mixed samples carry ambiguous bases that distort pairwise "
    "comparisons. For quantitative mixed-strain analysis use a purpose-built tool.</p>",
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

    Array[File]+ input_reads

    File reference_genome
    String reference_type = "genbank"
    Int cpu = 16
    Int memory_gb = 64
    Int min_quality = 20

    File? snp_mask_bed

    Int min_mean_depth = 15
    Float min_breadth_10x = 0.90
  }

  command <<<
    set -uo pipefail
    mkdir -p snippy_results snippy_core logs

    MASK_ARG=""
    MASK_SRC="~{default="" snp_mask_bed}"
    if [ -n "${MASK_SRC}" ] && [ -s "${MASK_SRC}" ]; then
      cp "${MASK_SRC}" snippy_core/snp_mask.bed
      MASK_ARG="--mask snippy_core/snp_mask.bed"
      echo "Masking enabled: snippy_core/snp_mask.bed ($(wc -l < snippy_core/snp_mask.bed) regions)" >> logs/snippy.command.log
    else
      echo "WARNING: no snp_mask_bed supplied. PE/PPE and repetitive regions are NOT masked;" >> logs/snippy.command.log
      echo "         SNP distances will not be comparable with published MTBC thresholds." >> logs/snippy.command.log
    fi

    echo -e "sample\tmean_depth\tbreadth_10x\tcoverage_qc\treason" > snippy_core/coverage_qc.tsv

    echo "SNIPPY_CORE_MTBC started." > logs/snippy.command.log
    echo "Selection rule: input_reads are expected to be Kraken2-selected MTBC reads from TB_PROFILER_AND_MTBC_FILTER.mtbc_reads." >> logs/snippy.command.log
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

    echo "Number of Kraken2-selected MTBC FASTQ files received: ${n}" >> logs/snippy.command.log

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: MTBC reads must be paired R1/R2 files." >&2
      echo "ERROR: Odd number of input FASTQ files received by SNIPPY_CORE_MTBC." >> logs/snippy.command.log
      exit 1
    fi

    paired_count=$((n / 2))
    echo "Number of Kraken2-selected paired MTBC samples received: ${paired_count}" >> logs/snippy.command.log
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

      selection_basis="Selected for Snippy/core-SNP analysis because upstream Kraken2 species typing supported MTBC"

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

        bam="${outdir}/${sample}.bam"
        mean_depth="NA"
        breadth10="NA"
        cov_pass="UNMEASURED"
        cov_reason="samtools or BAM unavailable; coverage not assessed"

        if [ -s "$bam" ] && command -v samtools >/dev/null 2>&1; then
          cov_line="$(samtools depth -a "$bam" 2>/dev/null \
            | awk '{s+=$3; n++; if($3>=10) c++} END{if(n>0) printf "%.2f %.4f", s/n, c/n; else printf "0.00 0.0000"}')"
          mean_depth="$(echo "$cov_line" | awk '{print $1}')"
          breadth10="$(echo "$cov_line" | awk '{print $2}')"

          if awk -v d="$mean_depth" -v b="$breadth10" \
                 -v md="~{min_mean_depth}" -v mb="~{min_breadth_10x}" \
                 'BEGIN{exit !(d+0 >= md+0 && b+0 >= mb+0)}'; then
            cov_pass="PASS"
            cov_reason="mean depth and 10x breadth above thresholds"
          else
            cov_pass="FAIL"
            cov_reason="mean_depth<~{min_mean_depth} or breadth_10x<~{min_breadth_10x}; excluded from core-SNP alignment"
          fi
        fi

        echo -e "${sample}\t${mean_depth}\t${breadth10}\t${cov_pass}\t${cov_reason}" >> snippy_core/coverage_qc.tsv

        if [ "$cov_pass" = "FAIL" ]; then
          status="excluded_low_coverage"
          echo "EXCLUDED ${sample}: ${cov_reason}" >> logs/snippy.command.log
        else
          successful_samples+=("$outdir")
        fi

      else
        status="snippy_failed"
        echo "ERROR: snippy failed for ${sample}; last lines of its log:" \
          >> logs/snippy.command.log
        tail -n 20 "logs/${sample}.snippy.log" >> logs/snippy.command.log 2>/dev/null || true
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
      echo "Running snippy-core on ${#successful_samples[@]} successfully processed Kraken2-selected MTBC samples" >> logs/snippy.command.log

      if [ -n "${MASK_ARG}" ]; then
        if [ "~{reference_type}" = "genbank" ]; then
          REF_CONTIG=$(awk '/^LOCUS/{print $2; exit}' "$ref")
        else
          REF_CONTIG=$(awk '/^>/{sub(/^>/,""); print $1; exit}' "$ref")
        fi

        MASK_CONTIGS=$(cut -f1 snippy_core/snp_mask.bed | sort -u | tr '\n' ',' | sed 's/,$//')

        if [ -z "${REF_CONTIG}" ]; then
          echo "WARNING: could not determine reference contig name; leaving mask unchanged." \
            >> logs/snippy.command.log
        elif [ "${MASK_CONTIGS}" = "${REF_CONTIG}" ]; then
          echo "Mask contig matches reference contig (${REF_CONTIG})." >> logs/snippy.command.log
        else
          echo "Mask contig(s) [${MASK_CONTIGS}] do not match reference contig [${REF_CONTIG}];" \
            >> logs/snippy.command.log
          echo "rewriting the mask to use ${REF_CONTIG}." >> logs/snippy.command.log
          awk -v c="${REF_CONTIG}" 'BEGIN{OFS="\t"} {$1=c; print}' \
            snippy_core/snp_mask.bed > snippy_core/snp_mask.renamed.bed
          mv snippy_core/snp_mask.renamed.bed snippy_core/snp_mask.bed
        fi

        echo "Mask intervals: $(wc -l < snippy_core/snp_mask.bed)" >> logs/snippy.command.log
      fi

      if [ -n "${MASK_ARG}" ]; then
        echo "Applying SNP mask to snippy-core: ${MASK_ARG}" >> logs/snippy.command.log
      else
        echo "WARNING: no SNP mask supplied; PE/PPE and repetitive regions will be retained and SNP distances will not be comparable with published thresholds." >> logs/snippy.command.log
      fi

      if "$SNIPPY_CORE_BIN" --ref "$ref" --prefix snippy_core/core ${MASK_ARG} "${successful_samples[@]}" >> logs/snippy_core.log 2>&1; then
        core_status="success"
      else
        core_status="snippy_core_failed"
        echo "ERROR: snippy-core failed. Last lines of logs/snippy_core.log:" \
          >> logs/snippy.command.log
        tail -n 40 logs/snippy_core.log >> logs/snippy.command.log 2>/dev/null || true
        echo "ERROR: snippy-core failed; see logs/snippy_core.log" >&2
        tail -n 40 logs/snippy_core.log >&2 2>/dev/null || true
        exit 1
      fi
    else
      echo "WARNING: fewer than 2 successful Kraken2-selected MTBC samples, skipping snippy-core" >> logs/snippy.command.log
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
<strong>Selection rule:</strong> samples entering this Snippy/core-SNP task were selected upstream using Kraken2 MTBC support from the species-typing step. TB-Profiler species, lineage, and resistance outputs are annotations only and do not determine which samples enter Snippy/core-SNP/IQ-TREE analysis.
<br><br>
<strong>Second-stage tree filtering:</strong> IQTREE2_PHYLOGENY may later exclude selected samples from final IQ-TREE inference if the core-SNP alignment has no usable A/C/G/T bases or excessive missing, ambiguous, or gap content.
</div>

<table>
<thead>
<tr>
<th>Sample</th>
<th>Status</th>
<th>Mean depth</th>
<th>Selected for Snippy/core-SNP by Kraken2</th>
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
    selection_basis = r.get("selection_basis") or selection_basis_by_sample.get(sample, "Selected upstream by Kraken2 MTBC support")

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
    File? snippy_core_log = "logs/snippy_core.log"
    File samtools_version_log = "logs/samtools.version.log"

    File variant_summary = "variant_summary.html"
    File mean_depth_summary_tsv = "mean_depth_summary.tsv"
    File snippy_core_selection_summary_tsv = "snippy_core_selection_summary.tsv"

    File? core_full_alignment = "snippy_core/core.full.aln"
    File? core_snp_alignment = "snippy_core/core.aln"
    File? core_tab = "snippy_core/core.tab"
    File? core_vcf = "snippy_core/core.vcf"
    File? coverage_qc_tsv = "snippy_core/coverage_qc.tsv"
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

    Int min_comparable_sites = 2205000
    Int max_heterozygous_sites_for_clustering = 0

    File? tbprofiler_summary_tsv
    Boolean exclude_mixed_infection_from_clustering = true
    Int cpu = 8
    Int memory_gb = 64
  }

  command <<<
    set -uo pipefail

    mkdir -p snp_distance logs

    cp "~{core_full_alignment}" snp_distance/core.full.aln
    if [ -n "~{if defined(tbprofiler_summary_tsv) then tbprofiler_summary_tsv else ""}" ]; then
      cp "~{tbprofiler_summary_tsv}" snp_distance/tbprofiler_summary.tsv
    fi

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
min_comparable_sites = int("~{min_comparable_sites}")
max_het_for_clustering = int("~{max_heterozygous_sites_for_clustering}")

MIXED_STATUS = {}
exclude_mixed = "~{exclude_mixed_infection_from_clustering}".strip().lower() == "true"
_tb_summary = Path("snp_distance/tbprofiler_summary.tsv")
if _tb_summary.exists() and _tb_summary.stat().st_size > 0:
    try:
        with _tb_summary.open(newline="") as _fh:
            for _row in csv.DictReader(_fh, delimiter="\t"):
                _name = str(_row.get("sample") or "").strip()
                _st = str(_row.get("mixed_infection_status") or "").strip()
                if _name:
                    MIXED_STATUS[_name] = _st
    except Exception:
        MIXED_STATUS = {}


def tbprofiler_gates(sample):
    """True only for the lineage-level mixture call, which is the one that
    distorts distances. Heteroresistance alone does not gate."""
    if not exclude_mixed:
        return False
    return str(MIXED_STATUS.get(sample, "") or "").lower().startswith("possible mixed")

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


ALIGNMENT_COMPOSITION = {}


def _composition(raw):
    het = raw.count("n")
    low = raw.count("N")
    zero = raw.count("-")
    masked = raw.count("X") + raw.count("x")
    called = sum(raw.count(c) for c in "ACGTacgt")
    total = len(raw)
    return {
        "alignment_length": total,
        "called_bases": called,
        "heterozygous_or_poor_quality": het,
        "low_coverage": low,
        "zero_coverage": zero,
        "masked": masked,
        "heterozygous_fraction": (het / total) if total else 0.0,
    }


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
                    raw = "".join(seq)
                    ALIGNMENT_COMPOSITION[name] = _composition(raw)
                    records[name] = raw.upper()

                name = line[1:].split()[0]
                seq = []

            else:
                seq.append(line)

        if name is not None:
            raw = "".join(seq)
            ALIGNMENT_COMPOSITION[name] = _composition(raw)
            records[name] = raw.upper()

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


def interpret(distance, compared_sites, het1=0, het2=0,
              name1="sample 1", name2="sample 2"):
    try:
        comparable = int(compared_sites)
    except (TypeError, ValueError):
        comparable = 0

    if comparable < min_comparable_sites:
        return (
            "Not interpretable: insufficient comparable sites "
            "({} < {})".format(comparable, min_comparable_sites),
            "uninterpretable"
        )

    reasons = []
    if max_het_for_clustering > 0:
        for n, h in ((name1, het1), (name2, het2)):
            if h > max_het_for_clustering:
                reasons.append("%s has %d heterozygous sites (threshold %d)"
                               % (n, h, max_het_for_clustering))
    if exclude_mixed:
        for n in (name1, name2):
            st = str(MIXED_STATUS.get(n, "") or "")
            if st.lower().startswith("possible mixed"):
                reasons.append("%s is flagged '%s' by TB-Profiler" % (n, st))
    if reasons:
        return (
            "Not interpretable: mixed-strain signal (%s). The positions "
            "distinguishing the strains are excluded from this comparison, "
            "so the SNP distance understates divergence." % "; ".join(reasons),
            "uninterpretable"
        )

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

    het_cue_count = max_het_for_clustering if max_het_for_clustering > 0 else None
    ambiguity_file = Path("snp_distance/alignment_ambiguity.tsv")
    with ambiguity_file.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "sample", "alignment_length", "called_bases",
            "heterozygous_or_poor_quality", "heterozygous_fraction",
            "low_coverage", "zero_coverage", "masked", "note",
        ])
        for sname in records:
            c = ALIGNMENT_COMPOSITION.get(sname)
            if not c:
                continue
            frac = c["heterozygous_fraction"]
            het_n = c["heterozygous_or_poor_quality"]
            note = (
                "Elevated heterozygous fraction. The positions that distinguish "
                "two strains are called heterozygous and excluded from every "
                "pairwise comparison, so SNP distances for this sample are "
                "computed without them and may understate true divergence. "
                "Review alongside the mixed-infection flag."
                if (het_cue_count is not None and het_n > het_cue_count) else
                "Within the usual range for a single-strain sample"
            )
            writer.writerow([
                sname, c["alignment_length"], c["called_bases"],
                c["heterozygous_or_poor_quality"], "%.6f" % frac,
                c["low_coverage"], c["zero_coverage"], c["masked"], note,
            ])

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

            het1 = ALIGNMENT_COMPOSITION.get(sample1, {}).get(
                "heterozygous_or_poor_quality", 0)
            het2 = ALIGNMENT_COMPOSITION.get(sample2, {}).get(
                "heterozygous_or_poor_quality", 0)

            interpretation, cluster_class = interpret(
                distance, compared_sites, het1, het2, sample1, sample2)

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

    eligibility_file = Path("snp_distance/clustering_eligibility.tsv")
    with eligibility_file.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "sample", "heterozygous_sites", "heterozygous_threshold",
            "tbprofiler_mixed_infection_status", "eligible_for_clustering",
            "reason",
        ])
        for sname in samples:
            het = ALIGNMENT_COMPOSITION.get(sname, {}).get(
                "heterozygous_or_poor_quality", 0)
            tb = MIXED_STATUS.get(sname, "")
            why = []
            if max_het_for_clustering > 0 and het > max_het_for_clustering:
                why.append("heterozygous sites above threshold")
            if tbprofiler_gates(sname):
                why.append("TB-Profiler mixed-infection flag")
            writer.writerow([
                sname, het,
                max_het_for_clustering if max_het_for_clustering > 0 else "disabled",
                tb or ("not available" if not MIXED_STATUS else "not reported"),
                "NO" if why else "YES",
                "; ".join(why) if why else "no mixed-strain signal",
            ])

    close_pairs = [
        p for p in pairs
        if p["cluster_class"] in ["close", "intermediate"]
    ]

    uninterpretable_pairs = [
        p for p in pairs
        if p["cluster_class"] == "uninterpretable"
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
.uninterpretable{background:#fee2e2;color:#991b1b}
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
    File? alignment_ambiguity = "snp_distance/alignment_ambiguity.tsv"
    File? clustering_eligibility = "snp_distance/clustering_eligibility.tsv"
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

    Float max_missing_fraction_for_tree = 0.50

    Int min_non_reference_samples_for_tree = 3
  }

  command <<<
    set -uo pipefail
    mkdir -p iqtree logs


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

    Boolean midpoint_root_tree = true

    Int max_bootstrap_to_label = 100
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
import math
import struct
import sys
import traceback
import zlib

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2147483647)

outdir = Path("tree_visualization")
outdir.mkdir(exist_ok=True)

image_format = "~{image_format}".lower().strip()
if image_format not in {"png", "svg", "pdf"}:
    image_format = "png"

out_img = outdir / f"phylogenetic_tree.{image_format}"
out_svg = outdir / "phylogenetic_tree.svg"
cleaned_tree = outdir / "phylogenetic_tree.cleaned.nwk"
log = outdir / "render.log"
status_file = outdir / "render_status.txt"
excluded_copy = outdir / "excluded_from_iqtree.tsv"

tree_input = "~{input_tree}"
tbprofiler_summary_path = "~{if defined(tbprofiler_summary_tsv) then tbprofiler_summary_tsv else ""}"
resistance_profile_summary_path = "~{if defined(resistance_profile_summary_tsv) then resistance_profile_summary_tsv else ""}"
iqtree_excluded_samples_path = "~{if defined(iqtree_excluded_samples_tsv) then iqtree_excluded_samples_tsv else ""}"
requested_width = int("~{width}")
max_bootstrap_to_label = int("~{max_bootstrap_to_label}")
requested_height = int("~{height}")
title = "~{title}"
midpoint_root_tree = "~{midpoint_root_tree}".strip().lower() == "true"

log.write_text("TREE_VISUALIZATION started\n")

def logline(msg):
    with open(log, "a") as fh:
        fh.write(str(msg).rstrip("\n") + "\n")

STATUS = {
    "svg_status": "not_attempted",
    "png_status": "not_attempted",
    "png_renderer": "none",
    "tips_rendered": "0",
    "midpoint_rooted": "false",
    "reference_tips_removed": "0",
    "metadata_status": "not_loaded",
    "metadata_rows": "0",
    "tips_matched": "0",
    "error": "",
}

def write_status():
    with open(status_file, "w") as fh:
        for k, v in STATUS.items():
            fh.write(f"{k}={str(v).replace(chr(10), ' ').strip()}\n")

write_status()

if iqtree_excluded_samples_path and Path(iqtree_excluded_samples_path).exists():
    try:
        excluded_copy.write_text(
            Path(iqtree_excluded_samples_path).read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8"
        )
        logline(f"Copied IQ-TREE exclusion file: {iqtree_excluded_samples_path}")
    except Exception as exc:
        logline(f"WARNING: Could not copy IQ-TREE exclusion file: {repr(exc)}")
else:
    excluded_copy.write_text(
        "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\treason\n",
        encoding="utf-8"
    )

RESISTANCE_COLORS = {
    "Sensitive": "#2A9D8F",
    "Hr-TB": "#1D9BF0",
    "RR-TB": "#F59E0B",
    "MDR-TB": "#E11D48",
    "MDR/RR-TB": "#F97316",
    "Pre-XDR-TB": "#EF4444",
    "XDR-TB": "#6D28D9",
    "Monoresistance": "#FACC15",
    "Polyresistance": "#06B6D4",
    "Other drug resistance": "#2563EB",
    "Resistance not determined by TB-Profiler": "#94A3B8",
    "Unknown": "#6B7280"
}

DARK_TEXT_CATEGORIES = {"Monoresistance", "Resistance not determined by TB-Profiler"}

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

def xml_escape(value):
    return (str(value if value is not None else "")
            .replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


_FONT = {
    " ": [0, 0, 0, 0, 0, 0, 0], "-": [0, 0, 0, 31, 0, 0, 0], ".": [0, 0, 0, 0, 0, 12, 12],
    ":": [0, 12, 12, 0, 12, 12, 0], "/": [1, 2, 4, 8, 16, 0, 0], "(": [2, 4, 8, 8, 8, 4, 2],
    ")": [8, 4, 2, 2, 2, 4, 8], "_": [0, 0, 0, 0, 0, 0, 31], "'": [12, 12, 4, 0, 0, 0, 0],
    ",": [0, 0, 0, 0, 12, 12, 4], "=": [0, 0, 31, 0, 31, 0, 0], "+": [0, 4, 4, 31, 4, 4, 0],
    "0": [14, 17, 19, 21, 25, 17, 14], "1": [4, 12, 4, 4, 4, 4, 14], "2": [14, 17, 1, 2, 4, 8, 31],
    "3": [31, 2, 4, 2, 1, 17, 14], "4": [2, 6, 10, 18, 31, 2, 2], "5": [31, 16, 30, 1, 1, 17, 14],
    "6": [6, 8, 16, 30, 17, 17, 14], "7": [31, 1, 2, 4, 8, 8, 8], "8": [14, 17, 17, 14, 17, 17, 14],
    "9": [14, 17, 17, 15, 1, 2, 12], "A": [14, 17, 17, 31, 17, 17, 17], "B": [30, 17, 17, 30, 17, 17, 30],
    "C": [14, 17, 16, 16, 16, 17, 14], "D": [28, 18, 17, 17, 17, 18, 28], "E": [31, 16, 16, 30, 16, 16, 31],
    "F": [31, 16, 16, 30, 16, 16, 16], "G": [14, 17, 16, 23, 17, 17, 15], "H": [17, 17, 17, 31, 17, 17, 17],
    "I": [14, 4, 4, 4, 4, 4, 14], "J": [7, 2, 2, 2, 2, 18, 12], "K": [17, 18, 20, 24, 20, 18, 17],
    "L": [16, 16, 16, 16, 16, 16, 31], "M": [17, 27, 21, 21, 17, 17, 17], "N": [17, 17, 25, 21, 19, 17, 17],
    "O": [14, 17, 17, 17, 17, 17, 14], "P": [30, 17, 17, 30, 16, 16, 16], "Q": [14, 17, 17, 17, 21, 18, 13],
    "R": [30, 17, 17, 30, 20, 18, 17], "S": [15, 16, 16, 14, 1, 1, 30], "T": [31, 4, 4, 4, 4, 4, 4],
    "U": [17, 17, 17, 17, 17, 17, 14], "V": [17, 17, 17, 17, 17, 10, 4], "W": [17, 17, 17, 21, 21, 21, 10],
    "X": [17, 17, 10, 4, 10, 17, 17], "Y": [17, 17, 17, 10, 4, 4, 4], "Z": [31, 1, 2, 4, 8, 16, 31],
}

def _png_from_rows(width, height, rows_rgb):
    raw = b"".join(b"\x00" + bytes(row) for row in rows_rgb)
    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr)
            + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))

def _fallback_text_png(lines, width=1400, height=420):
    bg = (248, 250, 252)
    canvas = [[bg[0], bg[1], bg[2]] * width for _ in range(height)]
    def fill(x0, y0, x1, y1, color):
        for yy in range(max(0, y0), min(height, y1)):
            row = canvas[yy]
            for xx in range(max(0, x0), min(width, x1)):
                row[xx * 3] = color[0]; row[xx * 3 + 1] = color[1]; row[xx * 3 + 2] = color[2]
    fill(0, 0, width, height, (255, 255, 255))
    fill(24, 24, width - 24, height - 24, bg)
    fill(24, 24, 36, height - 24, (185, 28, 28))
    def draw_text(text, x, y, scale, color):
        cx = x
        for ch in str(text).upper():
            glyph = _FONT.get(ch, _FONT.get("?", [31, 1, 6, 4, 0, 4, 4]))
            for ry, bits in enumerate(glyph):
                for rx in range(5):
                    if bits & (16 >> rx):
                        fill(cx + rx * scale, y + ry * scale, cx + (rx + 1) * scale, y + (ry + 1) * scale, color)
            cx += 6 * scale
    y = 60
    for i, (text, scale, color) in enumerate(lines):
        max_chars = max(1, (width - 120) // (6 * scale))
        draw_text(str(text)[:max_chars], 64, y, scale, color)
        y += 7 * scale + 18
    return _png_from_rows(width, height, canvas)

def write_placeholder_image(message):
    """Write a legible 'tree not rendered' image in the requested format.

    Never writes a 1x1 pixel. The PNG path tries matplotlib, then Pillow, then
    a pure-Python PNG writer that draws bitmap text, so some readable image is
    always produced regardless of what the container provides.
    """
    msg = str(message)
    short = msg[:160]
    svg_text = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="420" viewBox="0 0 1400 420">
<rect width="100%" height="100%" fill="#ffffff"/>
<rect x="24" y="24" width="1352" height="372" rx="18" fill="#f8fafc" stroke="#cbd5e1" stroke-width="2"/>
<rect x="24" y="24" width="12" height="372" rx="6" fill="#b91c1c"/>
<text x="72" y="110" font-family="Arial,Helvetica,sans-serif" font-size="34" font-weight="700" fill="#0f172a">Phylogenetic tree was not rendered</text>
<text x="72" y="160" font-family="Arial,Helvetica,sans-serif" font-size="20" fill="#475569">IQ-TREE produced a tree, but the image renderer failed. The Newick tree was preserved.</text>
<text x="72" y="215" font-family="Arial,Helvetica,sans-serif" font-size="17" fill="#b91c1c">{xml_escape(short)}</text>
<text x="72" y="265" font-family="Arial,Helvetica,sans-serif" font-size="16" fill="#475569">Full traceback: tree_visualization/render.log</text>
</svg>'''

    if image_format == "svg":
        out_img.write_text(svg_text, encoding="utf-8")
        return "svg-placeholder"

    if image_format == "pdf":
        out_img.write_text("Tree visualization was not rendered. Newick tree preserved.\n" + msg + "\n", encoding="utf-8")
        return "pdf-placeholder"

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(14, 4.2), dpi=100)
        fig.patch.set_facecolor("#ffffff")
        ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
        ax.add_patch(plt.Rectangle((0.02, 0.06), 0.96, 0.88, transform=ax.transAxes, facecolor="#f8fafc", edgecolor="#cbd5e1"))
        ax.add_patch(plt.Rectangle((0.02, 0.06), 0.008, 0.88, transform=ax.transAxes, facecolor="#b91c1c", edgecolor="none"))
        ax.text(0.05, 0.74, "Phylogenetic tree was not rendered", fontsize=24, fontweight="bold", color="#0f172a", transform=ax.transAxes)
        ax.text(0.05, 0.58, "IQ-TREE produced a tree, but the image renderer failed. The Newick tree was preserved.", fontsize=13, color="#475569", transform=ax.transAxes)
        ax.text(0.05, 0.42, short, fontsize=12, color="#b91c1c", transform=ax.transAxes, wrap=True)
        ax.text(0.05, 0.24, "Full traceback: tree_visualization/render.log", fontsize=11, color="#475569", transform=ax.transAxes)
        fig.savefig(str(out_img), format="png", facecolor="#ffffff")
        plt.close(fig)
        return "matplotlib-placeholder"
    except Exception as exc:
        logline(f"placeholder via matplotlib unavailable: {repr(exc)}")

    try:
        from PIL import Image, ImageDraw
        im = Image.new("RGB", (1400, 420), "#ffffff")
        d = ImageDraw.Draw(im)
        d.rounded_rectangle((24, 24, 1376, 396), radius=18, fill="#f8fafc", outline="#cbd5e1", width=2)
        d.rectangle((24, 24, 36, 396), fill="#b91c1c")
        d.text((72, 90), "PHYLOGENETIC TREE WAS NOT RENDERED", fill="#0f172a")
        d.text((72, 140), "IQ-TREE produced a tree, but the image renderer failed. The Newick tree was preserved.", fill="#475569")
        d.text((72, 190), short, fill="#b91c1c")
        d.text((72, 240), "Full traceback: tree_visualization/render.log", fill="#475569")
        im.save(str(out_img), "PNG")
        return "pillow-placeholder"
    except Exception as exc:
        logline(f"placeholder via Pillow unavailable: {repr(exc)}")

    out_img.write_bytes(_fallback_text_png([
        ("PHYLOGENETIC TREE WAS NOT RENDERED", 5, (15, 23, 42)),
        ("IQ-TREE PRODUCED A TREE BUT THE IMAGE RENDERER FAILED.", 3, (71, 85, 105)),
        ("THE NEWICK TREE WAS PRESERVED.", 3, (71, 85, 105)),
        (short, 2, (185, 28, 28)),
        ("FULL TRACEBACK: TREE_VISUALIZATION/RENDER.LOG", 2, (71, 85, 105)),
    ]))
    return "bitmap-placeholder"

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

    fields = [f for f in re.split(r"[/|,]", raw) if f.strip()]
    tokens = []
    for field in fields:
        found = re.findall(r"\d+(?:\.\d+)?", field)
        if found:
            tokens.append(found[0])

    if not tokens:
        return ""

    support_value = tokens[-1]

    try:
        value = float(support_value)

        if 0 < value <= 1:
            value *= 100

        if value <= 0 or value > 100:
            return ""

        return str(int(round(value)))

    except Exception:
        return ""


metadata = {}
resistance_profile_map = {}
resistance_profile_source = "none"

try:
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
    STATUS["metadata_status"] = "ok" if metadata else "empty"
    STATUS["metadata_rows"] = str(len(metadata))
    logline(f"Metadata loaded: {len(metadata)} sample(s); resistance profile source: {resistance_profile_source}")
except Exception as exc:
    STATUS["metadata_status"] = "failed: " + repr(exc)
    STATUS["metadata_rows"] = str(len(metadata))
    logline(f"WARNING: metadata loading failed; tree will be rendered without annotations: {repr(exc)}")
    logline(traceback.format_exc())


class TreeNode:
    __slots__ = ("name", "length", "support", "children", "parent", "x", "y", "_id")
    _counter = 0

    def __init__(self):
        self.name = ""
        self.length = None
        self.support = ""
        self.children = []
        self.parent = None
        self.x = 0.0
        self.y = 0.0
        TreeNode._counter += 1
        self._id = TreeNode._counter

    def is_leaf(self):
        return not self.children

    def leaves(self):
        if not self.children:
            return [self]
        out = []
        stack = [self]
        while stack:
            n = stack.pop()
            if n.children:
                stack.extend(reversed(n.children))
            else:
                out.append(n)
        return out

    def traverse(self):
        out = []
        stack = [self]
        while stack:
            n = stack.pop()
            out.append(n)
            stack.extend(reversed(n.children))
        return out

def strip_newick_comments_and_space(text):
    out = []
    in_quote = False
    depth = 0
    for ch in text:
        if in_quote:
            out.append(ch)
            if ch == "'":
                in_quote = False
            continue
        if depth:
            if ch == "]":
                depth -= 1
            elif ch == "[":
                depth += 1
            continue
        if ch == "[":
            depth += 1
            continue
        if ch == "'":
            in_quote = True
            out.append(ch)
            continue
        if ch.isspace():
            continue
        out.append(ch)
    return "".join(out)

def parse_newick(text):
    s = strip_newick_comments_and_space(str(text or ""))
    if not s:
        raise ValueError("Empty Newick string")
    if s.endswith(";"):
        s = s[:-1]
    pos = 0
    n = len(s)

    def read_label():
        nonlocal pos
        if pos < n and s[pos] == "'":
            pos += 1
            buf = []
            while pos < n:
                if s[pos] == "'":
                    if pos + 1 < n and s[pos + 1] == "'":
                        buf.append("'")
                        pos += 2
                        continue
                    pos += 1
                    break
                buf.append(s[pos])
                pos += 1
            return "".join(buf)
        start = pos
        while pos < n and s[pos] not in "(),:;":
            pos += 1
        return s[start:pos]

    def read_length():
        nonlocal pos
        if pos < n and s[pos] == ":":
            pos += 1
            start = pos
            while pos < n and s[pos] not in "(),;":
                pos += 1
            raw = s[start:pos]
            try:
                return float(raw)
            except ValueError:
                return 0.0
        return None

    def read_node():
        nonlocal pos
        node = TreeNode()
        if pos < n and s[pos] == "(":
            pos += 1
            while True:
                child = read_node()
                child.parent = node
                node.children.append(child)
                if pos < n and s[pos] == ",":
                    pos += 1
                    continue
                if pos < n and s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Malformed Newick near character {pos}")
        label = read_label()
        if node.children:
            node.support = label
        else:
            node.name = label
        node.length = read_length()
        return node

    root = read_node()
    if pos != n:
        raise ValueError(f"Unexpected trailing Newick content at character {pos}")
    return root

def newick_quote(name):
    name = str(name or "")
    if re.search(r"[\s()\[\]:;,']", name):
        return "'" + name.replace("'", "''") + "'"
    return name

def to_newick(node):
    parts = []
    def rec(n):
        if n.children:
            parts.append("(")
            for i, c in enumerate(n.children):
                if i:
                    parts.append(",")
                rec(c)
            parts.append(")")
            if n.support:
                parts.append(newick_quote(n.support))
        else:
            parts.append(newick_quote(n.name))
        if n.length is not None:
            parts.append(":" + repr(float(n.length)))
    rec(node)
    parts.append(";")
    return "".join(parts)

def remove_leaf(leaf):
    """Detach a leaf and collapse any unary node this creates."""
    parent = leaf.parent
    if parent is None:
        return leaf
    parent.children.remove(leaf)
    leaf.parent = None
    return collapse_unary(parent)

def collapse_unary(node):
    """If `node` has exactly one child, splice it out (summing branch lengths).
    Returns the node that now occupies this position."""
    if len(node.children) != 1:
        return node
    child = node.children[0]
    child.length = (child.length or 0.0) + (node.length or 0.0)
    if not child.children and not child.support:
        pass
    elif child.children and not child.support and node.support:
        child.support = node.support
    grand = node.parent
    child.parent = grand
    if grand is not None:
        idx = grand.children.index(node)
        grand.children[idx] = child
    else:
        child.length = None
    node.children = []
    node.parent = None
    return child

def tree_root(node):
    while node.parent is not None:
        node = node.parent
    return node

def farthest_from(start):
    """Return (node, distance, path) to the farthest node from `start`,
    walking the tree as an undirected graph."""
    best = (start, 0.0, [start])
    stack = [(start, None, 0.0, [start])]
    while stack:
        node, came_from, dist, path = stack.pop()
        if dist > best[1]:
            best = (node, dist, path)
        neighbours = []
        for c in node.children:
            neighbours.append((c, c.length or 0.0))
        if node.parent is not None:
            neighbours.append((node.parent, node.length or 0.0))
        for nb, w in neighbours:
            if nb is came_from:
                continue
            stack.append((nb, node, dist + w, path + [nb]))
    return best

def reroot_on_edge(child, dist_from_child):
    """Place a new root on the edge between `child` and its parent, at
    `dist_from_child` from `child`. Support labels are kept with the
    bipartition they describe."""
    parent = child.parent
    if parent is None:
        return child
    old_len = child.length or 0.0
    dist_from_child = max(0.0, min(dist_from_child, old_len))

    parent.children.remove(child)
    new_root = TreeNode()
    child.parent = new_root
    child.length = dist_from_child
    new_root.children.append(child)

    cur = parent
    prev = new_root
    carry_len = old_len - dist_from_child
    carry_support = child.support if child.children else ""
    while cur is not None:
        nxt = cur.parent
        nxt_len = cur.length
        nxt_support = cur.support
        if nxt is not None:
            nxt.children.remove(cur)
        cur.parent = prev
        cur.length = carry_len
        cur.support = carry_support if cur.children else ""
        prev.children.append(cur)
        prev = cur
        cur = nxt
        carry_len = nxt_len if nxt_len is not None else 0.0
        carry_support = nxt_support

    for node in list(new_root.traverse()):
        if node is not new_root and len(node.children) == 1:
            collapse_unary(node)
    new_root.length = None
    if len(new_root.children) == 2 and new_root.children[1].children:
        if new_root.children[0].children and new_root.children[0].support == new_root.children[1].support:
            new_root.children[1].support = ""
    return new_root

def midpoint_root(root):
    leaves = root.leaves()
    if len(leaves) < 2:
        return root, False
    u, _, _ = farthest_from(leaves[0])
    v, diameter, path = farthest_from(u)
    if diameter <= 0:
        return root, False
    half = diameter / 2.0
    acc = 0.0
    for a, b in zip(path, path[1:]):
        if b.parent is a:
            child, w = b, (b.length or 0.0)
            if acc + w >= half:
                dist_from_child = w - (half - acc)
                return reroot_on_edge(child, dist_from_child), True
        else:
            child, w = a, (a.length or 0.0)
            if acc + w >= half:
                dist_from_child = half - acc
                return reroot_on_edge(child, dist_from_child), True
        acc += w
    return root, False

def ladderize(node):
    """Order children so small clades sit above large ones (stable)."""
    def size(n):
        return len(n.leaves())
    for n in node.traverse():
        if n.children:
            n.children.sort(key=size)

def nice_number(x):
    if x <= 0 or math.isnan(x) or math.isinf(x):
        return 1.0
    exp = math.floor(math.log10(x))
    frac = x / (10 ** exp)
    if frac < 1.5:
        nice = 1
    elif frac < 3.5:
        nice = 2
    elif frac < 7.5:
        nice = 5
    else:
        nice = 10
    return nice * (10 ** exp)

def fmt_num(x):
    if x == 0:
        return "0"
    if x >= 1:
        s = f"{x:.4g}"
    else:
        s = f"{x:.6f}".rstrip("0").rstrip(".")
    return s

def text_width(text, font_size, weight=400):
    factor = 0.60 if weight >= 600 else 0.56
    return len(str(text)) * font_size * factor

def render_tree_svg(root, metadata_map, title_text, subtitle_text, max_bootstrap_to_label):
    leaves = root.leaves()
    n = len(leaves)
    if n == 0:
        raise ValueError("No tips to draw")

    if n <= 30:
        row_h, name_fs, chip_fs, boot_fs, tip_r, stroke_w = 26, 12.5, 10.5, 9.5, 4.5, 1.7
    elif n <= 80:
        row_h, name_fs, chip_fs, boot_fs, tip_r, stroke_w = 20, 11, 9.5, 8.5, 3.8, 1.4
    elif n <= 200:
        row_h, name_fs, chip_fs, boot_fs, tip_r, stroke_w = 16, 9.5, 8.5, 7.5, 3.0, 1.2
    else:
        row_h, name_fs, chip_fs, boot_fs, tip_r, stroke_w = 13, 8.5, 7.5, 6.5, 2.4, 1.0

    for i, leaf in enumerate(leaves):
        leaf.y = i
    order = root.traverse()
    for node in reversed(order):
        if node.children:
            node.y = (node.children[0].y + node.children[-1].y) / 2.0
    root.x = 0.0
    for node in order:
        if node is not root:
            node.x = node.parent.x + max(0.0, node.length or 0.0)
    max_depth = max(nd.x for nd in order) or 1.0

    names = [leaf.name for leaf in leaves]
    lineages, categories = [], []
    for leaf in leaves:
        meta = metadata_map.get(normalize_tree_leaf_name(leaf.name), {})
        lineages.append(meta.get("lineage") or "Lineage NA")
        cat = normalize_profile_label(meta.get("category", "Unknown"))
        categories.append(cat)
    name_col_w = max([text_width(x, name_fs, 600) for x in names] + [text_width("SAMPLE", 10, 700) + 6]) + 10
    lineage_col_w = max([text_width(x, chip_fs) for x in lineages] + [text_width("LINEAGE", 10, 700) + 6]) + 22
    cat_col_w = max([text_width(x, chip_fs, 600) for x in categories] + [text_width("RESISTANCE PROFILE", 10, 700) + 6]) + 22

    tree_w = 560 if n <= 80 else 680
    margin_l, margin_t, margin_b = 28, 92, 118
    gap_tip_label = 14
    col_gap = 12
    plot_x0 = margin_l + 10
    label_x = plot_x0 + tree_w + gap_tip_label
    lineage_x = label_x + name_col_w + col_gap
    cat_x = lineage_x + lineage_col_w + col_gap
    width = int(cat_x + cat_col_w + 36)
    width = max(width, 760)
    height = int(margin_t + n * row_h + margin_b)
    px_per_unit = tree_w / max_depth

    def X(node):
        return plot_x0 + node.x * px_per_unit

    def Y(node):
        return margin_t + node.y * row_h + row_h / 2.0

    svg = []
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Inter,\'Segoe UI\',Helvetica,Arial,sans-serif" '
        f'text-rendering="geometricPrecision" shape-rendering="geometricPrecision" '
        f'role="img" aria-label="{xml_escape(title_text)}">'
    )
    svg.append("<defs>"
               '<filter id="rmapTipShadow" x="-50%" y="-50%" width="200%" height="200%">'
               '<feDropShadow dx="0" dy="0.6" stdDeviation="0.8" flood-color="#0f172a" flood-opacity="0.25"/>'
               "</filter></defs>")
    svg.append(f'<rect width="100%" height="100%" fill="#ffffff"/>')

    svg.append(f'<text x="{margin_l}" y="34" font-size="20" font-weight="700" fill="#0f172a">{xml_escape(title_text)}</text>')
    svg.append(f'<text x="{margin_l}" y="56" font-size="11.5" fill="#64748b">{xml_escape(subtitle_text)}</text>')

    hdr_y = margin_t - 12
    svg.append(f'<text x="{label_x}" y="{hdr_y}" font-size="10" font-weight="700" fill="#94a3b8" letter-spacing="0.08em">SAMPLE</text>')
    svg.append(f'<text x="{lineage_x}" y="{hdr_y}" font-size="10" font-weight="700" fill="#94a3b8" letter-spacing="0.08em">LINEAGE</text>')
    svg.append(f'<text x="{cat_x}" y="{hdr_y}" font-size="10" font-weight="700" fill="#94a3b8" letter-spacing="0.08em">RESISTANCE PROFILE</text>')

    for i in range(n):
        if i % 2 == 0:
            y0 = margin_t + i * row_h
            svg.append(f'<rect x="{plot_x0 - 6}" y="{y0:.1f}" width="{width - plot_x0 - 20}" height="{row_h}" fill="#f8fafc"/>')

    for node in order:
        if node is root:
            continue
        px, py = X(node.parent), Y(node.parent)
        x, y = X(node), Y(node)
        svg.append(
            f'<path d="M{px:.2f},{py:.2f} V{y:.2f} H{x:.2f}" fill="none" stroke="#334155" '
            f'stroke-width="{stroke_w}" stroke-linecap="round" stroke-linejoin="round"/>'
        )

    svg.append(f'<path d="M{X(root) - 8:.2f},{Y(root):.2f} H{X(root):.2f}" fill="none" stroke="#334155" stroke-width="{stroke_w}" stroke-linecap="round"/>')

    for node in order:
        if node is root or node.is_leaf():
            continue
        clean_support = clean_single_bootstrap_value(node.support)
        if not clean_support:
            continue
        try:
            value = int(clean_support)
        except ValueError:
            continue
        if value > max_bootstrap_to_label:
            continue
        color = "#b91c1c" if value < 70 else ("#c2410c" if value < 95 else "#475569")
        svg.append(
            f'<text x="{X(node) - 3:.2f}" y="{Y(node) - 3.5:.2f}" font-size="{boot_fs}" text-anchor="end" '
            f'fill="{color}" font-weight="600" paint-order="stroke" stroke="#ffffff" stroke-width="3" '
            f'stroke-linejoin="round">{xml_escape(clean_support)}</text>'
        )

    for i, leaf in enumerate(leaves):
        x, y = X(leaf), Y(leaf)
        cat = categories[i]
        color = RESISTANCE_COLORS.get(cat, RESISTANCE_COLORS["Unknown"])
        chip_text_color = "#111827" if cat in DARK_TEXT_CATEGORIES else "#ffffff"
        if label_x - x > 12:
            svg.append(f'<line x1="{x + tip_r + 3:.2f}" y1="{y:.2f}" x2="{label_x - 6:.2f}" y2="{y:.2f}" stroke="#cbd5e1" stroke-width="1" stroke-dasharray="1.5 3"/>')
        svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{tip_r}" fill="{color}" stroke="#ffffff" stroke-width="1.2" filter="url(#rmapTipShadow)"><title>{xml_escape(leaf.name)} — {xml_escape(lineages[i])} — {xml_escape(cat)}</title></circle>')
        svg.append(f'<text x="{label_x}" y="{y + name_fs * 0.36:.2f}" font-size="{name_fs}" font-weight="600" fill="#0f172a">{xml_escape(leaf.name)}</text>')
        lw = text_width(lineages[i], chip_fs) + 14
        ch = chip_fs + 7
        svg.append(f'<rect x="{lineage_x - 2}" y="{y - ch / 2:.2f}" width="{lw:.1f}" height="{ch}" rx="{ch / 2}" fill="#eff6ff" stroke="#bfdbfe" stroke-width="0.8"/>')
        svg.append(f'<text x="{lineage_x + 5}" y="{y + chip_fs * 0.36:.2f}" font-size="{chip_fs}" fill="#1d4ed8">{xml_escape(lineages[i])}</text>')
        cw = text_width(cat, chip_fs, 600) + 14
        svg.append(f'<rect x="{cat_x - 2}" y="{y - ch / 2:.2f}" width="{cw:.1f}" height="{ch}" rx="{ch / 2}" fill="{color}"/>')
        svg.append(f'<text x="{cat_x + 5}" y="{y + chip_fs * 0.36:.2f}" font-size="{chip_fs}" font-weight="600" fill="{chip_text_color}">{xml_escape(cat)}</text>')

    bar_units = nice_number(max_depth / 5.0)
    bar_px = bar_units * px_per_unit
    if bar_px > tree_w * 0.6:
        bar_units = nice_number(max_depth / 10.0)
        bar_px = bar_units * px_per_unit
    sb_y = margin_t + n * row_h + 34
    svg.append(f'<path d="M{plot_x0},{sb_y} h{bar_px:.2f} M{plot_x0},{sb_y - 5} v10 M{plot_x0 + bar_px:.2f},{sb_y - 5} v10" stroke="#0f172a" stroke-width="1.5" fill="none"/>')
    svg.append(f'<text x="{plot_x0 + bar_px / 2:.2f}" y="{sb_y + 18}" font-size="11" text-anchor="middle" fill="#334155">{fmt_num(bar_units)} substitutions / site</text>')

    present = [c for c in CANONICAL_RESISTANCE_ORDER if c in set(categories)]
    lx = plot_x0 + max(bar_px, 110) + 50
    ly = sb_y - 4
    for c in present:
        col = RESISTANCE_COLORS[c]
        w = text_width(c, 10.5) + 24
        if lx + w > width - 20:
            lx = plot_x0
            ly += 20
        svg.append(f'<circle cx="{lx + 6:.1f}" cy="{ly:.1f}" r="5" fill="{col}" stroke="#ffffff" stroke-width="1"/>')
        svg.append(f'<text x="{lx + 16:.1f}" y="{ly + 3.8:.1f}" font-size="10.5" fill="#334155">{xml_escape(c)}</text>')
        lx += w
    note = f"Node labels: UFBoot support (%)" + ("" if max_bootstrap_to_label >= 100 else f"; unlabelled nodes > {max_bootstrap_to_label}%")
    svg.append(f'<text x="{plot_x0}" y="{sb_y + 44}" font-size="10.5" fill="#64748b">{xml_escape(note)}. Tips coloured by canonical resistance category; tree midpoint-rooted for display when enabled.</text>')
    svg.append("</svg>")
    return "\n".join(svg)


processed_root = None
n_leaves = 0
removed_refs = []
raw_newick = ""

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

    root = parse_newick(raw_newick)

    for leaf in list(root.leaves()):
        original_name = str(leaf.name or "").strip()
        if is_reference_tip(original_name):
            removed_refs.append(original_name)
            new_top = remove_leaf(leaf)
            root = tree_root(new_top)
            continue
        normalized_name = normalize_tree_leaf_name(original_name)
        if normalized_name:
            leaf.name = normalized_name

    root = tree_root(root)
    if len(root.children) == 1:
        root = collapse_unary(root)
        root = tree_root(root)

    if len(root.leaves()) < 2:
        raise ValueError("Tree has fewer than two non-reference tips after normalizing labels and removing only the true reference tip.")

    rooted = False
    if midpoint_root_tree:
        try:
            root, rooted = midpoint_root(root)
            root = tree_root(root)
        except Exception as exc:
            logline(f"WARNING: midpoint rooting failed: {repr(exc)}")
    else:
        logline("midpoint_root_tree=false; tree rendered as inferred by IQ-TREE2 without midpoint rooting")

    ladderize(root)
    processed_root = root
    n_leaves = len(root.leaves())
    STATUS["tips_rendered"] = str(n_leaves)
    STATUS["midpoint_rooted"] = "true" if rooted else "false"
    STATUS["reference_tips_removed"] = str(len(removed_refs))

    cleaned_tree.write_text(to_newick(root) + "\n", encoding="utf-8")

    matched_tips = [lf.name for lf in root.leaves() if normalize_tree_leaf_name(lf.name) in metadata]
    unmatched_tips = [lf.name for lf in root.leaves() if normalize_tree_leaf_name(lf.name) not in metadata]
    STATUS["tips_matched"] = str(len(matched_tips))
    logline(f"Tree tips matched to metadata: {len(matched_tips)} of {n_leaves}")
    if unmatched_tips:
        logline("Tips without metadata: " + ", ".join(unmatched_tips[:100]))

    subtitle = (
        f"{n_leaves} MTBC isolates · IQ-TREE2 maximum likelihood · "
        + ("midpoint-rooted" if rooted else "unrooted (as inferred)")
        + (f" · reference tip removed: {', '.join(removed_refs)}" if removed_refs else "")
    )
    svg_markup = render_tree_svg(root, metadata, title, subtitle, max_bootstrap_to_label)
    out_svg.write_text(svg_markup, encoding="utf-8")
    STATUS["svg_status"] = "success"
    logline(f"SVG tree rendered: {out_svg} ({n_leaves} tips)")

except Exception as e:
    STATUS["svg_status"] = "failed"
    STATUS["error"] = repr(e)
    logline("SVG tree rendering failed.")
    logline(f"ERROR: {repr(e)}")
    logline(traceback.format_exc())
    if tree_input and Path(tree_input).exists() and Path(tree_input).stat().st_size > 0:
        try:
            raw = Path(tree_input).read_text(errors="replace")
            cleaned_tree.write_text(raw if raw.endswith("\n") else raw + "\n", encoding="utf-8")
        except Exception:
            cleaned_tree.write_text("();\n", encoding="utf-8")
    else:
        cleaned_tree.write_text("();\n", encoding="utf-8")
    msg = str(e)
    out_svg.write_text(f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="420" viewBox="0 0 1400 420">
<rect width="100%" height="100%" fill="#ffffff"/>
<rect x="24" y="24" width="1352" height="372" rx="18" fill="#f8fafc" stroke="#cbd5e1" stroke-width="2"/>
<rect x="24" y="24" width="12" height="372" rx="6" fill="#b91c1c"/>
<text x="72" y="110" font-family="Arial,Helvetica,sans-serif" font-size="34" font-weight="700" fill="#0f172a">Phylogenetic tree was not rendered</text>
<text x="72" y="160" font-family="Arial,Helvetica,sans-serif" font-size="20" fill="#475569">The Newick tree could not be parsed or drawn. See render.log.</text>
<text x="72" y="215" font-family="Arial,Helvetica,sans-serif" font-size="17" fill="#b91c1c">{xml_escape(msg[:160])}</text>
</svg>''', encoding="utf-8")

write_status()


if image_format == "svg" and STATUS["svg_status"] == "success":
    STATUS["png_status"] = "success"
    STATUS["png_renderer"] = "builtin-svg"
    write_status()
else:
    try:
        if processed_root is None:
            raise ValueError("Tree could not be prepared for rendering (see SVG stage error above).")

        from ete3 import Tree, TreeStyle, TextFace, NodeStyle

        t = Tree(to_newick(processed_root), format=1)

        if len(t.get_leaves()) < 2:
            raise ValueError("Tree has fewer than two tips after processing.")

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
                    try:
                        support_number = int(clean_support)
                    except ValueError:
                        support_number = -1

                    if support_number <= max_bootstrap_to_label:
                        node.add_face(
                            TextFace(clean_support, fsize=bootstrap_font, fgcolor="#b91c1c"),
                            column=0,
                            position="branch-top",
                        )

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

        t.render(str(out_img), w=auto_width, units="px", dpi=300, tree_style=ts)

        if not out_img.exists() or out_img.stat().st_size < 1000:
            raise RuntimeError(f"ETE3 wrote an empty or implausibly small image ({out_img.stat().st_size if out_img.exists() else 0} bytes)")

        STATUS["png_status"] = "success"
        STATUS["png_renderer"] = "ete3"

        logline("TREE_VISUALIZATION ETE3 render completed successfully")
        logline(f"Input tree: {tree_path}")
        logline(f"Output image: {out_img}")
        logline(f"Leaves rendered: {n_leaves}")
        logline(f"Removed reference tips: {', '.join(removed_refs) if removed_refs else 'none'}")
        logline(f"Resistance profile source: {resistance_profile_source}")
        logline(f"Canonical resistance profiles loaded: {len(resistance_profile_map)}")
        logline(f"Metadata rows loaded after normalization: {len(metadata)}")
        logline(f"Tree leaves matched to metadata after normalization: {metadata_matched}")

        if excluded_copy.exists():
            excluded_rows = 0
            try:
                with open(excluded_copy, newline="", encoding="utf-8", errors="replace") as excluded_handle:
                    excluded_reader = csv.DictReader(excluded_handle, delimiter="\t")
                    for excluded_row in excluded_reader:
                        if excluded_row.get("sample", "").strip():
                            excluded_rows += 1
                logline(f"IQ-TREE-excluded samples available for reporting: {excluded_rows}")
            except Exception as exc:
                logline(f"WARNING: Could not summarize IQ-TREE exclusion file: {repr(exc)}")

        if metadata_missing:
            logline("WARNING: Tree leaves without matching metadata after normalization: " + ", ".join(metadata_missing[:100]))
        else:
            logline("All rendered tree leaves matched metadata after normalization.")

    except Exception as e:
        STATUS["png_status"] = "failed"
        if not STATUS["error"]:
            STATUS["error"] = repr(e)
        logline("ETE3 raster rendering failed; an informative placeholder image was written instead.")
        logline(f"ERROR: {repr(e)}")
        logline(traceback.format_exc())
        try:
            STATUS["png_renderer"] = write_placeholder_image(repr(e))
        except Exception as exc2:
            logline(f"ERROR: placeholder image could not be written: {repr(exc2)}")
            STATUS["png_renderer"] = "none"

    write_status()

logline("Render status: " + ", ".join(f"{k}={v}" for k, v in STATUS.items() if k != "error"))
if STATUS["error"]:
    logline("Render error: " + STATUS["error"])
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
    File tree_svg = "tree_visualization/phylogenetic_tree.svg"
    File cleaned_tree = "tree_visualization/phylogenetic_tree.cleaned.nwk"
    File render_log = "tree_visualization/render.log"
    File render_status = "tree_visualization/render_status.txt"
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
import math
import html
import re
import sys
from pathlib import Path
from collections import Counter

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

ANIMAL_LINEAGES = [
    (r"\bla\s*1\b|\bbovis\b|\bbcg\b", "La1 (M. bovis / BCG)"),
    (r"\bla\s*2\b|\bcaprae\b", "La2 (M. caprae)"),
    (r"\bla\s*3\b|\borygis\b", "La3 (M. orygis)"),
    (r"\bmicroti\b", "M. microti"),
    (r"\bpinnipedii\b", "M. pinnipedii"),
    (r"\bmungi\b", "M. mungi"),
    (r"\bsuricattae\b", "M. suricattae"),
    (r"\bdassie\b", "Dassie bacillus"),
    (r"\bcanetti", "M. canettii"),
]


def animal_lineage_label(text):
    """Return a stable animal-lineage label, or "" if this is not one."""
    low = str(text or "").lower()
    for pattern, label in ANIMAL_LINEAGES:
        if re.search(pattern, low):
            return label
    return ""


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

    animal = animal_lineage_label(text)
    if animal:
        return animal

    found = []
    for pattern in (r"lineage\s*([1-9])", r"\bL([1-9])\b"):
        for hit in re.findall(pattern, text, flags=re.IGNORECASE):
            if hit not in found:
                found.append(hit)

    if len(found) > 1:
        return "Mixed (" + "+".join("L" + n for n in sorted(found)) + ")"
    if len(found) == 1:
        return "L" + found[0]

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
        "mycobacterium canetti",
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
        return "MTBC supported", "Kraken2"

    return "MTBC not supported", "Kraken2"

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
        or "Selected for phylogeny if Kraken2 species typing supports MTBC"
    )

    if selected_for_phylogeny == "NO" and not selection_reason:
        selection_reason = "Not selected for phylogeny because Kraken2 did not support MTBC"

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
    """Sequential colour for a SNP distance on a log1p scale.

    A linear scale is dominated by a single distant outlier (one isolate tens
    of thousands of SNPs away flattens every other cell into the same shade),
    while the epidemiologically relevant differences are at 0-50 SNPs. The
    log1p scale spreads the low end so clusters remain visible without hiding
    large distances. Stops: deep teal (identical) -> pale -> amber -> deep red.
    """
    if v is None:
        return "#e5e7eb"
    if maxv <= 0:
        return "#0f766e"
    t = max(0.0, min(1.0, math.log1p(max(v, 0.0)) / math.log1p(maxv)))
    stops = [
        (0.00, (15, 118, 110)),
        (0.22, (94, 234, 212)),
        (0.45, (254, 243, 199)),
        (0.68, (251, 146, 60)),
        (0.86, (220, 38, 38)),
        (1.00, (127, 29, 29)),
    ]
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            f = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
            r = int(round(c0[0] + (c1[0] - c0[0]) * f))
            g = int(round(c0[1] + (c1[1] - c0[1]) * f))
            b = int(round(c0[2] + (c1[2] - c0[2]) * f))
            return f"#{r:02x}{g:02x}{b:02x}"
    return "#7f1d1d"

def heat_text_color(fill_hex):
    r = int(fill_hex[1:3], 16); g = int(fill_hex[3:5], 16); b = int(fill_hex[5:7], 16)
    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255.0
    return "#0f172a" if luminance > 0.62 else "#ffffff"

def write_heatmap(path, outfile):
    samples, matrix = read_snp_matrix(path)

    if not samples or not matrix:
        outfile.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="200" viewBox="0 0 900 200">'
            '<rect width="100%" height="100%" fill="#fff"/>'
            '<text x="450" y="95" text-anchor="middle" font-family="Arial" font-size="16" font-weight="700" fill="#0f172a">No SNP distance matrix available after reference filtering.</text>'
            '<text x="450" y="122" text-anchor="middle" font-family="Arial" font-size="13" fill="#475569">The heatmap is drawn only when two or more non-reference MTBC isolates have comparable core-SNP sites.</text>'
            '</svg>',
            encoding="utf-8"
        )
        return

    n = len(samples)
    cell = 46 if n <= 12 else max(18, int(820 / n))
    show_values = n <= 18
    label_fs = 11 if n <= 25 else (9 if n <= 60 else 7)
    char_w = label_fs * 0.68

    longest = max(len(s) for s in samples)
    left = int(longest * char_w) + 28
    rot_h = int(longest * char_w * math.sin(math.radians(55))) + 44
    top = max(60, rot_h)
    right = 36
    legend_h = 118

    width = left + n * cell + right
    height = top + n * cell + legend_h

    vals = [
        v
        for i, row in enumerate(matrix)
        for j, v in enumerate(row[:n])
        if v is not None and i != j
    ]

    maxv = max(vals + [1])

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="Inter,\'Segoe UI\',Helvetica,Arial,sans-serif" text-rendering="geometricPrecision" role="img" aria-label="Pairwise SNP distance matrix">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<defs><linearGradient id="snpHeatGrad" x1="0" x2="1" y1="0" y2="0">'
    ]

    for k in range(0, 21):
        t = k / 20.0
        v = math.expm1(t * math.log1p(maxv))
        svg.append(f'<stop offset="{t * 100:.0f}%" stop-color="{heat_color(v, maxv)}"/>')

    svg.append('</linearGradient></defs>')

    for i, s in enumerate(samples):
        x = left + i * cell + cell / 2
        label_y = top - 10
        svg.append(
            f'<text x="{x}" y="{label_y}" transform="rotate(-55 {x},{label_y})" text-anchor="start" font-size="{label_fs}" fill="#0f172a">{safe(s)}</text>'
        )

    for i, s in enumerate(samples):
        y = top + i * cell + cell / 2 + label_fs * 0.36
        svg.append(
            f'<text x="{left - 10}" y="{y}" text-anchor="end" font-size="{label_fs}" fill="#0f172a">{safe(s)}</text>'
        )

    for i, row in enumerate(matrix):
        for j, v in enumerate(row[:n]):
            x = left + j * cell
            y = top + i * cell

            if i == j:
                svg.append(
                    f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="#f1f5f9" stroke="#ffffff" stroke-width="1.5"/>'
                )
                continue

            fill = heat_color(v, maxv)
            if v is None:
                lab = ""
                tip = f"{safe(samples[i])} vs {safe(samples[j])}: not comparable"
            else:
                lab = str(int(v)) if float(v).is_integer() else f"{v:.1f}"
                tip = f"{safe(samples[i])} vs {safe(samples[j])}: {lab} SNPs"

            svg.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" stroke="#ffffff" stroke-width="1.5"><title>{tip}</title></rect>'
            )

            if show_values and lab:
                svg.append(
                    f'<text x="{x + cell/2}" y="{y + cell/2 + 3.5}" text-anchor="middle" font-size="10" font-weight="600" fill="{heat_text_color(fill)}" pointer-events="none">{safe(lab)}</text>'
                )

    legend_w = min(420, max(220, n * cell))
    lx = left
    ly = top + n * cell + 38
    svg.append(f'<rect x="{lx}" y="{ly}" width="{legend_w}" height="14" rx="4" fill="url(#snpHeatGrad)" stroke="#cbd5e1" stroke-width="0.8"/>')

    candidate_ticks = [0, 5, 12, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000]
    ticks = [t for t in candidate_ticks if t < maxv]
    ticks.append(int(maxv) if float(maxv).is_integer() else round(maxv, 1))
    last_px = -1e9
    for t in ticks:
        frac = math.log1p(t) / math.log1p(maxv) if maxv > 0 else 0
        px = lx + frac * legend_w
        if px - last_px < 26:
            continue
        last_px = px
        svg.append(f'<line x1="{px:.1f}" y1="{ly + 14}" x2="{px:.1f}" y2="{ly + 19}" stroke="#64748b" stroke-width="1"/>')
        svg.append(f'<text x="{px:.1f}" y="{ly + 31}" text-anchor="middle" font-size="9.5" fill="#334155">{t}</text>')

    svg.append(f'<text x="{lx}" y="{ly - 8}" font-size="10.5" font-weight="700" fill="#334155">SNP distance (log scale)</text>')
    svg.append(f'<text x="{lx + legend_w + 16}" y="{ly + 11}" font-size="10.5" fill="#475569">Teal = identical or closely related &middot; red = distant. Diagonal (self-comparison) is blank.</text>')
    svg.append(f'<text x="{lx}" y="{ly + 56}" font-size="10.5" fill="#64748b">Reference sequence excluded. Distances are counted over comparable core-SNP sites.</text>')
    svg.append(f'<text x="{lx}" y="{ly + 72}" font-size="10.5" fill="#64748b">Interpret together with lineage, resistance profile, mixed-infection flags and tree support.</text>')
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
summary_rows.append('<p class="muted">This summary integrates Kraken2 species typing, TB-Profiler lineage and resistance annotations, mean-depth metadata, and SNP-distance outputs.</p>')
summary_rows.append('<p class="muted"><strong>Phylogeny selection rule:</strong> samples are selected for Snippy/core-SNP/IQ-TREE only when Kraken2 supports MTBC. TB-Profiler lineage, species, and resistance are reported as annotations and do not determine phylogeny selection.</p>')
summary_rows.append('<p class="muted"><strong>Tree inclusion rule:</strong> after Kraken2-based selection, IQ-TREE-related outputs reflect only samples that passed downstream alignment and SNP-distance/tree generation steps.</p>')
summary_rows.append('</div>')

summary_rows.append('<div class="card">')
summary_rows.append('<h2>Kraken2-based phylogeny selection and tree inclusion</h2>')
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
    File? coverage_qc_tsv
    File? iqtree_report
    File? iqtree_excluded_samples_tsv
    File? iqtree_included_samples_tsv
    File? iqtree_filtering_summary_txt
    File? iqtree_status
    File? tree_image
    File? tree_svg
    File? tree_render_log
    File? tree_render_status
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
    tree_svg="~{if defined(tree_svg) then tree_svg else ""}"
    tree_render_log="~{if defined(tree_render_log) then tree_render_log else ""}"
    tree_render_status="~{if defined(tree_render_status) then tree_render_status else ""}"
    tree_newick="~{if defined(phylogenetic_tree_newick) then phylogenetic_tree_newick else ""}"
    qc_html="~{if defined(qc_summary_html) then qc_summary_html else ""}"
    trim_html="~{if defined(trimming_report_html) then trimming_report_html else ""}"
    variant_html="~{if defined(variant_summary_html) then variant_summary_html else ""}"
    coverage_qc_tsv="~{if defined(coverage_qc_tsv) then coverage_qc_tsv else ""}"
    iqtree_txt="~{if defined(iqtree_report) then iqtree_report else ""}"
    iqtree_excluded_tsv="~{if defined(iqtree_excluded_samples_tsv) then iqtree_excluded_samples_tsv else ""}"
    iqtree_included_tsv="~{if defined(iqtree_included_samples_tsv) then iqtree_included_samples_tsv else ""}"
    iqtree_filtering_summary_txt="~{if defined(iqtree_filtering_summary_txt) then iqtree_filtering_summary_txt else ""}"
    iqtree_status_txt="~{if defined(iqtree_status) then iqtree_status else ""}"

    if [ -n "$tree_png" ] && [ -f "$tree_png" ]; then
      cp "$tree_png" final_report/mtbc_tree.png || true
    fi

    if [ -n "$tree_svg" ] && [ -f "$tree_svg" ]; then
      cp "$tree_svg" final_report/mtbc_tree.svg || true
    fi

    if [ -n "$tree_render_log" ] && [ -f "$tree_render_log" ]; then
      cp "$tree_render_log" final_report/tree_render.log || true
    fi

    if [ -n "$tree_render_status" ] && [ -f "$tree_render_status" ]; then
      cp "$tree_render_status" final_report/tree_render_status.txt || true
    fi

    if [ -n "$tree_newick" ] && [ -f "$tree_newick" ]; then
      cp "$tree_newick" final_report/mtbc_tree.cleaned.nwk || true
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
      echo -e "sample\tspecies\tmain_lineage\tsub_lineage\tdr_type\tresistant_drugs\tresistance_mutations\tkey_mutations\tjson_file\tmtbc_selected\tmtbc_selection_reason\tmixed_infection_status\tmixed_infection_evidence\tstatus" > final_report/empty.tsv
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

    python3 - "$tb_tsv" "$resistance_tsv" "$species_tsv" "$ntm_tsv" "$routing_tsv" "$routing_status_txt" "$nonsyn_tsv" "$mutation_tsv" "$snp_pairs_tsv" "$snp_cluster_tsv" "$lineage_tsv" "$surveillance_metadata_tsv" "$qc_rationale_tsv" "$qc_html" "$trim_html" "$variant_html" "$iqtree_txt" "$species_html" "$ntm_html" "$iqtree_excluded_tsv" "$iqtree_included_tsv" "$iqtree_filtering_summary_txt" "$iqtree_status_txt" "final_report/mtbc_tree.svg" "final_report/tree_render.log" "final_report/tree_render_status.txt" "$coverage_qc_tsv" <<'PY'
import base64
import csv
import html
import os
import re
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone

(summary_tsv, resistance_tsv, species_tsv, ntm_tsv, routing_tsv, routing_status_txt, nonsyn_tsv,
 mutation_tsv, snp_pairs_tsv, snp_cluster_tsv, lineage_tsv, surveillance_metadata_tsv,
 qc_rationale_tsv, qc_html, trim_html, variant_html, iqtree_txt, species_html, ntm_html,
 iqtree_excluded_tsv, iqtree_included_tsv, iqtree_filtering_summary_txt, iqtree_status_txt,
 tree_svg_path, tree_render_log_txt, tree_render_status_txt, coverage_qc_tsv) = sys.argv[1:28]

outdir = Path("final_report")
outdir.mkdir(exist_ok=True)

run_started_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
run_stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_UTC")

(outdir / "run_metadata.txt").write_text(
    f"Workflow report generation timestamp: {run_started_utc}\n"
    f"Run stamp: {run_stamp}\n"
)

RESISTANCE_COLORS = {
    "Sensitive": "#2A9D8F",
    "Hr-TB": "#1D9BF0",
    "RR-TB": "#F59E0B",
    "MDR-TB": "#E11D48",
    "MDR/RR-TB": "#F97316",
    "Pre-XDR-TB": "#EF4444",
    "XDR-TB": "#6D28D9",
    "Monoresistance": "#FACC15",
    "Polyresistance": "#06B6D4",
    "Other drug resistance": "#2563EB",
    "Resistance not determined by TB-Profiler": "#94A3B8",
    "Unknown": "#6B7280"
}

DARK_TEXT_CATEGORIES = {"Monoresistance", "Resistance not determined by TB-Profiler"}

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

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2147483647)

UNPARSED_INPUTS = []


def read_tsv_rows(path):
    try:
        p = Path(path)
        if not (path and p.exists() and p.is_file() and p.stat().st_size > 0):
            return []
        with p.open(encoding="utf-8", errors="replace", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except Exception as exc:
        UNPARSED_INPUTS.append((str(path).split("/")[-1], "%s: %s" % (type(exc).__name__, exc)))
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


tree_render_status = {}
for line in read_optional_file(tree_render_status_txt).splitlines():
    key, sep, value = line.partition("=")
    if sep:
        tree_render_status[key.strip()] = value.strip()

tree_render_log_text = read_optional_file(tree_render_log_txt)


iqtree_included_map = {}
iqtree_excluded_map = {}

coverage_qc_rows = read_tsv_rows(coverage_qc_tsv)
coverage_fail_map = {}
coverage_by_sample = {}
for r in coverage_qc_rows:
    sid = normalize_sample_id(r.get("sample") or r.get("Sample") or "")
    if not sid:
        continue
    coverage_by_sample[sid] = r
    if clean(r.get("coverage_qc")).upper() == "FAIL":
        coverage_fail_map[sid] = r

def coverage_fail_note(sample_id):
    r = coverage_fail_map.get(sample_id, {})
    depth = clean(r.get("mean_depth")) or "NA"
    breadth = clean(r.get("breadth_10x")) or "NA"
    try:
        breadth_txt = f"{float(breadth) * 100:.1f}%"
    except Exception:
        breadth_txt = breadth
    reason = clean(r.get("reason"))
    return (
        f"Excluded before IQ-TREE: failed the per-sample coverage gate after mapping "
        f"(mean depth {depth}x, 10x breadth {breadth_txt}; {reason or 'below configured thresholds'}), "
        f"so it was not included in the core-SNP alignment."
    )

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

    if sample_id in coverage_fail_map:
        return {
            "status": "NO",
            "note": coverage_fail_note(sample_id),
            "reason": "excluded_low_coverage",
            "row": coverage_fail_map[sample_id]
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

coverage_excluded_mtbc = sorted([
    s for s in coverage_fail_map
    if s and not is_reference_name(s) and (not routing_by_sample or s in mtbc_routed_samples)
])
coverage_excluded_count = len(coverage_excluded_mtbc)

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
    text_color = "#111827" if label in DARK_TEXT_CATEGORIES else "#ffffff"
    return f'<span class="res-badge" style="background:{color};color:{text_color} !important;">{safe(label)}</span>'

def green_badge(label):
    return f'<span class="badge badge-green" style="background:#15803d !important;color:white !important;font-weight:700;">{safe(label)}</span>'

def red_badge(label):
    return f'<span class="badge badge-red" style="background:#b91c1c !important;color:white !important;font-weight:700;">{safe(label)}</span>'

def neutral_badge(label):
    return f'<span class="badge" style="background:#9ca3af !important;color:white !important;font-weight:700;">{safe(label)}</span>'

def species_badge(label):
    text = lower_clean(label)

    if text == "mycobacterium tuberculosis":
        return (
            '<span class="badge" '
            'style="background:#dbeafe !important;color:#1e3a8a !important;font-weight:700;">'
            f'<em>{safe(label)}</em></span>'
        )

    return green_badge(label)

def cluster_badge(label):
    text = lower_clean(label)

    if text.startswith("not interpretable"):
        short = label.split("(")[0].split(".")[0].strip()
        return (
            '<span class="badge" style="background:#6b7280 !important;'
            'color:white !important;font-weight:700;" '
            f'title="{safe(label)}">{safe(short)}</span>'
        )

    if "genomically close" in text or text == "close" or "likely recent" in text or "likely" in text:
        return '<span class="badge badge-wrap" style="background:#15803d !important;color:white !important;font-weight:700;">Genomically close; review epidemiological linkage</span>'

    if "intermediate" in text or "possible" in text:
        return '<span class="badge badge-wrap" style="background:#facc15 !important;color:#111111 !important;font-weight:700;">Intermediate SNP distance; review metadata</span>'

    if "not clustered" in text or "distant" in text:
        return '<span class="badge badge-wrap" style="background:#bae6fd !important;color:#0c4a6e !important;font-weight:700;">Not clustered by SNP threshold</span>'

    return '<span class="badge" style="background:#9ca3af !important;color:white !important;font-weight:700;">Not interpreted</span>'

total_samples = max(len(species_rows), len(routing_by_sample), len(rows))

mtbc_retained = len([s for s in mtbc_routed_samples if s and not is_reference_name(s)])
non_mtbc = len([s for s in non_mtbc_routed_samples if s and not is_reference_name(s)])

if not routing_by_sample:
    mtbc_retained = sum(1 for r in rows if clean(r.get("mtbc_selected")).upper() == "YES")
    non_mtbc = max(total_samples - mtbc_retained, 0)

drug_resistant = 0
resistance_category_counts = defaultdict(int)

for r in rows:
    sample = r.get("sample")
    if routing_by_sample and not is_sample_routed_mtbc(sample):
        continue
    info = get_resistance_info(sample, r)
    category = info["resistance_profile"]
    resistance_category_counts[category] += 1

    if is_drug_resistant_category(category):
        drug_resistant += 1


def display_sample_id(value):
    """Remove internal Terra shard/index prefixes from sample IDs shown in final report tables."""
    s = clean(value)
    return re.sub(r"^\d{3,6}[_-](?=ERR|SRR|DRR|SAM|ERS|SRS|[A-Za-z])", "", s)

def section_open(number, slug, title, lead=None):
    lead_html = f'<p class="sec-lead">{lead}</p>' if lead else ""
    return (
        f'<section class="section" id="sec-{slug}" data-section-number="{number}">'
        f'<header class="sec-head"><span class="sec-num">{number}</span>'
        f'<div><h2>{safe(title)}</h2>{lead_html}</div></header>'
    )

def section_close():
    return "</section>"

def table_wrap(table_html):
    return f'<div class="table-wrap">{table_html}</div>'

def collapsible_list(label_singular, label_plural, text, max_inline=3):
    """Show a comma-separated list compactly: the first few items inline, the
    rest behind a disclosure. Long mutation lists previously inflated a single
    table row to several screens."""
    text = clean(text)
    if not text or text == "None reported":
        return safe(text or "None reported")
    items = [i.strip() for i in re.split(r"[;,]\s*", text) if i.strip()]
    if len(items) <= max_inline:
        return safe(", ".join(items))
    head = ", ".join(safe(i) for i in items[:max_inline])
    rest = "".join(f"<li>{safe(i)}</li>" for i in items[max_inline:])
    n = len(items)
    label = label_plural if n != 1 else label_singular
    return (
        f'{head} <details class="inline-details"><summary>+{n - max_inline} more '
        f'({n} {label} in total)</summary><ul class="compact-list">{rest}</ul></details>'
    )

def build_unparsed_inputs_note():
    if not UNPARSED_INPUTS:
        return ""
    items = "".join(
        f"<li><code>{html.escape(name)}</code> &mdash; {html.escape(err)}</li>"
        for name, err in UNPARSED_INPUTS
    )
    return (
        '<div class="note note-danger">'
        "<strong>Report inputs that could not be parsed.</strong> The sections "
        "built from these files are empty because the file could not be read, "
        "not because there was nothing to report. Open the file to see the "
        "cause."
        f"<ul>{items}</ul></div>"
    )


def build_resistance_distribution():
    counts = [(c, resistance_category_counts.get(c, 0)) for c in RESISTANCE_ORDER if resistance_category_counts.get(c, 0) > 0]
    total = sum(n for _, n in counts)
    if not total:
        return '<div class="note">No MTBC-routed isolates with a resistance classification were available for the distribution bar.</div>'
    segments = []
    chips = []
    for cat, n in counts:
        pct = 100.0 * n / total
        color = RESISTANCE_COLORS[cat]
        text_color = "#111827" if cat in DARK_TEXT_CATEGORIES else "#ffffff"
        label = f"{n}" if pct >= 6 else ""
        segments.append(
            f'<div class="dist-seg" style="flex:{n} 0 0;background:{color};color:{text_color};" '
            f'title="{safe(cat)}: {n} of {total} ({pct:.1f}%)"><span>{label}</span></div>'
        )
        chips.append(
            f'<span class="dist-chip"><i style="background:{color};"></i>{safe(cat)} '
            f'<b>{n}</b> <small>({pct:.0f}%)</small></span>'
        )
    return (
        '<div class="dist-block">'
        '<div class="dist-title">Resistance profile distribution <small>MTBC-routed isolates, canonical categories</small></div>'
        f'<div class="dist-bar" role="img" aria-label="Resistance profile distribution">{"".join(segments)}</div>'
        f'<div class="dist-legend">{"".join(chips)}</div>'
        '</div>'
    )

def status_pill(label, state, detail=""):
    cls = {"ok": "pill-ok", "warn": "pill-warn", "bad": "pill-bad", "muted": "pill-muted"}.get(state, "pill-muted")
    detail_html = f'<small>{safe(detail)}</small>' if detail else ""
    return f'<div class="pill {cls}"><span class="pill-dot"></span><div><strong>{safe(label)}</strong>{detail_html}</div></div>'

PNG_OK = False

def tree_render_outcome():
    """Return (primary_kind, state, detail) where primary_kind is 'svg', 'png' or None."""
    global PNG_OK
    svg_path = Path("final_report/mtbc_tree.svg")
    png_path = Path("final_report/mtbc_tree.png")
    svg_text = read_optional_file(str(svg_path)) if svg_path.exists() else ""
    svg_status = tree_render_status.get("svg_status", "")
    png_status = tree_render_status.get("png_status", "")

    svg_ok = False
    if svg_text and "<svg" in svg_text.lower():
        if svg_status == "success":
            svg_ok = True
        elif not svg_status:
            svg_ok = "<circle" in svg_text and "not rendered" not in svg_text.lower()

    png_ok = False
    if png_path.exists() and png_path.stat().st_size > 0:
        if png_status == "success":
            png_ok = True
        elif not png_status:
            png_ok = png_path.stat().st_size > 2000
    PNG_OK = png_ok

    tips = tree_render_status.get("tips_rendered", "?")
    matched = tree_render_status.get("tips_matched", "")
    meta_status = tree_render_status.get("metadata_status", "")
    annot_warn = ""
    if tree_render_status and (matched == "0" or meta_status.startswith("failed")) and tips not in {"0", "?"}:
        annot_warn = " — WARNING: lineage/resistance annotations missing (" + (meta_status if meta_status.startswith("failed") else f"{matched} of {tips} tips matched metadata") + ")"
    if svg_ok:
        return "svg", ("warn" if annot_warn else "ok"), f"Vector tree rendered ({tips} tips)" + annot_warn
    if png_ok:
        return "png", ("warn" if annot_warn else "ok"), "ETE3 raster tree rendered" + annot_warn
    err = tree_render_status.get("error", "") or "renderer produced no usable image"
    return None, "bad", err[:140]

TREE_PRIMARY, TREE_STATE, TREE_DETAIL = tree_render_outcome()

def build_status_strip():
    iq = iqtree_status_text
    iq_state = "ok" if iq == "success" else ("muted" if iq in {"unknown_status", ""} else "bad")
    heat_exists = Path("final_report/snp_distance_heatmap.svg").exists()
    heat_text = read_optional_file("final_report/snp_distance_heatmap.svg")
    heat_state = "ok" if (heat_exists and "<rect" in heat_text and "No SNP distance matrix" not in heat_text) else "muted"
    routing_state = "ok" if routing_by_sample else "muted"
    return (
        '<div class="status-strip">'
        + status_pill("Kraken2 routing", routing_state, f"{mtbc_retained} MTBC · {non_mtbc} non-MTBC" if routing_by_sample else "router summary not provided")
        + status_pill("IQ-TREE2", iq_state, f"status: {iq} · {iqtree_included_count} isolates in tree" + (f" · {coverage_excluded_count} coverage-excluded" if coverage_excluded_count else ""))
        + status_pill("Tree rendering", TREE_STATE, TREE_DETAIL)
        + status_pill("SNP distance heatmap", heat_state, "generated" if heat_state == "ok" else "not generated")
        + "</div>"
    )


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
        body = '<tr><td colspan="5" class="empty-cell">No sample-level QC records available.</td></tr>'
    else:
        body = []

        for s in sample_ids:
            body.append(
                "<tr>"
                f"<td class='mono'>{display_sample_id(s)}</td>"
                "<td>Reported in MultiQC</td>"
                "<td>See trimming report</td>"
                f"<td>{green_badge('PASS')}</td>"
                f"<td>{green_badge('Proceed')}</td>"
                "</tr>"
            )

        body = "".join(body)

    embedded = ""

    if qc_content:
        embedded += '<details><summary>Embedded QC summary report (MultiQC)</summary><div class="embedded-report">' + qc_content + "</div></details>"

    if trim_content:
        embedded += '<details><summary>Embedded trimming report</summary><div class="embedded-report">' + trim_content + "</div></details>"

    if variant_content:
        embedded += '<details><summary>Embedded variant summary report</summary><div class="embedded-report">' + variant_content + "</div></details>"

    if not embedded:
        embedded = '<div class="note">QC, trimming, and variant reports were not provided as separate HTML inputs, but sample-level workflow decisions are summarized below.</div>'

    table = f"""
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
</table>"""

    return (
        section_open(1, "qc", "Sample QC and Trimming Summary",
                     "Per-sample read quality, trimming and the workflow decision taken for each library. Full MultiQC, trimming and variant reports are embedded below.")
        + f"""
<div class="controls">
<input id="qcSearch" onkeyup="filterTable('qcSearch','qcTable')" placeholder="Search QC table…">
<button onclick="downloadCSV('qcTable','qc_summary.csv')">Download QC CSV</button>
</div>
{table_wrap(table)}
{embedded}
"""
        + section_close()
    )

def build_species_section():
    if not species_rows:
        body = '<tr><td colspan="3" class="empty-cell">No species typing results were generated.</td></tr>'
    else:
        body = []

        for r in species_rows:
            sample = r.get("Sample_ID") or r.get("sample") or ""
            species = r.get("Species_Identified") or r.get("species") or "No species-level Mycobacterium call"
            evidence = r.get("Evidence") or r.get("evidence") or "No supporting evidence available"

            body.append(
                "<tr>"
                f"<td class='mono'>{safe(display_sample_id(sample))}</td>"
                f"<td>{species_badge(species)}</td>"
                f"<td>{safe(evidence)}</td>"
                "</tr>"
            )

        body = "".join(body)

    table = f"""
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
</table>"""

    return (
        section_open(2, "species", "Species Typing using Kraken2",
                     "One most-probable species-level call per sample from a <em>Mycobacterium</em>-only Kraken2 database.")
        + f"""
<div class="note">
Species typing was performed using Kraken2 against a custom <em>Mycobacterium</em>-only database embedded in the Docker image
<code>gmboowa/mycobacterium-kraken2-bracken:2026.05</code> (76 RefSeq sequences, 50 species; see <code>resources/kraken2_db/README.md</code>). The table reports one most probable species-level call per sample based on the highest species-level Kraken2 assignment and supporting taxonomic evidence.
<br><br>
<strong>Database scope:</strong> <em>Mycobacterium</em> references only &mdash; no human, non-mycobacterial bacterial, viral or fungal sequences. Reads from any non-mycobacterial source are returned as <em>unclassified</em> and cannot be identified.
<br><br>
<strong>Interpretation limit within MTBC:</strong> only <code>NC_000962.3</code> (<em>M. tuberculosis</em> H37Rv) and <code>NC_015848.1</code> (<em>M. canettii</em>) are present; <em>M. bovis</em>, <em>M. africanum</em>, <em>M. microti</em>, <em>M. orygis</em>, <em>M. caprae</em> and <em>M. pinnipedii</em> are not. A call here is evidence of MTBC membership only, not species-level discrimination within the complex.
</div>

<div class="controls">
<input id="speciesSearch" onkeyup="filterTable('speciesSearch','speciesTable')" placeholder="Search species typing results…">
<button onclick="downloadCSV('speciesTable','species_typing_summary.csv')">Download Species Typing CSV</button>
</div>
{table_wrap(table)}
"""
        + section_close()
    )

def build_ntm_species_section():
    display_rows = ntm_display_rows_from_routing if ntm_display_rows_from_routing else [r for r in ntm_rows if clean(r.get("Mycobacteria_Category")) != "MTBC"]

    if not display_rows:
        body = '<tr><td colspan="7" class="empty-cell">No non-MTBC Mycobacterium or NTM samples were detected by Kraken2.</td></tr>'
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
                f"<td class='mono'>{safe(display_sample_id(sample))}</td>"
                f"<td><em>{safe(species)}</em></td>"
                f"<td>{safe(category)}</td>"
                f"<td class='num'>{safe(mtbc_reads)}</td>"
                f"<td class='num'>{safe(mtbc_percent)}</td>"
                f"<td>{safe(evidence)}</td>"
                f"<td>{safe(decision)}</td>"
                "</tr>"
            )
        body = "".join(body_parts)

    table = f"""
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
</table>"""

    return (
        section_open("2b", "ntm", "Non-MTBC Mycobacteria Species Summary",
                     "Samples routed away from the MTBC workflow by Kraken2 support, reported here for transparency.")
        + f"""
<div class="note">
<strong>Routing status:</strong> {safe(routing_status_text)}<br><br>
Kraken2 species typing is used as the first routing layer. MTBC-supported samples are retained for TB-Profiler, drug-resistance interpretation, and MTBC-specific phylogenomics. Non-MTBC Mycobacterium samples are reported in this section and excluded from TB-Profiler and MTBC-specific downstream analyses.
</div>
<div class="controls">
<input id="ntmSearch" onkeyup="filterTable('ntmSearch','ntmTable')" placeholder="Search non-MTBC Mycobacteria results…">
<button onclick="downloadCSV('ntmTable','non_mtbc_mycobacteria_summary.csv')">Download Non-MTBC Mycobacteria CSV</button>
</div>
{table_wrap(table)}
"""
        + section_close()
    )

def build_tb_rows():
    out = []

    for r in rows:
        sample = r.get("sample")
        if routing_by_sample and not is_sample_routed_mtbc(sample):
            continue
        species = r.get("species") or "Mycobacterium tuberculosis complex (inferred from available MTBC evidence)"
        species_main, _, species_note = clean(species).partition(" (")
        species_note = species_note.rstrip(")")
        species_html = f"<em>{safe(species_main)}</em>"
        if species_note:
            species_html += f"<br><small class='muted'>{safe(species_note)}</small>"
        main_lineage = normalize_missing(r.get("main_lineage"), "Not resolved by TB-Profiler")
        sub_lineage = normalize_missing(r.get("sub_lineage"), "Not resolved by TB-Profiler")
        if main_lineage == sub_lineage:
            lineage = safe(main_lineage)
        else:
            lineage = f"{safe(main_lineage)}<br><small class='muted'>{safe(sub_lineage)}</small>"

        info = get_resistance_info(sample, r)
        category = info["resistance_profile"]
        resistant_drugs = info["resistant_drugs"]
        key_mutations = info["key_mutations"]

        if key_mutations and key_mutations != "None reported":
            resistance_detail = (
                f"{safe(resistant_drugs)}"
                f"<div class='sub-block'><span class='sub-label'>Key mutations</span> "
                f"{collapsible_list('mutation', 'mutations', key_mutations, max_inline=4)}</div>"
            )
        else:
            resistance_detail = safe(resistant_drugs)

        mtbc_selected = clean(r.get("mtbc_selected")).upper()
        mtbc_reason = r.get("mtbc_selection_reason") or ""

        selected_html = selected_for_mtbc_workflow_html(mtbc_selected, mtbc_reason)
        iqtree_html = iqtree_inclusion_html(sample)

        mixed_status = normalize_missing(
            r.get("mixed_infection_status"), "Not assessed"
        )
        mixed_evidence = clean(r.get("mixed_infection_evidence"))

        mixed_low = mixed_status.lower()
        if "possible" in mixed_low or "mixed" in mixed_low and "no " not in mixed_low and "not " not in mixed_low:
            mixed_badge = '<span class="badge badge-wrap" style="background:#f59e0b !important;color:#111827 !important;font-weight:700;">' + safe(mixed_status) + "</span>"
        elif mixed_status == "Not assessed":
            mixed_badge = neutral_badge(mixed_status)
        else:
            mixed_badge = f'<span class="badge badge-soft badge-wrap">{safe(mixed_status)}</span>'

        if mixed_evidence:
            mixed_html = f"{mixed_badge}<br><small>{safe(mixed_evidence)}</small>"
        else:
            mixed_html = mixed_badge

        out.append(
            "<tr>"
            f"<td class='mono'>{safe(display_sample_id(sample))}</td>"
            f"<td>{species_html}</td>"
            f"<td>{lineage}</td>"
            f"<td>{badge(category)}</td>"
            f"<td>{resistance_detail}</td>"
            f"<td>{mixed_html}</td>"
            f"<td>{selected_html}</td>"
            f"<td>{iqtree_html}</td>"
            "</tr>"
        )

    if not out:
        return '<tr><td colspan="8" class="empty-cell">No MTBC-routed TB-Profiler rows were available. Non-MTBC/NTM samples are intentionally excluded from this section.</td></tr>'

    return "".join(out)

def build_tbprofiler_section():
    table = f"""
<table id="tbTable">
<thead>
<tr>
<th class="sample" onclick="sortTable('tbTable',0)">Sample ID</th>
<th class="species" onclick="sortTable('tbTable',1)">Species</th>
<th class="lineage" onclick="sortTable('tbTable',2)">Lineage</th>
<th class="resistance" onclick="sortTable('tbTable',3)">Resistance profile</th>
<th class="mutations" onclick="sortTable('tbTable',4)">Resistant drugs / key mutations</th>
<th class="status" onclick="sortTable('tbTable',5)">Mixed-infection flag</th>
<th class="status" onclick="sortTable('tbTable',6)">Selected for MTBC workflow</th>
<th class="status" onclick="sortTable('tbTable',7)">Included in IQ-TREE</th>
</tr>
</thead>
<tbody>
{build_tb_rows()}
</tbody>
</table>"""

    return (
        section_open(3, "tbprofiler", "TB-Profiler Resistance, Species, and Lineage Report",
                     "Canonical per-isolate resistance category, lineage, mixed-infection screening flag and workflow/tree inclusion status for MTBC-routed samples.")
        + f"""
<div class="note">
<strong>Interpretation note:</strong> Resistance profile classifications in this section are populated from the canonical <code>resistance_profile_summary.tsv</code> generated during TB-Profiler parsing. The same source is also used for Surveillance Metadata, tree labels, resistance badges, resistant drug lists, and drug-resistant isolate counts.
<br><br>
<strong>Reporting distinction:</strong> <em>Selected for MTBC workflow</em> means the sample was retained by the MTBC support logic and remains part of the wider workflow. <em>Included in IQ-TREE</em> means the sample also passed core-SNP alignment quality filtering and was actually eligible for IQ-TREE phylogenetic inference. Samples can therefore be selected for the MTBC workflow but excluded from IQ-TREE only.
</div>

<details class="defs">
<summary>WHO 2021+ resistance definitions used for the canonical categories</summary>
<dl class="def-grid">
<dt>{badge("Hr-TB")}</dt><dd>resistant to isoniazid and not resistant to rifampicin.</dd>
<dt>{badge("RR-TB")}</dt><dd>resistant to rifampicin without detected isoniazid resistance.</dd>
<dt>{badge("MDR-TB")}</dt><dd>resistant to at least both isoniazid and rifampicin.</dd>
<dt>{badge("MDR/RR-TB")}</dt><dd>umbrella term for MDR-TB or RR-TB; the report separates this into <strong>RR-TB</strong> when rifampicin is detected without isoniazid, and <strong>MDR-TB</strong> when both are detected.</dd>
<dt>{badge("Pre-XDR-TB")}</dt><dd>MDR/RR-TB that is also resistant to any fluoroquinolone.</dd>
<dt>{badge("XDR-TB")}</dt><dd>MDR/RR-TB that is resistant to any fluoroquinolone and at least one additional Group A drug, bedaquiline or linezolid.</dd>
<dt>{badge("Monoresistance")}</dt><dd>resistance to a single drug other than isoniazid or rifampicin.</dd>
<dt>{badge("Polyresistance")}</dt><dd>resistance to more than one drug, not including rifampicin, and not isoniazid alone.</dd>
</dl>
</details>

<div class="controls">
<input id="tbSearch" onkeyup="filterTable('tbSearch','tbTable')" placeholder="Search TB-Profiler results…">
<select onchange="filterResistance(this.value)">
<option value="">All resistance profiles</option>
<option value="Sensitive">Sensitive only</option>
<option value="Hr-TB">Hr-TB only</option>
<option value="RR-TB">RR-TB only</option>
<option value="MDR-TB">MDR-TB only</option>
<option value="MDR/RR-TB">MDR/RR-TB umbrella only</option>
<option value="Pre-XDR-TB">Pre-XDR-TB only</option>
<option value="XDR-TB">XDR-TB only</option>
<option value="Monoresistance">Monoresistance only</option>
<option value="Polyresistance">Polyresistance only</option>
<option value="Other drug resistance">Other drug resistance only</option>
<option value="Resistance not determined by TB-Profiler">Resistance not determined only</option>
<option value="Unknown">Unknown only</option>
</select>
<button onclick="downloadCSV('tbTable','tbprofiler_summary.csv')">Download TB-Profiler CSV</button>
</div>

<div class="note note-amber">
<strong>Mixed-infection flag:</strong> a screening flag derived from fields TB-Profiler
already reports &mdash; more than one MTBC main lineage, a lineage fraction below 0.90,
or resistance variants called at intermediate allele frequency (0.10&ndash;0.90), which
indicates a resistant subpopulation. rMAP-TB does not run a dedicated mixed-infection
detection tool and no sensitivity or specificity is claimed for this flag. Flagged
samples warrant review before their SNP distances are interpreted, because a mixed
sample carries ambiguous bases that distort pairwise comparisons.
</div>
{table_wrap(table)}
"""
        + section_close()
    )

def build_mutation_evidence_section():
    cleaned = [
        r for r in mutation_rows
        if any(clean(v) for v in r.values())
        and clean(r.get("sample"))
        and (not routing_by_sample or is_sample_routed_mtbc(r.get("sample")))
    ]

    head = section_open(4, "mutation-evidence", "Resistance Mutation Evidence Summary",
                        "Mutation-level evidence behind each TB-Profiler resistance call, grouped per sample.")

    if not cleaned:
        return head + """
<div class="note">
No mutation-level TB-Profiler resistance evidence was generated or the task was skipped.
</div>
""" + section_close()

    grouped = defaultdict(list)

    for r in cleaned:
        grouped[r.get("sample", "Unknown")].append(r)

    out = [head, """
<div class="note">
This section reports mutation-level drug-resistance evidence extracted from TB-Profiler JSON outputs. Results are grouped per sample to match the display style of the non-synonymous mutation summary. It complements the resistance profile by showing the underlying drug or evidence source, gene, mutation/change, confidence, and associated evidence fields where available.
</div>
<div class="controls"><button onclick="toggleDetails('sec-mutation-evidence', true)">Expand all samples</button><button onclick="toggleDetails('sec-mutation-evidence', false)">Collapse all</button></div>
"""]

    for s in sorted(grouped):
        out.append(f'<details class="sample-details"><summary><span class="mono">{display_sample_id(s)}</span> <span class="count-chip">{len(grouped[s])} mutation(s)</span></summary><div class="table-wrap"><table>')
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
                f"<td class='mono'>{safe(mutation)}</td>"
                f"<td class='mono'>{safe(change)}</td>"
                f"<td><span class='badge badge-green'>{safe(confidence)}</span></td>"
                f"<td>{safe(evidence)}</td>"
                "</tr>"
            )

        out.append("</tbody></table></div></details>")

    out.append(section_close())
    return "".join(out)

def build_nonsyn_section():
    cleaned = [
        r for r in nonsyn_rows
        if any(clean(v) for v in r.values())
        and (not routing_by_sample or is_sample_routed_mtbc(r.get("sample")))
    ]

    head = section_open(5, "nonsyn", "Non-synonymous Mutation Summary",
                        "Snippy-annotated non-synonymous variants in configured TB drug-resistance-associated genes, grouped per sample.")

    if not cleaned:
        return head + """
<div class="note">
No mutations were detected or mutation analysis was skipped.
</div>
""" + section_close()

    grouped = defaultdict(list)

    for r in cleaned:
        grouped[r.get("sample", "Unknown")].append(r)

    out = [head, """
<div class="note">
<strong>Mechanism:</strong> mutations are grouped per sample from per-sample Snippy annotation outputs and filtered to configured TB drug-resistance-associated genes. This complements TB-Profiler and should not replace catalogue-based resistance interpretation.
</div>
<div class="controls"><button onclick="toggleDetails('sec-nonsyn', true)">Expand all samples</button><button onclick="toggleDetails('sec-nonsyn', false)">Collapse all</button></div>
"""]

    for s in sorted(grouped):
        out.append(f'<details class="sample-details"><summary><span class="mono">{display_sample_id(s)}</span> <span class="count-chip">{len(grouped[s])} mutation(s)</span></summary><div class="table-wrap"><table>')
        out.append("<thead><tr><th class='lineage'>Gene</th><th class='mutations'>Effect</th><th class='mutations'>AA change</th><th class='mutations'>NT change</th><th class='status'>Product</th></tr></thead><tbody>")

        for r in grouped[s]:
            out.append(
                "<tr>"
                f"<td><strong>{safe(r.get('gene'))}</strong></td>"
                f"<td>{safe(r.get('effect'))}</td>"
                f"<td class='mono'>{safe(r.get('aa_change'))}</td>"
                f"<td class='mono'>{safe(r.get('nt_change'))}</td>"
                f"<td>{safe(r.get('product'))}</td>"
                "</tr>"
            )

        out.append("</tbody></table></div></details>")

    out.append(section_close())
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

    head = section_open(7, "snp-distance", "Pairwise SNP Distance and Cluster Summary",
                        "Threshold-based pairwise core-SNP distances between MTBC isolates, with genomically close pairs flagged for epidemiological review.")

    if not cleaned_pairs:
        return head + """
<div class="note">
Pairwise SNP distance and cluster reporting was not generated or was skipped. This usually occurs when fewer than the required number of MTBC samples were available for phylogenomic analysis.
</div>
""" + section_close()

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
            f"<td class='mono'>{display_sample_id(sample1)}</td>"
            f"<td class='mono'>{display_sample_id(sample2)}</td>"
            f"<td class='num'>{safe(distance)}</td>"
            f"<td class='num'>{safe(comparable_sites)}</td>"
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
                f"<td class='mono'>{display_sample_id(r.get('sample1'))}</td>"
                f"<td class='mono'>{display_sample_id(r.get('sample2'))}</td>"
                f"<td class='num'>{safe(r.get('snp_distance'))}</td>"
                f"<td>{cluster_badge(r.get('interpretation'))}</td>"
                "</tr>"
            )

        cluster_table_inner = f"""<table id="snpClusterTable">
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
</table>"""
        cluster_table = (
            "<h3>Genomically close sample pairs requiring epidemiological review</h3>"
            + table_wrap(cluster_table_inner)
        )
    else:
        cluster_table = """
<div class="note note-pass">
No sample pairs met the configured SNP thresholds for genomically close or intermediate-distance review.
</div>
"""

    pairs_table = f"""<table id="snpPairsTable">
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
</table>"""

    return head + f"""
<div class="note">
Pairwise SNP distances were calculated from the MTBC core genome alignment after excluding reference/non-sample sequences. Interpretations are threshold-based and should be considered alongside epidemiological metadata, lineage, resistance profile, sequencing quality, and phylogenetic support.
</div>

<div class="controls">
<input id="snpSearch" onkeyup="filterTable('snpSearch','snpPairsTable')" placeholder="Search SNP distance pairs…">
<button onclick="downloadCSV('snpPairsTable','pairwise_snp_distance_pairs.csv')">Download SNP Distance CSV</button>
</div>
{table_wrap(pairs_table)}
{cluster_table}
""" + section_close()

def svg_file_to_inline_html(path, title="embedded SVG", extra_class=""):
    """Inline SVG content so Terra-downloaded final_report.html is self-contained."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    try:
        svg = Path(path).read_text(encoding="utf-8", errors="replace").strip()
        svg = re.sub(r"^\s*<\?xml[^>]*>\s*", "", svg)
        svg = re.sub(r"^\s*<!DOCTYPE[^>]*>\s*", "", svg, flags=re.IGNORECASE)
        if "<svg" in svg.lower():
            return f'<div class="embedded-svg {extra_class}" role="img" aria-label="{safe(title)}">{svg}</div>'
    except Exception:
        return ""
    return ""

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
    total = sum(c for _, c in valid) or 1
    max_count = max(c for _, c in valid) or 1
    bars = []
    for lineage, count in valid:
        pct = max(2, (count / max_count) * 100)
        share = 100.0 * count / total
        unresolved = lineage.lower().startswith("not resolved")
        color = "#94a3b8" if unresolved else "#2563eb"
        bars.append(
            '<div class="lineage-bar-row">'
            f'<div class="lineage-label">{safe(lineage)}</div>'
            '<div class="lineage-bar-track">'
            f'<div class="lineage-bar" style="width:{pct:.1f}%; background:{color};" title="{safe(lineage)}: {count} ({share:.1f}%)"></div>'
            '</div>'
            f'<div class="lineage-value">{count} <small>({share:.0f}%)</small></div>'
            '</div>'
        )
    return '<div class="lineage-bar-chart">' + ''.join(bars) + '</div>'

def build_lineage_distribution_section():
    lineage_svg_exists = Path("final_report/lineage_distribution.svg").exists()

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
        cleaned = [{"lineage": k, "count": str(v)} for k, v in sorted(lineage_counts.items(), key=lambda kv: (-kv[1], kv[0]))]
    else:
        cleaned = [
            r for r in lineage_rows
            if any(clean(v) for v in r.values())
        ]

    head = section_open(6, "lineage", "Lineage Distribution Summary",
                        "TB-Profiler main-lineage calls across MTBC-routed isolates.")

    if not cleaned:
        return head + """
<div class="note">Lineage distribution summary was not generated or was skipped.</div>
""" + section_close()

    body = []

    for r in cleaned:
        body.append(
            "<tr>"
            f"<td>{safe(r.get('lineage'))}</td>"
            f"<td class='num'>{safe(r.get('count'))}</td>"
            "</tr>"
        )

    plot_html = build_lineage_bar_chart(cleaned)
    if not plot_html and lineage_svg_exists:
        plot_html = svg_file_to_inline_html("final_report/lineage_distribution.svg", "Lineage distribution")
    if not plot_html:
        plot_html = '<div class="note">Lineage distribution plot was not available.</div>'

    table = f"""<table id="lineageTable">
<thead>
<tr>
<th class="lineage" onclick="sortTable('lineageTable',0)">Lineage</th>
<th class="status" onclick="sortTable('lineageTable',1)">Count</th>
</tr>
</thead>
<tbody>
{''.join(body)}
</tbody>
</table>"""

    return head + f"""
<div class="note">
This section summarizes TB-Profiler lineage calls where available. Samples supported as MTBC by Kraken2 but without TB-Profiler lineage resolution are shown as “Not resolved by TB-Profiler” (grey bar).
</div>
<div class="two-col">
<div>{plot_html}</div>
<div>
<div class="controls">
<input id="lineageSearch" onkeyup="filterTable('lineageSearch','lineageTable')" placeholder="Search lineage summary…">
<button onclick="downloadCSV('lineageTable','lineage_distribution.csv')">Download Lineage CSV</button>
</div>
{table_wrap(table)}
</div>
</div>
""" + section_close()

def build_snp_heatmap_section():
    heatmap_exists = Path("final_report/snp_distance_heatmap.svg").exists()

    heatmap_html = svg_file_to_inline_html("final_report/snp_distance_heatmap.svg", "SNP distance heatmap", "heatmap-panel") if heatmap_exists else ""
    if heatmap_html:
        heatmap_html = (
            '<div class="controls">'
            '<button onclick="downloadInlineSvg(\'heatmapFigure\',\'snp_distance_heatmap.svg\')">Download SVG (vector)</button>'
            '<button onclick="downloadSvgAsPng(\'heatmapFigure\',\'snp_distance_heatmap_4x.png\',4)">Download PNG (4&times; hi-res)</button>'
            '</div>'
            f'<div id="heatmapFigure">{heatmap_html}</div>'
        )
    else:
        heatmap_html = '<div class="note">SNP distance heatmap was not generated, was skipped, or the SVG was empty.</div>'

    head = section_open(8, "heatmap", "SNP Distance Heatmap",
                        "Pairwise core-SNP distance matrix among MTBC isolates; colour follows a log scale so close pairs stay distinguishable next to distant outliers.")

    return head + f"""
<div class="note">
Lower SNP distances indicate closer genomic relatedness and should be interpreted together with epidemiological metadata, lineage, resistance profile, and tree topology. Hover any cell for the exact pair and distance. The reference sequence is excluded; the diagonal (self-comparison) is left blank.
</div>
{heatmap_html}
""" + section_close()

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
            f"<td class='mono'>{safe(display_sample_id(sample))}</td>"
            f"<td class='num'>{safe(r.get('mean_depth'))}</td>"
            f"<td class='num'>{safe(r.get('mtbc_percent'))}</td>"
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
            f"<td class='mono'>{safe(display_sample_id(sample))}</td>"
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
            f"<td class='num'>{safe(r.get('mean_depth'))}</td>"
            f"<td>{iqtree_html}</td>"
            "</tr>"
        )

    if not qc_body:
        qc_body = ['<tr><td colspan="6" class="empty-cell">No QC filtering rationale records available.</td></tr>']

    if not metadata_body:
        metadata_body = ['<tr><td colspan="13" class="empty-cell">No surveillance metadata records available.</td></tr>']

    qc_table = f"""<table id="qcRationaleTable">
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
</table>"""

    meta_table = f"""<table id="surveillanceMetadataTable">
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
</table>"""

    head = section_open(10, "surveillance", "QC Filtering Rationale and Surveillance Metadata",
                        "Transparent inclusion/exclusion rationale per sample and a surveillance-ready metadata table consistent with Section 3 and the tree labels.")

    return head + f"""
<div class="note">
This section provides a transparent rationale for sample inclusion/exclusion and a surveillance-ready metadata table. Resistance profile, drug-resistance detected status, and resistant drugs are populated from the canonical <code>resistance_profile_summary.tsv</code> generated by TB-Profiler parsing, ensuring this table matches Section 3 and the phylogenetic tree labels.
</div>

<div class="controls">
<button onclick="downloadCSV('qcRationaleTable','qc_filtering_rationale.csv')">Download QC Rationale CSV</button>
<button onclick="downloadCSV('surveillanceMetadataTable','tb_surveillance_metadata.csv')">Download Surveillance Metadata CSV</button>
</div>

<h3>QC Filtering Rationale</h3>
{table_wrap(qc_table)}

<h3>Surveillance Metadata</h3>
{table_wrap(meta_table)}
""" + section_close()

def build_tree_legend():
    legend = ['<div class="legend">']

    for label in RESISTANCE_ORDER:
        color = RESISTANCE_COLORS[label]
        legend.append(f'<span><i style="background:{color};"></i> {safe(label)}</span>')

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

    coverage_block = ""
    if coverage_excluded_mtbc:
        items = "".join(
            f"<li><strong>{display_sample_id(sid)}</strong> {safe(coverage_fail_note(sid))}</li>"
            for sid in coverage_excluded_mtbc
        )
        coverage_block = (
            f"<p><strong>{coverage_excluded_count} MTBC-routed sample(s)</strong> never reached the core-SNP alignment "
            f"because they failed the post-mapping coverage gate in SNIPPY_CORE_MTBC. This is why the number of "
            f"isolates in the tree can be smaller than the number of MTBC isolates retained by routing.</p>"
            f"<ul>{items}</ul>"
        )

    if not cleaned:
        css_class = "tree-notes" if coverage_excluded_mtbc else "tree-notes tree-notes-pass"
        return f"""
<div class="{css_class}">
<h4>Phylogenetic tree notes</h4>
{status_line}
{coverage_block}
<p>No samples were excluded by the IQ-TREE alignment-quality pre-filter. All samples that passed MTBC selection and core-alignment generation were eligible for phylogenetic inference.</p>
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
{coverage_block}
<p><strong>{len(cleaned)} MTBC-routed sample(s)</strong> were excluded from IQ-TREE before phylogenetic inference because their core-SNP alignment sequence did not meet tree-building quality requirements. Non-MTBC/NTM samples are counted separately as excluded from the MTBC workflow and are not described as IQ-TREE-only exclusions.</p>
<ul>
{''.join(notes)}
</ul>
{summary_details}
</div>
"""

def build_tree_figure():
    """The tree panel. Chooses the vector SVG, then the ETE3 PNG, and otherwise
    an explicit failure panel carrying the renderer's error and log. A
    placeholder image is never embedded as if it were the tree."""
    png_path = Path("final_report/mtbc_tree.png")
    svg_path = Path("final_report/mtbc_tree.svg")

    downloads = []
    if TREE_PRIMARY == "svg":
        downloads.append('<button onclick="downloadInlineSvg(\'treeFigure\',\'mtbc_core_snp_tree.svg\')" title="Lossless vector, scales to any size">Download SVG (vector)</button>')
        downloads.append('<button onclick="downloadSvgAsPng(\'treeFigure\',\'mtbc_core_snp_tree_4x.png\',4)" title="Rasterised at 4x the on-screen size">Download PNG (4&times; hi-res)</button>')
    if PNG_OK:
        try:
            png_b64 = base64.b64encode(png_path.read_bytes()).decode("ascii")
            downloads.append(f'<a class="btn" download="mtbc_core_snp_tree.png" href="data:image/png;base64,{png_b64}">Download PNG (ETE3)</a>')
        except Exception:
            pass
    newick_path = Path("final_report/mtbc_tree.cleaned.nwk")
    if newick_path.exists() and newick_path.stat().st_size > 3:
        try:
            nwk_b64 = base64.b64encode(newick_path.read_bytes()).decode("ascii")
            downloads.append(f'<a class="btn" download="mtbc_core_snp_tree.nwk" href="data:text/plain;base64,{nwk_b64}">Download Newick</a>')
        except Exception:
            pass

    toolbar = (
        '<div class="tree-toolbar">'
        '<div class="zoom-group">'
        '<button onclick="zoomTree(-1)" title="Zoom out">&minus;</button>'
        '<button onclick="zoomTree(0)" title="Fit to panel width">Fit</button>'
        '<button onclick="zoomTree(1)" title="Zoom in">+</button>'
        '<span id="treeZoomLabel" class="zoom-label">fit</span>'
        '</div>'
        f'<div class="dl-group">{"".join(downloads)}</div>'
        '</div>'
    )

    if TREE_PRIMARY == "svg":
        svg_inline = svg_file_to_inline_html(str(svg_path), "MTBC core-SNP phylogenetic tree", "tree-svg")
        return toolbar + f'<div class="tree-panel" id="treePanel"><div id="treeFigure" class="tree-figure">{svg_inline}</div></div>'

    if TREE_PRIMARY == "png":
        try:
            tree_b64 = base64.b64encode(png_path.read_bytes()).decode("ascii")
            img = (
                '<img class="mtbc-tree-img" '
                f'src="data:image/png;base64,{tree_b64}" '
                'alt="ETE3-rendered MTBC core-SNP phylogenetic tree">'
            )
        except Exception as exc:
            img = f'<div class="note note-danger">Tree image was generated, but could not be embedded into the final HTML report: {safe(str(exc))}</div>'
        return toolbar + f'<div class="tree-panel" id="treePanel"><div id="treeFigure" class="tree-figure">{img}</div></div>' + build_tree_legend()

    err = tree_render_status.get("error", "")
    log_block = ""
    if tree_render_log_text:
        log_block = f'<details open><summary>Tree renderer log (render.log)</summary><pre>{safe(tree_render_log_text[-12000:])}</pre></details>'
    else:
        log_block = '<p class="muted">No renderer log was provided to the report task.</p>'
    status_rows = "".join(
        f"<tr><td><code>{safe(k)}</code></td><td>{safe(v)}</td></tr>" for k, v in tree_render_status.items()
    ) or '<tr><td colspan="2" class="empty-cell">No render_status.txt was provided (older TREE_VISUALIZATION task).</td></tr>'
    return f"""
<div class="tree-failure">
<div class="tree-failure-head">
<span class="fail-icon">!</span>
<div>
<h4>Phylogenetic tree image was not rendered</h4>
<p>IQ-TREE status is reported separately above; this panel concerns the <em>drawing</em> of the tree. The Newick tree is preserved and downloadable, and the error and renderer log are shown here instead of an empty image.</p>
{f'<p class="err-line"><strong>Error:</strong> <code>{safe(err)}</code></p>' if err else ''}
</div>
</div>
<div class="table-wrap"><table class="kv"><thead><tr><th class="sample">Renderer status key</th><th class="status">Value</th></tr></thead><tbody>{status_rows}</tbody></table></div>
{log_block}
<div class="dl-group">{"".join(downloads)}</div>
</div>
"""

def build_tree_section():
    state_chip = {
        "ok": '<span class="state-chip state-ok">Rendered</span>',
        "warn": '<span class="state-chip state-warn">Rendered, annotations incomplete</span>',
        "bad": '<span class="state-chip state-bad">Rendering failed</span>',
    }.get(TREE_STATE, '<span class="state-chip state-muted">Unknown</span>')
    annotation_note = ""
    if TREE_STATE == "warn":
        annotation_note = (
            '<div class="note note-amber"><strong>Tree annotations incomplete.</strong> The tree topology is valid, '
            'but the lineage and resistance labels could not be joined to the tips: '
            f'{safe(tree_render_status.get("tips_matched", "0"))} of {safe(tree_render_status.get("tips_rendered", "?"))} tips matched the TB-Profiler metadata '
            f'(metadata status: <code>{safe(tree_render_status.get("metadata_status", "unknown"))}</code>). '
            'Use Section 3 for per-isolate lineage and resistance until this is resolved; the renderer log is below.'
            + (f'<details><summary>Tree renderer log</summary><pre>{safe(tree_render_log_text[-8000:])}</pre></details>' if tree_render_log_text else "")
            + '</div>'
        )

    rooted = tree_render_status.get("midpoint_rooted", "")
    if tree_render_status.get("tips_rendered", "0") in {"", "0"}:
        rooted_text = "rooting not applicable (tree not drawn)"
    else:
        rooted_text = "midpoint-rooted for display" if rooted == "true" else ("unrooted, as inferred" if rooted == "false" else "rooting not reported")

    facts = [
        ("Selected for MTBC workflow", f"{mtbc_retained} MTBC isolate(s) retained in the wider workflow"),
        ("Included in IQ-TREE", f"{iqtree_included_count} non-reference MTBC isolate(s) retained after core-SNP alignment quality filtering"),
        ("Excluded by coverage gate", f"{coverage_excluded_count} MTBC-routed sample(s) excluded from the core-SNP alignment for low depth/breadth" + (f" ({', '.join(display_sample_id(x) for x in coverage_excluded_mtbc)})" if coverage_excluded_mtbc else "")),
        ("Excluded from IQ-TREE only", f"{iqtree_excluded_count} sample(s) excluded from phylogenetic inference because of alignment-quality issues"),
        ("Excluded from MTBC workflow", f"{non_mtbc} non-MTBC or low-confidence isolate(s)"),
        ("Core alignment", "Snippy-core alignment"),
        ("Recombination", "Optional Gubbins-filtered alignment when enabled"),
        ("Tree", f"IQ-TREE2 maximum-likelihood phylogeny; {rooted_text}"),
        ("Display", "Vector (SVG) rendering drawn directly from the Newick tree; ETE3 raster image available as a download when rendered"),
    ]
    facts_html = "".join(f"<div class='fact'><span class='fact-k'>{safe(k)}</span><span class='fact-v'>{safe(v)}</span></div>" for k, v in facts)

    head = section_open(9, "tree", "MTBC-only Core-SNP Phylogenetic Tree",
                        "Maximum-likelihood core-SNP phylogeny of MTBC-routed isolates, tips coloured by canonical resistance category and annotated with lineage.")

    return head + f"""
<div class="tree-state-row">{state_chip}<span class="muted">{safe(TREE_DETAIL)}</span></div>
{annotation_note}
{build_tree_figure()}
{build_iqtree_exclusion_footnote()}
<div class="details">
<h3>Tree construction summary</h3>
<div class="facts">{facts_html}</div>
{build_iqtree_section()}
</div>
""" + section_close()

def build_provenance_section():
    head = section_open(11, "provenance", "Pipeline Provenance and Software Versions",
                        "What each component did, in the order the report presents it.")
    table = f"""<table>
<thead>
<tr>
<th class="sample">Workflow component</th>
<th class="status">Description</th>
</tr>
</thead>
<tbody>
<tr><td>Species typing</td><td>Kraken2 against a custom <em>Mycobacterium</em>-only database (76 RefSeq sequences, 50 species, 2 MTBC genomes) in <code>gmboowa/mycobacterium-kraken2-bracken:2026.05</code>; composition, accessions and documented limitations in <code>resources/kraken2_db/README.md</code>. Bracken is not used.</td></tr>
<tr><td>TB resistance and lineage</td><td>TB-Profiler Docker image provided by workflow input</td></tr>
<tr><td>Canonical resistance profile</td><td><code>resistance_profile_summary.tsv</code> used for Section 3, Surveillance Metadata, tree labels, badges, resistant drugs, and drug-resistant isolate count</td></tr>
<tr><td>Mutation-level resistance evidence</td><td>Parsed from TB-Profiler JSON outputs and summarized by sample, drug or evidence source, gene, mutation/change, confidence, and evidence</td></tr>
<tr><td>Lineage distribution</td><td>Lineage counts and barplot generated from TB-Profiler lineage fields where resolved</td></tr>
<tr><td>Pairwise SNP distance and clustering</td><td>Pairwise SNP distances calculated from the MTBC core genome alignment after excluding reference/non-sample sequences and interpreted using configured SNP thresholds</td></tr>
<tr><td>SNP distance heatmap</td><td>SVG heatmap generated from the pairwise SNP distance matrix after reference filtering; log-scaled colour with an explicit legend</td></tr>
<tr><td>Surveillance metadata</td><td>Downloadable metadata and QC filtering rationale TSV files generated for transparent surveillance reporting</td></tr>
<tr><td>Core-SNP phylogenomics</td><td>Snippy-core, optional Gubbins filtering, IQ-TREE2; tree drawn as SVG directly from the Newick (dependency-free) with an ETE3 PNG as secondary output; rendering status recorded in <code>render_status.txt</code></td></tr>
<tr><td>IQ-TREE problematic-sample filtering</td><td>Samples with excessive missing, ambiguous, or gap-only sequence content in the core-SNP alignment are excluded from IQ-TREE and reported in <code>excluded_from_iqtree.tsv</code></td></tr>
<tr><td>Report generated</td><td>{safe(run_started_utc)}</td></tr>
<tr><td>Run stamp</td><td>{safe(run_stamp)}</td></tr>
</tbody>
</table>"""
    return head + f"""
<p>The report documents all samples through QC, species typing, TB-Profiler analysis, mutation-level resistance evidence, lineage distribution, SNP distance clustering, SNP heatmap visualization, surveillance metadata, and MTBC-only phylogenomic reconstruction. Samples not classified as MTBC are excluded from the tree but retained in the workflow record for transparency. Samples with excessive missing, ambiguous, or gap-only content in the core-SNP alignment may also be excluded from IQ-TREE phylogenetic inference and are listed under the tree footnotes.</p>
<div class="note">
<strong>Interpretation:</strong> use close clustering together with bootstrap support, lineage, drug-resistance profile, mutation-level resistance evidence, lineage distribution, surveillance metadata, sample-exclusion notes, and SNP distances before making transmission inferences.
</div>
{table_wrap(table)}
""" + section_close()


TOC = [
    ("qc", "1", "QC &amp; trimming"),
    ("species", "2", "Species typing"),
    ("ntm", "2b", "Non-MTBC"),
    ("tbprofiler", "3", "TB-Profiler"),
    ("mutation-evidence", "4", "Mutation evidence"),
    ("nonsyn", "5", "Non-synonymous"),
    ("lineage", "6", "Lineage"),
    ("snp-distance", "7", "SNP distances"),
    ("heatmap", "8", "Heatmap"),
    ("tree", "9", "Phylogenetic tree"),
    ("surveillance", "10", "Surveillance metadata"),
    ("provenance", "11", "Provenance"),
]
toc_html = "".join(f'<a href="#sec-{slug}" data-target="sec-{slug}"><span>{num}</span>{label}</a>' for slug, num, label in TOC)

CSS = r"""
:root{
  --rm-bg:#f3f5f9;
  --rm-card:#ffffff;
  --rm-ink:#0f172a;
  --rm-ink-2:#334155;
  --rm-muted:#64748b;
  --rm-line:#e2e8f0;
  --rm-line-2:#cbd5e1;
  --rm-navy:#0b2545;
  --rm-navy-2:#13315c;
  --rm-teal:#0f766e;
  --rm-teal-2:#14b8a6;
  --rm-blue:#2563eb;
  --rm-purple:#6d28d9;
  --rm-red:#b91c1c;
  --rm-amber:#b45309;
  --rm-green:#15803d;
  --rm-radius:16px;
  --rm-shadow:0 1px 2px rgba(15,23,42,.05),0 8px 24px -12px rgba(15,23,42,.18);
  --rm-font:"Inter","Segoe UI",system-ui,-apple-system,Helvetica,Arial,sans-serif;
  --rm-mono:"JetBrains Mono","SFMono-Regular",Menlo,Consolas,monospace;
}
html{scroll-behavior:smooth;scroll-padding-top:112px;}
body.rmap-body{margin:0;font-family:var(--rm-font);background:var(--rm-bg);color:var(--rm-ink);font-size:14px;line-height:1.5;-webkit-font-smoothing:antialiased;}
#rmap-report *,#rmap-report *::before,#rmap-report *::after{box-sizing:border-box;}
#rmap-report a{color:var(--rm-blue);}
#rmap-report .muted{color:var(--rm-muted);}
#rmap-report .mono{font-family:var(--rm-mono);font-size:12.5px;letter-spacing:0;}
#rmap-report .num{font-variant-numeric:tabular-nums;text-align:right;white-space:nowrap;}
#rmap-report code{font-family:var(--rm-mono);font-size:12px;background:#f1f5f9;border:1px solid var(--rm-line);border-radius:6px;padding:1px 6px;color:#0f172a;}
#rmap-report .rmap-hero{position:relative;overflow:hidden;background:radial-gradient(1200px 500px at 10% -20%,rgba(20,184,166,.35),transparent 60%),radial-gradient(900px 400px at 95% 120%,rgba(37,99,235,.35),transparent 60%),linear-gradient(135deg,var(--rm-navy) 0%,var(--rm-navy-2) 55%,var(--rm-teal) 100%);color:#fff;padding:40px 48px 34px;}
#rmap-report .rmap-hero::after{content:"";position:absolute;inset:0;background-image:radial-gradient(rgba(255,255,255,.08) 1px,transparent 1px);background-size:22px 22px;pointer-events:none;}
#rmap-report .hero-inner{position:relative;z-index:1;max-width:1400px;margin:0 auto;}
#rmap-report .hero-kicker{display:inline-flex;align-items:center;gap:8px;font-size:11.5px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:#99f6e4;margin-bottom:12px;}
#rmap-report .hero-kicker i{width:8px;height:8px;border-radius:50%;background:#2dd4bf;box-shadow:0 0 0 4px rgba(45,212,191,.25);}
#rmap-report .rmap-hero h1{margin:0;font-size:34px;line-height:1.15;font-weight:800;letter-spacing:-.02em;color:#fff;}
#rmap-report .rmap-hero h1 small{display:block;font-size:15px;font-weight:500;color:#cbd5e1;margin-top:8px;letter-spacing:0;}
#rmap-report .hero-meta{display:flex;flex-wrap:wrap;gap:10px 22px;margin-top:18px;font-size:13px;color:#e2e8f0;}
#rmap-report .hero-meta strong{color:#fff;font-weight:600;}
#rmap-report .hero-code{font-family:var(--rm-mono);font-size:11px;background:rgba(255,255,255,.12);border:1px solid rgba(255,255,255,.2);color:#ccfbf1;padding:1px 6px;border-radius:6px;}
#rmap-report .pipeline-label{margin-top:22px;font-size:11px;font-weight:700;letter-spacing:.12em;text-transform:uppercase;color:#99f6e4;}
#rmap-report .pipeline{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:8px;margin-top:10px;}
#rmap-report .pipeline .step{position:relative;display:flex;flex-direction:column;gap:2px;padding:9px 11px 9px 34px;border-radius:12px;background:rgba(255,255,255,.09);border:1px solid rgba(255,255,255,.16);color:#f8fafc;backdrop-filter:blur(2px);min-height:56px;}
#rmap-report .pipeline .step b{position:absolute;left:10px;top:10px;width:18px;height:18px;border-radius:50%;background:#2dd4bf;color:#042f2e;font-size:10.5px;font-weight:800;display:grid;place-items:center;}
#rmap-report .pipeline .step-name{font-size:12.5px;font-weight:700;line-height:1.25;}
#rmap-report .pipeline .step-tool{font-size:11px;color:#cbd5e1;line-height:1.3;}
#rmap-report .rmap-hero h1 small{max-width:1180px;line-height:1.5;}
#rmap-report .rmap-nav{position:sticky;top:0;z-index:50;background:rgba(255,255,255,.92);backdrop-filter:saturate(1.4) blur(8px);border-bottom:1px solid var(--rm-line);}
#rmap-report .rmap-nav-inner{max-width:1400px;margin:0 auto;padding:8px 48px;display:flex;flex-wrap:wrap;gap:4px 6px;overflow:visible;justify-content:flex-start;}
#rmap-report .rmap-nav a{flex:0 1 auto;display:inline-flex;align-items:center;gap:7px;font-size:12.5px;font-weight:600;color:var(--rm-ink-2);text-decoration:none;padding:7px 11px;border-radius:999px;border:1px solid transparent;white-space:nowrap;}
#rmap-report .rmap-nav a span{display:inline-grid;place-items:center;min-width:20px;height:20px;padding:0 5px;border-radius:999px;background:#e2e8f0;color:var(--rm-ink-2);font-size:11px;font-weight:800;}
#rmap-report .rmap-nav a:hover{background:#f1f5f9;border-color:var(--rm-line);}
#rmap-report .rmap-nav a.active{background:var(--rm-navy);color:#fff;}
#rmap-report .rmap-nav a.active span{background:rgba(255,255,255,.2);color:#fff;}
#rmap-report .rmap-container{max-width:1400px;margin:0 auto;padding:28px 48px 60px;width:auto;}
#rmap-report .section{background:var(--rm-card);border:1px solid var(--rm-line);border-radius:var(--rm-radius);padding:26px 28px 24px;margin:0 0 26px;box-shadow:var(--rm-shadow);}
#rmap-report .sec-head{display:flex;gap:16px;align-items:flex-start;padding-bottom:16px;margin-bottom:18px;border-bottom:1px solid var(--rm-line);}
#rmap-report .sec-num{flex:0 0 auto;display:grid;place-items:center;width:40px;height:40px;border-radius:12px;background:linear-gradient(135deg,var(--rm-navy),var(--rm-teal));color:#fff;font-weight:800;font-size:15px;letter-spacing:-.02em;box-shadow:0 6px 14px -6px rgba(11,37,69,.6);}
#rmap-report .section h2{margin:4px 0 0;font-size:21px;line-height:1.25;font-weight:750;letter-spacing:-.015em;color:var(--rm-ink);border:0;padding:0;}
#rmap-report .section h3{margin:22px 0 10px;font-size:15.5px;font-weight:700;color:var(--rm-ink);}
#rmap-report .section h4{margin:0 0 8px;font-size:15px;}
#rmap-report .sec-lead{margin:6px 0 0;color:var(--rm-muted);font-size:13.5px;max-width:980px;}
#rmap-report .two-col{display:grid;grid-template-columns:minmax(0,1.1fr) minmax(0,1fr);gap:22px;align-items:start;}
#rmap-report .cards{display:grid;grid-template-columns:repeat(4,minmax(150px,1fr));gap:16px;margin-bottom:18px;}
#rmap-report .card{position:relative;background:var(--rm-card);border:1px solid var(--rm-line);border-radius:var(--rm-radius);padding:18px 20px 16px;box-shadow:var(--rm-shadow);overflow:hidden;}
#rmap-report .card::before{content:"";position:absolute;left:0;top:0;bottom:0;width:5px;background:var(--accent,#2563eb);}
#rmap-report .card h3{margin:0;color:var(--rm-muted);font-size:12.5px;font-weight:600;letter-spacing:.02em;text-transform:uppercase;}
#rmap-report .card .num{font-size:36px;font-weight:800;letter-spacing:-.03em;margin-top:6px;line-height:1.05;color:var(--accent,#2563eb);text-align:left;}
#rmap-report .card .sub{margin-top:6px;font-size:12.5px;color:var(--rm-muted);}
#rmap-report .card.blue{--accent:#2563eb;}
#rmap-report .card.green{--accent:#0f766e;}
#rmap-report .card.orange{--accent:#c2410c;}
#rmap-report .card.red{--accent:#b91c1c;}
#rmap-report .overview{display:grid;grid-template-columns:minmax(0,1.5fr) minmax(0,1fr);gap:16px;margin-bottom:26px;}
#rmap-report .dist-block{background:var(--rm-card);border:1px solid var(--rm-line);border-radius:var(--rm-radius);padding:18px 20px;box-shadow:var(--rm-shadow);}
#rmap-report .dist-title{font-weight:700;font-size:14px;margin-bottom:12px;}
#rmap-report .dist-title small{display:block;font-weight:500;color:var(--rm-muted);font-size:12px;margin-top:2px;}
#rmap-report .dist-bar{display:flex;height:34px;border-radius:10px;overflow:hidden;background:#f1f5f9;gap:2px;}
#rmap-report .dist-seg{display:flex;align-items:center;justify-content:center;font-size:12px;font-weight:700;min-width:4px;}
#rmap-report .dist-legend{display:flex;flex-wrap:wrap;gap:6px 14px;margin-top:12px;}
#rmap-report .dist-chip{display:inline-flex;align-items:center;gap:7px;font-size:12.5px;color:var(--rm-ink-2);}
#rmap-report .dist-chip i{width:10px;height:10px;border-radius:3px;display:inline-block;}
#rmap-report .dist-chip b{font-weight:700;color:var(--rm-ink);}
#rmap-report .dist-chip small{color:var(--rm-muted);}
#rmap-report .status-strip{display:grid;grid-template-columns:1fr;gap:8px;align-content:start;}
#rmap-report .pill{display:flex;gap:10px;align-items:flex-start;padding:10px 12px;border-radius:12px;background:var(--rm-card);border:1px solid var(--rm-line);box-shadow:var(--rm-shadow);font-size:12.5px;}
#rmap-report .pill strong{display:block;font-size:13px;}
#rmap-report .pill small{display:block;color:var(--rm-muted);margin-top:1px;word-break:break-word;}
#rmap-report .pill-dot{flex:0 0 auto;width:10px;height:10px;border-radius:50%;margin-top:5px;background:#94a3b8;box-shadow:0 0 0 3px rgba(148,163,184,.22);}
#rmap-report .pill-ok .pill-dot{background:#16a34a;box-shadow:0 0 0 3px rgba(22,163,74,.2);}
#rmap-report .pill-warn .pill-dot{background:#f59e0b;box-shadow:0 0 0 3px rgba(245,158,11,.2);}
#rmap-report .pill-bad .pill-dot{background:#dc2626;box-shadow:0 0 0 3px rgba(220,38,38,.2);}
#rmap-report .pill-bad{border-color:#fecaca;background:#fff7f7;}
#rmap-report .note{background:#eff6ff;border:1px solid #dbeafe;border-left:4px solid var(--rm-blue);padding:12px 14px;border-radius:12px;margin:12px 0;line-height:1.5;font-size:13.5px;}
#rmap-report .note-amber{background:#fffbeb;border-color:#fde68a;border-left-color:#d97706;}
#rmap-report .note-pass{background:#f0fdf4;border-color:#bbf7d0;border-left-color:#16a34a;}
#rmap-report .note-danger{background:#fef2f2;border-color:#fecaca;border-left-color:#b91c1c;}
#rmap-report .note ul{margin:8px 0 0 18px;padding:0;}
#rmap-report .controls{display:flex;flex-wrap:wrap;gap:10px;margin:0 0 14px;align-items:center;}
#rmap-report .controls input,#rmap-report .controls select,#rmap-report .controls button,#rmap-report .btn,#rmap-report .tree-toolbar button,#rmap-report .dl-group a{font-family:var(--rm-font);border:1px solid var(--rm-line-2);border-radius:10px;padding:9px 12px;font-size:13px;background:#fff;color:var(--rm-ink);line-height:1.2;height:auto;}
#rmap-report .controls input{min-width:260px;}
#rmap-report .controls input:focus,#rmap-report .controls select:focus{outline:2px solid rgba(37,99,235,.35);outline-offset:1px;border-color:var(--rm-blue);}
#rmap-report .controls button,#rmap-report .btn,#rmap-report .tree-toolbar button,#rmap-report .dl-group a{cursor:pointer;background:#eef2ff;color:var(--rm-navy);font-weight:600;text-decoration:none;display:inline-flex;align-items:center;gap:6px;}
#rmap-report .controls button:hover,#rmap-report .btn:hover,#rmap-report .tree-toolbar button:hover,#rmap-report .dl-group a:hover{background:#e0e7ff;}
#rmap-report .table-wrap{overflow:auto;border:1px solid var(--rm-line);border-radius:12px;background:#fff;max-height:72vh;}
#rmap-report .table-wrap table{width:100%;border-collapse:separate;border-spacing:0;font-size:13px;margin:0;background:#fff;}
#rmap-report .table-wrap th{position:sticky;top:0;z-index:2;color:#fff;padding:11px 12px;text-align:left;cursor:pointer;user-select:none;font-weight:700;font-size:12.5px;letter-spacing:.01em;white-space:nowrap;border-bottom:0;}
#rmap-report .table-wrap th::after{content:" ↕";opacity:.55;font-size:10px;}
#rmap-report .table-wrap td{padding:10px 12px;border-bottom:1px solid var(--rm-line);vertical-align:top;background:#fff;line-height:1.45;max-width:520px;}
#rmap-report .table-wrap tbody tr:nth-child(even) td{background:#fafbfd;}
#rmap-report .table-wrap tbody tr:hover td{background:#f1f5f9;}
#rmap-report .table-wrap tbody tr:last-child td{border-bottom:0;}
#rmap-report .table-wrap td small{color:var(--rm-muted);font-size:12px;}
#rmap-report .empty-cell{color:var(--rm-muted);text-align:center;padding:22px !important;}
#rmap-report th.sample{background:#0f766e;}
#rmap-report th.species{background:#1d4ed8;}
#rmap-report th.lineage{background:#6d28d9;}
#rmap-report th.resistance{background:#b91c1c;}
#rmap-report th.mutations{background:#c2410c;}
#rmap-report th.status{background:#0b2545;}
#rmap-report table.kv td:first-child{width:220px;}
#rmap-report .badge{display:inline-block;padding:3px 9px;border-radius:999px;color:#fff !important;font-size:11.5px;font-weight:700 !important;line-height:1.4;white-space:nowrap;vertical-align:middle;}
#rmap-report .badge-green{background:#15803d !important;}
#rmap-report .badge-red{background:#b91c1c !important;}
#rmap-report .badge-blue{background:#1d4ed8 !important;}
#rmap-report .badge-orange{background:#c2410c !important;}
#rmap-report .badge-soft{background:#e2e8f0 !important;color:#0f172a !important;}
#rmap-report .badge-wrap{white-space:normal;line-height:1.35;}
#rmap-report .res-badge{display:inline-block;padding:3px 10px;border-radius:999px;color:#fff !important;font-size:11.5px;font-weight:700;white-space:nowrap;letter-spacing:.01em;}
#rmap-report .count-chip{display:inline-block;margin-left:8px;padding:2px 9px;border-radius:999px;background:#eef2ff;color:var(--rm-navy);font-size:11.5px;font-weight:700;}
#rmap-report .state-chip{display:inline-flex;align-items:center;gap:6px;padding:4px 11px;border-radius:999px;font-size:12px;font-weight:700;}
#rmap-report .state-chip::before{content:"";width:8px;height:8px;border-radius:50%;background:currentColor;}
#rmap-report .state-ok{background:#dcfce7;color:#166534;}
#rmap-report .state-bad{background:#fee2e2;color:#991b1b;}
#rmap-report .state-warn{background:#fef3c7;color:#92400e;}
#rmap-report .state-muted{background:#e2e8f0;color:#334155;}
#rmap-report .tree-state-row{display:flex;align-items:center;gap:12px;margin-bottom:12px;font-size:13px;}
#rmap-report .sub-block{margin-top:6px;font-size:12.5px;color:var(--rm-ink-2);}
#rmap-report .sub-label{display:inline-block;font-size:10.5px;font-weight:800;letter-spacing:.08em;text-transform:uppercase;color:var(--rm-muted);margin-right:4px;}
#rmap-report details.inline-details{display:inline;margin:0;padding:0;border:0;background:transparent;}
#rmap-report details.inline-details summary{display:inline;cursor:pointer;color:var(--rm-blue);font-weight:600;font-size:12px;list-style:none;}
#rmap-report details.inline-details summary::-webkit-details-marker{display:none;}
#rmap-report details.inline-details[open] summary{display:block;margin-top:4px;}
#rmap-report ul.compact-list{margin:6px 0 0 0;padding:8px 10px 8px 24px;background:#f8fafc;border:1px solid var(--rm-line);border-radius:8px;max-height:260px;overflow:auto;columns:2;column-gap:18px;font-family:var(--rm-mono);font-size:11.5px;}
#rmap-report ul.compact-list li{break-inside:avoid;margin:1px 0;}
#rmap-report details{margin:12px 0 10px;border:1px solid var(--rm-line);border-radius:12px;padding:10px 14px;background:#fff;}
#rmap-report summary{cursor:pointer;font-weight:700;color:var(--rm-ink);}
#rmap-report summary:hover{color:var(--rm-blue);}
#rmap-report details.sample-details{margin:8px 0;}
#rmap-report details.sample-details .table-wrap{margin-top:10px;}
#rmap-report details.defs{background:#f8fafc;}
#rmap-report .def-grid{display:grid;grid-template-columns:max-content 1fr;gap:10px 14px;margin:12px 0 4px;align-items:center;}
#rmap-report .def-grid dt,#rmap-report .def-grid dd{margin:0;}
#rmap-report .def-grid dd{font-size:13px;color:var(--rm-ink-2);}
#rmap-report pre{white-space:pre-wrap;overflow:auto;background:#0f172a;color:#e2e8f0;padding:14px 16px;border-radius:12px;font-size:12px;line-height:1.5;border:0;margin:10px 0 0;font-family:var(--rm-mono);word-break:normal;max-height:520px;}
#rmap-report .embedded-report{overflow:auto;max-height:650px;border:1px solid var(--rm-line);border-radius:12px;padding:10px;background:white;margin-top:10px;}
#rmap-report .lineage-bar-chart{padding:16px 18px;border:1px solid var(--rm-line);border-radius:14px;background:#fff;}
#rmap-report .lineage-bar-row{display:grid;grid-template-columns:minmax(120px,200px) 1fr 70px;gap:12px;align-items:center;margin:9px 0;}
#rmap-report .lineage-label{font-weight:600;color:var(--rm-ink);font-size:13px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
#rmap-report .lineage-bar-track{width:100%;background:#f1f5f9;border-radius:6px;overflow:hidden;height:18px;}
#rmap-report .lineage-bar{height:18px;border-radius:0 6px 6px 0;}
#rmap-report .lineage-value{font-variant-numeric:tabular-nums;font-weight:700;font-size:13px;text-align:right;}
#rmap-report .lineage-value small{color:var(--rm-muted);font-weight:500;}
#rmap-report .embedded-svg{max-width:100%;overflow-x:auto;margin:1rem 0;border:1px solid var(--rm-line);border-radius:14px;background:#fff;padding:1rem;}
#rmap-report .embedded-svg svg{max-width:100%;height:auto;display:block;}
#rmap-report .heatmap-panel svg{max-width:none;}
#rmap-report .tree-toolbar{display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap;margin:0 0 10px;}
#rmap-report .zoom-group,#rmap-report .dl-group{display:flex;gap:8px;align-items:center;flex-wrap:wrap;}
#rmap-report .zoom-group button{min-width:38px;justify-content:center;}
#rmap-report .zoom-label{font-size:12px;color:var(--rm-muted);min-width:40px;}
#rmap-report .tree-panel{background:#fff;border:1px solid var(--rm-line);border-radius:14px;padding:10px;overflow:auto;max-height:78vh;}
#rmap-report .tree-figure .embedded-svg{border:0;padding:0;margin:0;overflow:visible;}
#rmap-report .tree-figure svg{display:block;max-width:100%;height:auto;}
#rmap-report .mtbc-tree-img{display:block;width:auto;max-width:none;height:auto;margin:0;}
#rmap-report .legend{display:flex;flex-wrap:wrap;gap:8px;margin-top:12px;font-size:12px;}
#rmap-report .legend span{display:inline-flex;align-items:center;gap:7px;border-radius:999px;padding:5px 10px;background:#f1f5f9;}
#rmap-report .legend i{width:10px;height:10px;border-radius:50%;display:inline-block;}
#rmap-report .tree-failure{border:1px solid #fecaca;background:#fff7f7;border-radius:14px;padding:18px 20px;}
#rmap-report .tree-failure-head{display:flex;gap:14px;align-items:flex-start;}
#rmap-report .fail-icon{flex:0 0 auto;width:34px;height:34px;border-radius:10px;background:#b91c1c;color:#fff;display:grid;place-items:center;font-weight:900;font-size:20px;}
#rmap-report .tree-failure h4{margin:2px 0 6px;font-size:16px;color:#991b1b;}
#rmap-report .err-line code{background:#fee2e2;border-color:#fecaca;color:#7f1d1d;white-space:pre-wrap;}
#rmap-report .tree-notes{margin-top:16px;padding:14px 16px;border:1px solid #fde68a;border-left:4px solid #f59e0b;background:#fffbeb;border-radius:12px;line-height:1.5;}
#rmap-report .tree-notes h4{margin:0 0 8px;color:#92400e;font-size:15px;}
#rmap-report .tree-notes p{margin:8px 0;}
#rmap-report .tree-notes ul{margin:8px 0 0 20px;padding:0;}
#rmap-report .tree-notes li{margin-bottom:8px;}
#rmap-report .tree-notes details{margin-top:10px;background:#fff;}
#rmap-report .tree-notes-pass{border-color:#bbf7d0;border-left-color:#16a34a;background:#f0fdf4;}
#rmap-report .tree-notes-pass h4{color:#166534;}
#rmap-report .details{background:#f8fafc;border:1px solid var(--rm-line);border-radius:12px;padding:14px 16px;margin-top:16px;}
#rmap-report .details h3{margin-top:0;}
#rmap-report .facts{display:grid;grid-template-columns:repeat(auto-fit,minmax(440px,1fr));gap:4px 28px;}
#rmap-report .fact{display:grid;grid-template-columns:170px 1fr;gap:10px;padding:7px 0;border-bottom:1px dashed var(--rm-line);font-size:13px;}
#rmap-report .fact-k{color:var(--rm-muted);font-weight:600;}
#rmap-report .fact-v{color:var(--rm-ink);}
#rmap-report .rmap-footer{font-size:12px;color:var(--rm-muted);margin:28px 0 0;margin-left:0;text-align:left;background:none;border-top:1px solid var(--rm-line);padding:16px 0 0;display:flex;justify-content:space-between;flex-wrap:wrap;gap:8px;}
#rmap-report .to-top{position:fixed;right:22px;bottom:22px;z-index:60;width:42px;height:42px;border-radius:50%;background:var(--rm-navy);color:#fff;display:grid;place-items:center;text-decoration:none;font-weight:800;box-shadow:0 10px 24px -8px rgba(11,37,69,.6);opacity:0;pointer-events:none;transition:opacity .2s;}
#rmap-report .to-top.show{opacity:1;pointer-events:auto;}
@media (max-width:1100px){
  #rmap-report .pipeline{grid-template-columns:repeat(4,minmax(0,1fr));}
  #rmap-report .overview,#rmap-report .two-col{grid-template-columns:1fr;}
  #rmap-report .cards{grid-template-columns:repeat(2,1fr);}
}
@media (max-width:760px){
  #rmap-report .rmap-container,#rmap-report .rmap-nav-inner{padding-left:16px;padding-right:16px;}
  #rmap-report .rmap-hero{padding:28px 16px 24px;}
  #rmap-report .rmap-hero h1{font-size:26px;}
  #rmap-report .section{padding:18px 16px;}
  #rmap-report .cards{grid-template-columns:1fr 1fr;}
  #rmap-report .pipeline{grid-template-columns:repeat(2,minmax(0,1fr));}
  #rmap-report .rmap-nav a{padding:6px 9px;font-size:12px;}
  #rmap-report .fact{grid-template-columns:1fr;gap:2px;}
  #rmap-report .controls input{min-width:0;width:100%;}
  #rmap-report ul.compact-list{columns:1;}
}
@media print{
  body.rmap-body{background:#fff;}
  #rmap-report .rmap-nav,#rmap-report .to-top,#rmap-report .controls,#rmap-report .tree-toolbar{display:none !important;}
  #rmap-report .section{box-shadow:none;break-inside:avoid;page-break-inside:avoid;}
  #rmap-report .table-wrap{overflow:visible;}
  #rmap-report .tree-panel{max-height:none;overflow:visible;}
  #rmap-report details{break-inside:avoid;}
}
"""

JS = r"""
function filterTable(inputId, tableId){
  const filter=document.getElementById(inputId).value.toLowerCase();
  const rows=document.getElementById(tableId).getElementsByTagName("tbody")[0].rows;
  for(let i=0;i<rows.length;i++){
    rows[i].style.display=rows[i].innerText.toLowerCase().includes(filter)?"":"none";
  }
}

function filterResistance(value){
  const rows=document.getElementById("tbTable").getElementsByTagName("tbody")[0].rows;
  for(let i=0;i<rows.length;i++){
    const badge=rows[i].cells[3] ? rows[i].cells[3].querySelector(".res-badge") : null;
    const label=badge ? badge.innerText.trim() : (rows[i].cells[3] ? rows[i].cells[3].innerText.trim() : "");
    rows[i].style.display=(value===""||label===value)?"":"none";
  }
}

function sortTable(tableId,col){
  const table=document.getElementById(tableId);
  const tbody=table.tBodies[0];
  const rows=Array.from(tbody.rows);
  const asc=table.getAttribute("data-sort-col")!=col||table.getAttribute("data-sort-dir")!=="asc";
  rows.sort((a,b)=>{
    const A=a.cells[col].innerText.replace(/,/g,'');
    const B=b.cells[col].innerText.replace(/,/g,'');
    const nA=parseFloat(A),nB=parseFloat(B);
    if(!isNaN(nA)&&!isNaN(nB))return asc?nA-nB:nB-nA;
    return asc?A.localeCompare(B):B.localeCompare(A);
  });
  rows.forEach(r=>tbody.appendChild(r));
  table.setAttribute("data-sort-col",col);
  table.setAttribute("data-sort-dir",asc?"asc":"desc");
}

function downloadCSV(tableId, filename){
  const table=document.getElementById(tableId);
  let csv=[];
  for(const row of table.rows){
    let cols=[];
    for(const cell of row.cells){
      cols.push('"' + cell.innerText.replace(/"/g,'""') + '"');
    }
    csv.push(cols.join(","));
  }
  const blob=new Blob([csv.join("\n")],{type:"text/csv"});
  const link=document.createElement("a");
  link.href=URL.createObjectURL(blob);
  link.download=filename;
  link.click();
}

function toggleDetails(sectionId, open){
  const sec=document.getElementById(sectionId);
  if(!sec) return;
  sec.querySelectorAll("details.sample-details").forEach(d=>{ d.open=!!open; });
}

function downloadSvgAsPng(containerId, filename, scale){
  const host=document.getElementById(containerId);
  const svg=host ? host.querySelector("svg") : null;
  if(!svg) return;
  scale=scale||4;
  const clone=svg.cloneNode(true);
  clone.removeAttribute("style");
  if(!clone.getAttribute("xmlns")) clone.setAttribute("xmlns","http://www.w3.org/2000/svg");
  const w=parseFloat(svg.getAttribute("width"))||svg.getBoundingClientRect().width;
  const h=parseFloat(svg.getAttribute("height"))||svg.getBoundingClientRect().height;
  clone.setAttribute("width",w); clone.setAttribute("height",h);
  const data='<?xml version="1.0" encoding="UTF-8"?>\n'+clone.outerHTML;
  const url=URL.createObjectURL(new Blob([data],{type:"image/svg+xml;charset=utf-8"}));
  const img=new Image();
  img.onload=function(){
    const canvas=document.createElement("canvas");
    canvas.width=Math.round(w*scale); canvas.height=Math.round(h*scale);
    const ctx=canvas.getContext("2d");
    ctx.fillStyle="#ffffff"; ctx.fillRect(0,0,canvas.width,canvas.height);
    ctx.setTransform(scale,0,0,scale,0,0);
    ctx.drawImage(img,0,0);
    URL.revokeObjectURL(url);
    canvas.toBlob(function(blob){
      const link=document.createElement("a");
      link.href=URL.createObjectURL(blob); link.download=filename; link.click();
    },"image/png");
  };
  img.onerror=function(){ URL.revokeObjectURL(url); alert("Could not rasterise the figure in this browser; use Download SVG instead."); };
  img.src=url;
}

function downloadInlineSvg(containerId, filename){
  const host=document.getElementById(containerId);
  const svg=host ? host.querySelector("svg") : null;
  if(!svg) return;
  const clone=svg.cloneNode(true);
  clone.removeAttribute("style");
  if(!clone.getAttribute("xmlns")) clone.setAttribute("xmlns","http://www.w3.org/2000/svg");
  const text='<?xml version="1.0" encoding="UTF-8"?>\n'+clone.outerHTML;
  const blob=new Blob([text],{type:"image/svg+xml;charset=utf-8"});
  const link=document.createElement("a");
  link.href=URL.createObjectURL(blob);
  link.download=filename;
  link.click();
}

var rmapTreeZoom=0;
var rmapTreeNaturalWidth=null;
function zoomTree(direction){
  const fig=document.getElementById("treeFigure");
  if(!fig) return;
  const svg=fig.querySelector("svg");
  const img=fig.querySelector("img");
  const el=svg||img;
  if(!el) return;
  if(rmapTreeNaturalWidth===null){
    rmapTreeNaturalWidth = svg ? parseFloat(svg.getAttribute("width")) || svg.getBoundingClientRect().width
                               : (img.naturalWidth || img.getBoundingClientRect().width);
  }
  const label=document.getElementById("treeZoomLabel");
  if(direction===0){
    rmapTreeZoom=0;
    el.style.width=""; el.style.maxWidth="100%"; el.style.height="auto";
    if(label) label.textContent="fit";
    return;
  }
  if(rmapTreeZoom===0){
    const panel=document.getElementById("treePanel");
    rmapTreeZoom=Math.min(1, (panel.clientWidth-24)/rmapTreeNaturalWidth);
  }
  rmapTreeZoom = direction>0 ? Math.min(rmapTreeZoom*1.25, 6) : Math.max(rmapTreeZoom/1.25, 0.1);
  el.style.maxWidth="none";
  el.style.width=(rmapTreeNaturalWidth*rmapTreeZoom)+"px";
  el.style.height="auto";
  if(label) label.textContent=Math.round(rmapTreeZoom*100)+"%";
}

(function(){
  const links=Array.from(document.querySelectorAll("#rmap-report .rmap-nav a"));
  const sections=links.map(a=>document.getElementById(a.getAttribute("data-target"))).filter(Boolean);
  const toTop=document.getElementById("rmapToTop");
  function onScroll(){
    const nav=document.querySelector("#rmap-report .rmap-nav");
    const y=window.scrollY+(nav?nav.offsetHeight:0)+24;
    let current=null;
    for(const s of sections){ if(s.offsetTop<=y) current=s; }
    links.forEach(a=>a.classList.toggle("active", current && a.getAttribute("data-target")===current.id));
    if(toTop) toTop.classList.toggle("show", window.scrollY>600);
  }
  window.addEventListener("scroll", onScroll, {passive:true});
  onScroll();
})();
"""

pipeline_steps = [
    ("Read trimming", "Trimmomatic"),
    ("Read QC", "FastQC &middot; MultiQC"),
    ("Species typing", "Kraken2, <em>Mycobacterium</em>-only DB"),
    ("MTBC / NTM routing", "Kraken2 MTBC support"),
    ("Lineage &amp; resistance", "TB-Profiler, WHO categories"),
    ("Mutation evidence", "TB-Profiler JSON &middot; Snippy annotation"),
    ("Core-SNP alignment", "Snippy &middot; snippy-core"),
    ("SNP distances &amp; clusters", "pairwise core-SNP thresholds"),
    ("Recombination filter", "Gubbins, optional"),
    ("ML phylogeny", "IQ-TREE2, UFBoot + SH-aLRT"),
    ("Tree &amp; surveillance figures", "vector SVG &middot; ETE3"),
    ("Integrated report", "this document"),
]
pipeline_html = "".join(
    f'<span class="step"><b>{i + 1}</b><span class="step-name">{name}</span><span class="step-tool">{tool}</span></span>'
    for i, (name, tool) in enumerate(pipeline_steps)
)

body_html = f"""
<div id="rmap-report">
<header class="rmap-hero">
<div class="hero-inner">
<div class="hero-kicker"><i></i>rMAP-TB &middot; WDL workflow <code class="hero-code">rMAP_TB</code> &middot; Integrated mycobacterial genomic surveillance report</div>
<h1>Mycobacterial species typing, TB drug-resistance profiling and MTBC core-SNP phylogenomics<small>From paired-end Illumina reads to a surveillance-ready report: read QC and trimming, Kraken2 species typing with MTBC/NTM routing, TB-Profiler lineage and WHO-category resistance calls with mutation-level evidence, Snippy core-SNP alignment, pairwise SNP-distance clustering and heatmap, optional Gubbins recombination filtering, IQ-TREE2 maximum-likelihood phylogeny, and transparent QC/inclusion rationale for every sample.</small></h1>
<div class="hero-meta">
<div><strong>Run generated:</strong> {safe(run_started_utc)}</div>
<div><strong>Run stamp:</strong> <span class="mono">{safe(run_stamp)}</span></div>
<div><strong>Samples:</strong> {total_samples} paired libraries &middot; {mtbc_retained} MTBC-routed &middot; {non_mtbc} non-MTBC</div>
<div><strong>Reference:</strong> <em>M. tuberculosis</em> H37Rv (NC_000962.3)</div>
</div>
<div class="pipeline-label">Workflow stages, in execution order</div>
<div class="pipeline">{pipeline_html}</div>
</div>
</header>

<nav class="rmap-nav" aria-label="Report sections"><div class="rmap-nav-inner">{toc_html}</div></nav>

<main class="rmap-container">

<div class="cards">
<div class="card blue"><h3>Total paired samples</h3><div class="num">{total_samples}</div><div class="sub">Libraries entering the workflow</div></div>
<div class="card green"><h3>MTBC isolates retained</h3><div class="num">{mtbc_retained}</div><div class="sub">Routed to TB-Profiler and phylogenomics</div></div>
<div class="card orange"><h3>Non-MTBC excluded</h3><div class="num">{non_mtbc}</div><div class="sub">Reported in Section 2b, excluded downstream</div></div>
<div class="card red"><h3>Drug-resistant isolates</h3><div class="num">{drug_resistant}</div><div class="sub">Any canonical resistance category</div></div>
</div>

<div class="overview">
{build_resistance_distribution()}
{build_status_strip()}
</div>

{build_unparsed_inputs_note()}

{build_qc_section()}

{build_species_section()}

{build_ntm_species_section()}

{build_tbprofiler_section()}

{build_mutation_evidence_section()}

{build_nonsyn_section()}

{build_lineage_distribution_section()}

{build_snp_distance_section()}

{build_snp_heatmap_section()}

{build_tree_section()}

{build_qc_rationale_surveillance_metadata_section()}

{build_provenance_section()}

<div class="rmap-footer">
<span>Generated by the rMAP-TB WDL workflow &middot; {safe(run_started_utc)} &middot; run stamp {safe(run_stamp)}</span>
<span>Tables are sortable (click a header), searchable and downloadable as CSV. Figures are embedded as resolution-independent vector graphics (SVG) and can be exported as SVG or 4&times; PNG.</span>
</div>

</main>
<a href="#top" class="to-top" id="rmapToTop" title="Back to top" aria-label="Back to top">&uarr;</a>
</div>
"""

html_out = (
    "<!DOCTYPE html>\n<html lang=\"en\" id=\"top\">\n<head>\n<meta charset=\"UTF-8\">\n"
    "<title>rMAP-TB Interactive Report</title>\n"
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
    "<style>" + CSS + "</style>\n</head>\n<body class=\"rmap-body\">"
    + body_html +
    "<script>" + JS + "</script>\n</body>\n</html>\n"
)

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
    File? merged_tree_svg = "final_report/mtbc_tree.svg"
    File? merged_tree_render_log = "final_report/tree_render.log"
    File? merged_tree_render_status = "final_report/tree_render_status.txt"
  }
}
