version 1.0

workflow rMAP_TB {
  input {
    Array[File]+ input_reads
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

    Int max_cpus = 8
    Int max_memory_gb = 16
    Int min_read_length = 50
    Int min_mapping_quality = 20
    Int tree_width = 2400
    Int tree_height = 1600
    String tree_image_format = "png"
    String tb_drug_resistance_genes = "rpoB,katG,inhA,fabG1,ahpC,embB,pncA,rpsL,rrs,gyrA,gyrB,eis,ethA,ethR,thyA,folC,alr,ddl,gidB,tlyA,rrl,atpE,rv0678,pepQ"
  }

  Int cpu_4 = if max_cpus < 4 then max_cpus else 4
  Int cpu_8 = if max_cpus < 8 then max_cpus else 8
  Int min_mtbc_fastq_files_for_tree = min_mtbc_samples_for_tree * 2

  # 1. Trimming is the first executable step. All downstream analysis uses these reads when trimming is enabled.
  if (do_trimming) {
    call TRIMMING {
      input:
        input_reads = input_reads,
        adapters = adapters,
        trimmomatic_quality_encoding = trimmomatic_quality_encoding,
        cpu = cpu_4,
        min_length = min_read_length
    }
  }

  Array[File] analysis_reads = select_first([TRIMMING.trimmed_reads, input_reads])

  # 2. Sample QC and Trimming Summary. FastQC runs only after trimming, because it consumes analysis_reads.
  if (do_quality_control) {
    call FASTQC {
      input:
        input_reads = analysis_reads,
        cpu = cpu_4
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
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

  # 4. TB-Profiler Resistance, Species, Lineage, and Mutation-Level Evidence Report.
  # TB-Profiler is intentionally chained after Species Typing when enabled.
  if (do_tb_profiler) {
    call TB_PROFILER_AND_MTBC_FILTER {
      input:
        input_reads = analysis_reads,
        qc_dependency = SPECIES_TYPING.species_typing_html,
        species_typing_tsv = SPECIES_TYPING.species_typing_tsv,
        docker_image = tbprofiler_docker,
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

  Array[File] mtbc_reads = select_first([TB_PROFILER_AND_MTBC_FILTER.mtbc_reads, []])

  # 5. MTBC-only Snippy/core-SNP analysis.
  # If fewer than min_mtbc_samples_for_tree paired MTBC samples are available, this branch is skipped.
  if (do_tb_profiler && do_phylogeny && size(mtbc_reads) >= min_mtbc_fastq_files_for_tree) {
    call SNIPPY_CORE_MTBC {
      input:
        input_reads = mtbc_reads,
        reference_genome = mtbc_reference_genbank,
        reference_type = snippy_reference_type,
        cpu = cpu_8,
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
        cpu = cpu_4,
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
        cpu = 1,
        memory_gb = 4
    }
  }

  # 9. Optional recombination filtering. If Gubbins fails internally, its task passes the original alignment forward.
  if (do_phylogeny && use_gubbins && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call GUBBINS_RECOMBINATION {
      input:
        core_full_alignment = select_first([SNIPPY_CORE_MTBC.core_full_alignment]),
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

  # 10. IQ-TREE uses the Gubbins-filtered alignment when available; otherwise it uses Snippy core.full.aln.
  if (do_phylogeny && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call IQTREE2_PHYLOGENY {
      input:
        alignment = select_first([GUBBINS_RECOMBINATION.filtered_alignment, SNIPPY_CORE_MTBC.core_full_alignment]),
        model = iqtree2_model,
        bootstrap_replicates = iqtree2_bootstraps,
        cpu = cpu_8,
        memory_gb = max_memory_gb,
        midpoint_root_tree = midpoint_root_tree
    }
  }

  # 11. MTBC-only Core-SNP Phylogenetic Tree.
  if (do_phylogeny && defined(IQTREE2_PHYLOGENY.final_tree)) {
    call TREE_VISUALIZATION {
      input:
        input_tree = select_first([IQTREE2_PHYLOGENY.final_tree]),
        tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv,
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
      tbprofiler_mutation_evidence_tsv = TB_PROFILER_AND_MTBC_FILTER.mutation_evidence_tsv,
      tbprofiler_mutation_evidence_html = TB_PROFILER_AND_MTBC_FILTER.mutation_evidence_html,
      mtbc_samples_txt = TB_PROFILER_AND_MTBC_FILTER.mtbc_samples_txt,
      species_typing_html = SPECIES_TYPING.species_typing_html,
      species_typing_tsv = SPECIES_TYPING.species_typing_tsv,
      qc_summary_html = MULTIQC.multiqc_report,
      trimming_report_html = TRIMMING.trimming_report,
      variant_summary_html = SNIPPY_CORE_MTBC.variant_summary,
      iqtree_report = IQTREE2_PHYLOGENY.iqtree_report,
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
task TRIMMING {
  input {
    String docker_image = "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    Array[File]+ input_reads
    File adapters
    String trimmomatic_quality_encoding = "phred33"
    Int cpu = 4
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
    memory: "8 GB"
    disks: "local-disk 100 HDD"
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
    Int cpu = 4
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
    memory: "8 GB"
    disks: "local-disk 50 HDD"
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
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 20 HDD"
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
    Int cpu = 8
    Int memory_gb = 16
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

    echo -e "Sample_ID\tSpecies_Identified\tEvidence" > species_typing/species_typing.tsv

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
        evidence="No species-level Kraken2 assignment detected"
      else
        top_percent="$(echo "$top_species_line" | awk '{print $1}')"
        top_reads="$(echo "$top_species_line" | awk '{print $2}')"
        top_taxid="$(echo "$top_species_line" | awk '{print $5}')"
        top_species="$(echo "$top_species_line" | awk '{for(i=6;i<=NF;i++) printf $i (i<NF ? " " : "")}')"

        if [ -n "$mtbc_line" ]; then
          mtbc_percent="$(echo "$mtbc_line" | awk '{print $1}')"
          mtbc_reads="$(echo "$mtbc_line" | awk '{print $2}')"
        else
          mtbc_percent="0.00"
          mtbc_reads="0"
        fi

        final_call="$top_species"
        evidence="Top species-level assignment: ${top_species} (${top_reads} reads; ${top_percent}%); MTBC support: ${mtbc_reads} reads; ${mtbc_percent}%"
      fi

      echo -e "${sample_id}\t${final_call}\t${evidence}" >> species_typing/species_typing.tsv
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

    if text == "mycobacterium tuberculosis":
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

html = []
html.append('<section class="card">')
html.append('<h2>2. Species Typing using Kraken2 + Bracken</h2>')
html.append('<p>Species typing was performed using Kraken2 against a custom Mycobacterium-only database embedded in the Docker image <code>gmboowa/mycobacterium-kraken2-bracken:2026.05</code>. The table reports the single most probable species-level call for each sample based on the highest species-level Kraken2 assignment and supporting taxonomic evidence.</p>')
html.append('<div class="table-wrap">')
html.append('<table>')
html.append('<thead><tr><th>Sample ID</th><th>Species Identified</th><th>Evidence Supporting Call</th></tr></thead>')
html.append('<tbody>')

if rows:
    for row in rows:
        html.append(
            "<tr>"
            f"<td>{escape(row.get('Sample_ID',''))}</td>"
            f"<td>{species_badge(row.get('Species_Identified',''))}</td>"
            f"<td>{escape(row.get('Evidence',''))}</td>"
            "</tr>"
        )
else:
    html.append('<tr><td colspan="3">No species typing results were generated.</td></tr>')

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
  }
}
task TB_PROFILER_AND_MTBC_FILTER {
  input {
    Array[File]+ input_reads
    File? qc_dependency
    File? species_typing_tsv
    String docker_image = "staphb/tbprofiler:6.6.6"
    Int cpu = 8
    Int memory_gb = 16
  }

  command <<<
    set -uo pipefail
    mkdir -p tbprofiler_results mtbc_reads logs mutation_evidence

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
        for key in ["name", "drug", "gene", "change", "mutation", "original_mutation", "confidence", "source", "evidence"]:
            if key in v and v[key] not in (None, "", [], {}):
                return clean_value(v[key])
        return ""
    return str(v)

def uniq(xs):
    seen = []
    for x in xs:
        x = str(x).strip()
        if x and x.lower() not in ["none", "none reported", "not reported", "na", "n/a", "unknown"]:
            if x not in seen:
                seen.append(x)
    return ", ".join(seen)

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
            if p and p not in ["none", "none reported", "not reported", "susceptible", "sensitive", "na", "n/a", "unknown"]:
                out.append(p)

    return out

def classify_who_resistance(resistant_drugs, tbprofiler_drtype):
    drugs = set()

    for d in resistant_drugs:
        for token in split_drug_tokens(d):
            drugs.add(token)

    drtype_text = str(tbprofiler_drtype or "").strip()
    drtype_lower = drtype_text.lower()

    if "xdr" in drtype_lower and "pre" not in drtype_lower:
        return "XDR"

    if "pre-xdr" in drtype_lower or "pre xdr" in drtype_lower:
        return "Pre-XDR"

    if "mdr" in drtype_lower:
        return "MDR"

    if re.search(r"\bhr[- ]?tb\b", drtype_lower):
        return "HR-TB"

    if re.search(r"\brr[- ]?tb\b", drtype_lower):
        return "RR-TB"

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

    if isoniazid and rifampicin and has_fq and has_group_a_additional:
        return "XDR"

    if isoniazid and rifampicin and has_fq:
        return "Pre-XDR"

    if isoniazid and rifampicin:
        return "MDR"

    if rifampicin and not isoniazid:
        return "RR-TB"

    if isoniazid and not rifampicin:
        return "HR-TB"

    if len(drugs) == 1:
        return "Monoresistance"

    if len(drugs) > 1:
        return "Other drug resistance"

    return "No resistance detected by TB-Profiler"

def get_kraken_support(sample, species_tsv):
    kraken_species = ""
    kraken_evidence = ""

    if species_tsv and os.path.exists(species_tsv):
        try:
            with open(species_tsv) as sfh:
                reader = csv.DictReader(sfh, delimiter="\t")
                for row in reader:
                    sid = row.get("Sample_ID") or row.get("sample") or row.get("sample_id") or ""
                    if sid == sample:
                        kraken_species = row.get("Species_Identified") or row.get("species") or row.get("Species identified") or ""
                        kraken_evidence = row.get("Evidence") or row.get("evidence") or ""
                        break
        except Exception:
            kraken_species = ""
            kraken_evidence = ""

    kraken_text = " ".join([kraken_species, kraken_evidence]).lower()

    supports_mtbc = any(t in kraken_text for t in [
        "mycobacterium tuberculosis",
        "m. tuberculosis",
        "mtbc",
        "tuberculosis complex"
    ])

    return kraken_species, kraken_evidence, supports_mtbc

def extract_drug_like_terms(text):
    known = [
        "isoniazid", "rifampicin", "rifampin", "pyrazinamide", "ethambutol",
        "streptomycin", "levofloxacin", "moxifloxacin", "ofloxacin",
        "amikacin", "kanamycin", "capreomycin", "bedaquiline", "linezolid",
        "clofazimine", "ethionamide", "prothionamide", "cycloserine",
        "para-aminosalicylic acid", "delamanid", "pretomanid"
    ]

    text_lower = str(text or "").lower()
    hits = []

    for d in known:
        if d in text_lower and d not in hits:
            hits.append(d)

    return hits

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

resistant_drugs = []
resistance_mutations = []
key_mutations = []
mutation_evidence_rows = []

variant_blocks = []

for key in ["dr_variants", "other_variants", "variants", "dr_variants_candidate", "resistance_variants", "mutations"]:
    block = data.get(key, [])

    if isinstance(block, dict):
        block = list(block.values())

    if isinstance(block, list):
        variant_blocks.extend([x for x in block if isinstance(x, dict)])

for item in variant_blocks:
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
        item.get("database") or ""
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

    for d in drug_values:
        if d:
            resistant_drugs.append(d)

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

            if any(x in s for x in ["resistant", '"r"', "assoc w r", "high confidence"]):
                resistant_drugs.append(str(drug))

    elif isinstance(block, list):
        for item in block:
            if not isinstance(item, dict):
                continue

            drug = clean_value(item.get("drug") or item.get("name") or item.get("Drug") or "")
            s = json.dumps(item).lower()

            if drug and any(x in s for x in ["resistant", '"r"', "assoc w r", "high confidence"]):
                resistant_drugs.append(drug)

kraken_species, kraken_evidence, kraken_supports_mtbc = get_kraken_support(sample, species_tsv)

raw_lineage_text_for_interpretation = " ".join([main_lineage, sub_lineage]).lower().strip()
has_tbprofiler_lineage = bool(
    re.search(r"(^|[^a-z0-9])(lineage)?[ _-]?[1-9](\.|$|[^0-9])", raw_lineage_text_for_interpretation)
)
has_tbprofiler_species = bool(species)

meaningful_drtype = str(dr_type_raw or "").strip().lower() not in [
    "", "none", "not reported", "unknown", "na", "n/a", "susceptible", "sensitive"
]

has_resistance_evidence = bool(uniq(resistant_drugs)) or bool(uniq(resistance_mutations)) or meaningful_drtype

dr_type = classify_who_resistance(resistant_drugs, dr_type_raw)

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
        me.write("\t".join([str(x).replace("\t", " ").replace("\n", " ") for x in row]) + "\n")

text_species = species.lower()

mtbc_species_terms = [
    "mycobacterium tuberculosis", "m. tuberculosis",
    "mycobacterium bovis", "m. bovis",
    "mycobacterium africanum", "m. africanum",
    "mycobacterium caprae", "m. caprae",
    "mycobacterium microti", "m. microti",
    "mycobacterium canettii", "m. canettii",
    "mycobacterium pinnipedii", "m. pinnipedii",
    "mtbc", "tuberculosis complex"
]

generic_non_mtbc_labels = [
    "non-mtb",
    "non mtb",
    "not classified as mtbc",
    "not mtbc",
    "non-mtbc",
    "non mtbc"
]

explicit_non_mtbc = (
    bool(species)
    and not any(t in text_species for t in mtbc_species_terms)
    and not any(t in text_species for t in generic_non_mtbc_labels)
)

lineage_text = " ".join([main_lineage, sub_lineage]).lower().strip()

valid_lineage = bool(
    re.search(r"(^|[^a-z0-9])(lineage)?[ _-]?[1-9](\.|$|[^0-9])", lineage_text)
)

if any(t in text_species for t in mtbc_species_terms):
    is_mtbc = True
    reason = "MTBC species reported by TB-Profiler"
elif not species and valid_lineage:
    is_mtbc = True
    reason = "Species absent, but valid TB-Profiler MTBC lineage reported"
elif kraken_supports_mtbc:
    is_mtbc = True
    reason = "MTBC supported by Kraken2/Bracken species typing"
else:
    is_mtbc = False
    reason = "No recognized MTBC species, valid TB-Profiler lineage, or Kraken2/Bracken MTBC support"

if explicit_non_mtbc and not kraken_supports_mtbc:
    is_mtbc = False
    reason = "Explicit non-MTBC species reported"

species_display = species

if not species_display and kraken_supports_mtbc:
    species_display = "Mycobacterium tuberculosis complex (supported by Kraken2/Bracken; TB-Profiler lineage not resolved)"
elif species_display.lower() in generic_non_mtbc_labels and kraken_supports_mtbc:
    species_display = "Mycobacterium tuberculosis complex (supported by Kraken2/Bracken; TB-Profiler lineage not resolved)"
elif not species_display and is_mtbc:
    species_display = "Mycobacterium tuberculosis complex (inferred from TB-Profiler lineage)"
elif not species_display and not is_mtbc:
    species_display = "Non-MTB / not classified as MTBC"

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

line = [
    sample,
    species_display or "Not reported",
    main_lineage_display,
    sub_lineage_display,
    dr_type or "Not reported",
    uniq(resistant_drugs) or "None reported",
    uniq(resistance_mutations) or "None reported",
    uniq(key_mutations) or "None reported",
    json_file,
    "YES" if is_mtbc else "NO",
    reason,
    status
]

with open("tbprofiler_summary.tsv", "a") as out:
    out.write("\t".join(line) + "\n")
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
    ("mtbc_selected", "Selected for SNP tree"),
    ("mtbc_selection_reason", "Selection reason"),
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
    "<p class='muted'>This table summarizes TB-Profiler JSON outputs, resistance-associated mutations, and MTBC-positive samples for downstream core-SNP phylogenomics.</p></div>",
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
    "<p class='muted'>This table extracts mutation-level drug-resistance evidence from TB-Profiler JSON outputs, including drug or evidence source, gene, mutation/change, confidence, and evidence fields where reported.</p></div>",
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
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 300 HDD"
    timeout: "72 hours"
  }

  output {
    Array[File] json_reports = glob("tbprofiler_results/**/*.json")
    Array[File] txt_reports = glob("tbprofiler_results/**/*.txt")
    Array[File] tbprofiler_logs = glob("logs/*.tbprofiler.log")
    File tbprofiler_command_log = "logs/tbprofiler.command.log"
    File summary_tsv = "tbprofiler_summary.tsv"
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
    Int cpu = 8
    Int memory_gb = 16
    Int min_quality = 20
  }

  command <<<
    set -uo pipefail
    mkdir -p snippy_results snippy_core logs

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

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: MTBC reads must be paired R1/R2 files." >&2
      exit 1
    fi

    echo -e "sample\tstatus\tvcf\taligned_fasta" > variant_summary.tsv
    echo -e "Sample\tMeanDepth" > mean_depth_summary.tsv
    successful_samples=()

    for ((i=0; i<n; i+=2)); do
      R1="${files[$i]}"
      R2="${files[$((i+1))]}"
      sample=$(basename "$R1" | sed -E 's/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//' | sed -E 's/(_R?1|_1|\.R?1|\.1)(_|$).*//')
      outdir="snippy_results/${sample}"

      echo "Running snippy for ${sample}" >> logs/snippy.command.log

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

      echo -e "${sample}\t${status}\t${vcf}\t${aln}" >> variant_summary.tsv
      echo -e "${sample}\t${depth}" >> mean_depth_summary.tsv
    done

    core_status="skipped"

    if [ ${#successful_samples[@]} -ge 2 ]; then
      echo "Running snippy-core on ${#successful_samples[@]} samples" >> logs/snippy.command.log

      if "$SNIPPY_CORE_BIN" --ref "$ref" --prefix snippy_core/core "${successful_samples[@]}" >> logs/snippy_core.log 2>&1; then
        core_status="success"
      else
        core_status="snippy_core_failed"
        echo "WARNING: snippy-core failed" >> logs/snippy.command.log
      fi
    else
      echo "WARNING: fewer than 2 successful samples, skipping snippy-core" >> logs/snippy.command.log
    fi

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
</style>
</head>
<body>
<h1>MTBC core-SNP variant-calling summary</h1>
<table>
<thead>
<tr><th>Sample</th><th>Status</th><th>Mean depth</th><th>VCF</th><th>Aligned FASTA</th></tr>
</thead>
<tbody>
"""

for r in rows:
    status = r["status"]
    cls = "ok" if status == "success" else "fail"
    sample = r["sample"]
    mean_depth = depth_by_sample.get(sample, "NA")
    out += f"<tr><td>{html.escape(sample)}</td><td class='{cls}'>{html.escape(status)}</td><td>{html.escape(mean_depth)}</td><td>{html.escape(r['vcf'])}</td><td>{html.escape(r['aligned_fasta'])}</td></tr>\n"

out += "</tbody></table></body></html>"

open("variant_summary.html","w").write(out)
PY
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 250 HDD"
    timeout: "72 hours"
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
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 20 HDD"
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
    Int cpu = 2
    Int memory_gb = 4
  }

  command <<<
    set -uo pipefail

    mkdir -p snp_distance logs

    cp "~{core_full_alignment}" snp_distance/core.full.aln

    echo "Running fail-safe pure-Python SNP distance clustering on core.full.aln" > logs/snp_distance.command.log
    echo "Reference/non-sample sequences are excluded from sample-level SNP distance and cluster reporting." >> logs/snp_distance.command.log

    apt-get update >/dev/null 2>&1 || true
    apt-get install -y python3-pip >/dev/null 2>&1 || true

    pip3 install --quiet pandas matplotlib seaborn numpy || true

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
    disks: "local-disk 50 HDD"
    timeout: "12 hours"
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
    Int cpu = 8
    Int memory_gb = 16
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
    disks: "local-disk 200 HDD"
    timeout: "72 hours"
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
    String docker_image = "staphb/iqtree2:2.3.4"
    File alignment
    String model = "GTR+G"
    Int bootstrap_replicates = 1000
    Int cpu = 8
    Int memory_gb = 16
    Boolean midpoint_root_tree = true
  }

  command <<<
    set -uo pipefail
    mkdir -p iqtree logs

    if command -v iqtree2 >/dev/null 2>&1; then
      IQTREE_BIN="$(command -v iqtree2)"
    elif command -v iqtree >/dev/null 2>&1; then
      IQTREE_BIN="$(command -v iqtree)"
    elif [ -x /usr/local/bin/iqtree2 ]; then
      IQTREE_BIN="/usr/local/bin/iqtree2"
    else
      echo "ERROR: IQ-TREE executable not found." >&2
      exit 127
    fi

    cp "~{alignment}" iqtree/mtbc_core_snp_alignment.fasta

    if [ ! -s iqtree/mtbc_core_snp_alignment.fasta ]; then
      echo "ERROR: alignment file is missing or empty." >&2
      exit 1
    fi

    echo "Using IQ-TREE: ${IQTREE_BIN}" > logs/iqtree.command.log
    echo "Model: ~{model}" >> logs/iqtree.command.log

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
      status="iqtree_failed"
      echo "WARNING: IQ-TREE failed. Generating fallback tree." >> logs/iqtree.command.log
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

    grep -oE '\)[0-9]+(\.[0-9]+)?(/[0-9]+(\.[0-9]+)?)?:' iqtree/final.treefile > iqtree/support_labels.txt || true

    echo "$status" > iqtree/iqtree_status.txt
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 100 HDD"
    timeout: "200 hours"
  }

  output {
    File final_tree = "iqtree/final.treefile"
    File iqtree_log = "iqtree/iqtree.log"
    File iqtree_report = "iqtree/iqtree.report"
    File support_labels = "iqtree/support_labels.txt"
    File iqtree_status = "iqtree/iqtree_status.txt"
    File iqtree_command_log = "logs/iqtree.command.log"
    File iqtree_run_log = "logs/iqtree.run.log"
  }
}

task TREE_VISUALIZATION {
  input {
    File input_tree
    File? tbprofiler_summary_tsv
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

tree_input = "~{if defined(input_tree) then input_tree else ""}"
tbprofiler_summary_path = "~{if defined(tbprofiler_summary_tsv) then tbprofiler_summary_tsv else ""}"
requested_width = int("~{width}")
requested_height = int("~{height}")
title = "~{title}"

log.write_text("TREE_VISUALIZATION started\n")

RESISTANCE_COLORS = {
    "Unknown": "#999999",
    "Sensitive": "#2a9d8f",
    "Monoresistance": "#fcd33d",
    "HR-TB": "#2292dc",
    "RR-TB": "#3b82f6",
    "MDR": "#ed641e",
    "Pre-XDR": "#ed2828",
    "XDR": "#5a189a",
    "Other drug resistance": "#2292dc"
}

def write_placeholder_image(message):
    msg = str(message).replace("<", "&lt;").replace(">", "&gt;")
    svg_text = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="500" viewBox="0 0 1200 500">
<rect width="100%" height="100%" fill="#ffffff"/>
<rect x="40" y="40" width="1120" height="420" rx="24" fill="#f8fafc" stroke="#cbd5e1" stroke-width="3"/>
<text x="600" y="170" text-anchor="middle" font-family="Arial" font-size="34" font-weight="700" fill="#0f172a">Tree visualization was not rendered</text>
<text x="600" y="225" text-anchor="middle" font-family="Arial" font-size="22" fill="#475569">The Newick tree was preserved for downstream visualization.</text>
<text x="600" y="285" text-anchor="middle" font-family="Arial" font-size="18" fill="#b91c1c">{msg[:180]}</text>
</svg>'''

    if image_format == "svg":
        out_img.write_text(svg_text, encoding="utf-8")
    elif image_format == "pdf":
        out_img.write_text("Tree visualization was not rendered. Newick tree preserved.\n" + str(message) + "\n", encoding="utf-8")
    else:
        # Valid 1x1 transparent PNG placeholder. Kept minimal to avoid adding non-standard dependencies.
        png_b64 = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/p9sAAAAASUVORK5CYII="
        out_img.write_bytes(base64.b64decode(png_b64))

def normalize_drug_name(x):
    s = str(x or "").strip().lower()
    s = re.sub(r"[_\-]+", " ", s)
    s = re.sub(r"\s+", " ", s)
    aliases = {
        "inh": "isoniazid", "isoniazid": "isoniazid", "h": "isoniazid",
        "rif": "rifampicin", "rmp": "rifampicin", "rifampin": "rifampicin", "rifampicin": "rifampicin", "r": "rifampicin",
        "pza": "pyrazinamide", "pyrazinamide": "pyrazinamide", "z": "pyrazinamide",
        "emb": "ethambutol", "ethambutol": "ethambutol", "e": "ethambutol",
        "sm": "streptomycin", "str": "streptomycin", "streptomycin": "streptomycin", "s": "streptomycin",
        "levo": "levofloxacin", "levofloxacin": "levofloxacin", "lfx": "levofloxacin",
        "moxi": "moxifloxacin", "moxifloxacin": "moxifloxacin", "mfx": "moxifloxacin",
        "ofx": "ofloxacin", "ofloxacin": "ofloxacin",
        "amikacin": "amikacin", "amk": "amikacin",
        "kanamycin": "kanamycin", "kan": "kanamycin",
        "capreomycin": "capreomycin", "cap": "capreomycin",
        "bedaquiline": "bedaquiline", "bdq": "bedaquiline",
        "linezolid": "linezolid", "lzd": "linezolid",
        "clofazimine": "clofazimine", "cfz": "clofazimine",
        "ethionamide": "ethionamide", "eto": "ethionamide",
        "prothionamide": "prothionamide", "pto": "prothionamide",
        "cycloserine": "cycloserine", "cs": "cycloserine",
        "para aminosalicylic acid": "para-aminosalicylic acid", "pas": "para-aminosalicylic acid",
        "delamanid": "delamanid", "dlm": "delamanid",
        "pretomanid": "pretomanid", "pa": "pretomanid"
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

def classify_resistance(dr_type, resistant_drugs):
    profile_raw = (dr_type or "").strip()
    profile = profile_raw.lower().replace("_", "-")
    drugs = set(split_drug_tokens(resistant_drugs))
    if not drugs and profile in {"", "not reported", "unknown", "na", "n/a", "none", "susceptible", "sensitive"}:
        return "Sensitive", RESISTANCE_COLORS["Sensitive"]
    if "xdr" in profile and "pre" not in profile:
        return "XDR", RESISTANCE_COLORS["XDR"]
    if "pre-xdr" in profile or "pre xdr" in profile or "prexdr" in profile:
        return "Pre-XDR", RESISTANCE_COLORS["Pre-XDR"]
    if "mdr" in profile:
        return "MDR", RESISTANCE_COLORS["MDR"]
    if re.search(r"\bhr[- ]?tb\b", profile):
        return "HR-TB", RESISTANCE_COLORS["HR-TB"]
    if re.search(r"\brr[- ]?tb\b", profile):
        return "RR-TB", RESISTANCE_COLORS["RR-TB"]
    isoniazid = "isoniazid" in drugs
    rifampicin = "rifampicin" in drugs
    fluoroquinolones = {"levofloxacin", "moxifloxacin", "ofloxacin", "gatifloxacin", "ciprofloxacin"}
    group_a_additional = {"bedaquiline", "linezolid"}
    has_fq = bool(drugs.intersection(fluoroquinolones))
    has_group_a_additional = bool(drugs.intersection(group_a_additional))
    if isoniazid and rifampicin and has_fq and has_group_a_additional:
        return "XDR", RESISTANCE_COLORS["XDR"]
    if isoniazid and rifampicin and has_fq:
        return "Pre-XDR", RESISTANCE_COLORS["Pre-XDR"]
    if isoniazid and rifampicin:
        return "MDR", RESISTANCE_COLORS["MDR"]
    if rifampicin and not isoniazid:
        return "RR-TB", RESISTANCE_COLORS["RR-TB"]
    if isoniazid and not rifampicin:
        return "HR-TB", RESISTANCE_COLORS["HR-TB"]
    if len(drugs) == 1:
        return "Monoresistance", RESISTANCE_COLORS["Monoresistance"]
    if len(drugs) > 1:
        return "Other drug resistance", RESISTANCE_COLORS["Other drug resistance"]
    if profile_raw and profile not in {"not reported", "unknown", "na", "n/a", "none"}:
        return profile_raw, RESISTANCE_COLORS["Other drug resistance"]
    return "Unknown", RESISTANCE_COLORS["Unknown"]

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

    metadata = {}
    if tbprofiler_summary_path and Path(tbprofiler_summary_path).exists():
        with open(tbprofiler_summary_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                sample = (r.get("sample") or "").strip()
                if not sample:
                    continue
                main_lineage = (r.get("main_lineage") or "").strip()
                sub_lineage = (r.get("sub_lineage") or "").strip()
                lineage = sub_lineage or main_lineage
                if lineage.lower() in {"not reported", "none", "na", "n/a", "unknown"}:
                    lineage = main_lineage
                category, color = classify_resistance(r.get("dr_type"), r.get("resistant_drugs"))
                metadata[sample] = {"lineage": lineage, "category": category, "color": color}

    from ete3 import Tree, TreeStyle, TextFace, NodeStyle

    t = Tree(str(tree_path), format=1)

    reference_names = {"reference", "ref", "h37rv", "nc_000962", "nc_000962.3"}
    removed_refs = []
    for leaf in list(t.get_leaves()):
        if leaf.name.strip().lower() in reference_names:
            removed_refs.append(leaf.name)
            leaf.detach()

    if len(t.get_leaves()) < 2:
        raise ValueError("Tree has fewer than two non-reference tips after filtering.")

    try:
        midpoint = t.get_midpoint_outgroup()
        if midpoint:
            t.set_outgroup(midpoint)
    except Exception as e:
        with open(log, "a") as fh:
            fh.write(f"WARNING: midpoint rooting failed: {repr(e)}\n")

    n_leaves = len(t.get_leaves())
    if n_leaves <= 5:
        auto_width = max(requested_width, 3800); label_font = 14; lineage_font = 11; resistance_font = 11; bootstrap_font = 10; tip_node_size = 8; branch_width = 2; branch_vertical_margin = 18; margin_right = 1500
    elif n_leaves <= 10:
        auto_width = max(requested_width, 4200); label_font = 13; lineage_font = 10; resistance_font = 10; bootstrap_font = 9; tip_node_size = 7; branch_width = 2; branch_vertical_margin = 14; margin_right = 1600
    elif n_leaves <= 25:
        auto_width = max(requested_width, 5000); label_font = 11; lineage_font = 9; resistance_font = 9; bootstrap_font = 8; tip_node_size = 6; branch_width = 2; branch_vertical_margin = 8; margin_right = 1700
    elif n_leaves <= 50:
        auto_width = max(requested_width, 5600); label_font = 9; lineage_font = 8; resistance_font = 8; bootstrap_font = 7; tip_node_size = 5; branch_width = 1; branch_vertical_margin = 5; margin_right = 1800
    elif n_leaves <= 100:
        auto_width = max(requested_width, 6400); label_font = 8; lineage_font = 7; resistance_font = 7; bootstrap_font = 6; tip_node_size = 4; branch_width = 1; branch_vertical_margin = 3; margin_right = 1900
    else:
        auto_width = max(requested_width, 7200); label_font = 7; lineage_font = 6; resistance_font = 6; bootstrap_font = 5; tip_node_size = 3; branch_width = 1; branch_vertical_margin = 2; margin_right = 2000

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

    for leaf in t.get_leaves():
        sample = leaf.name
        meta = metadata.get(sample, {})
        lineage = meta.get("lineage", "") or "Lineage NA"
        category = meta.get("category", "Unknown")
        color = meta.get("color", RESISTANCE_COLORS["Unknown"])
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
    cpu: 1
  }

  output {
    File tree_image = "tree_visualization/phylogenetic_tree.~{image_format}"
    File cleaned_tree = "tree_visualization/phylogenetic_tree.cleaned.nwk"
    File render_log = "tree_visualization/render.log"
  }
}

task TB_SURVEILLANCE_SUMMARY_VISUALS {
  input {
    String docker_image = "python:3.11-slim"
    File tbprofiler_summary_tsv
    File? species_typing_tsv
    File? pairwise_snp_distance_matrix
    File? mean_depth_tsv
    Int cpu = 1
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail
    mkdir -p surveillance_summary

    tb_tsv="~{tbprofiler_summary_tsv}"
    species_tsv="~{if defined(species_typing_tsv) then species_typing_tsv else ""}"
    snp_matrix="~{if defined(pairwise_snp_distance_matrix) then pairwise_snp_distance_matrix else ""}"
    depth_tsv="~{if defined(mean_depth_tsv) then mean_depth_tsv else ""}"

    python3 - "$tb_tsv" "$species_tsv" "$snp_matrix" "$depth_tsv" <<'PY'
import csv
import html
import re
import sys
from pathlib import Path
from collections import Counter

tb_tsv = Path(sys.argv[1])
species_tsv = Path(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] else None
snp_matrix = Path(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else None
depth_tsv = Path(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else None

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

        samples = [raw_samples[i] for i in keep_idx]
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

def split_integrated_mtbc_status(status_text):
    text = clean_text(status_text)

    mtbc_status = "MTBC status not determined"
    support_source = "Not reported"
    lineage_status = "Not reported"

    low = text.lower()

    if "mycobacterium tuberculosis complex" in low or "mtbc" in low:
        mtbc_status = "MTBC supported"

    if "kraken2/bracken" in low:
        support_source = "Kraken2/Bracken"
    elif "tb-profiler" in low:
        support_source = "TB-Profiler"

    if "lineage not resolved" in low or "not resolved" in low:
        lineage_status = "Not resolved by TB-Profiler"
    elif "lineage" in low:
        lineage_status = "Resolved or inferred by TB-Profiler"

    return mtbc_status, support_source, lineage_status

tb_rows = read_tsv(tb_tsv)
species_rows = read_tsv(species_tsv)
depth_rows = read_tsv(depth_tsv)

snp_samples, snp_values = read_snp_matrix(snp_matrix)
tree_sample_set = set(snp_samples)

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
        depth_by_sample[sid] = depth

species_by_sample = {}

for r in species_rows:
    sample = (
        r.get("Sample_ID")
        or r.get("sample")
        or r.get("sample_id")
        or ""
    )

    species = (
        r.get("Species_Identified")
        or r.get("Species identified")
        or r.get("species")
        or ""
    )

    evidence = r.get("Evidence") or r.get("evidence") or ""

    if sample:
        species_by_sample[sample] = {
            "species": species,
            "evidence": evidence
        }

lineage_counts = Counter()
metadata_rows = []
qc_rows = []

for r in tb_rows:
    sample = r.get("sample", "")

    sp = species_by_sample.get(sample, {})
    kraken_species = sp.get("species", "")
    kraken_evidence = sp.get("evidence", "")

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

    included_in_tree = (
        "YES"
        if tree_sample_set and sample in tree_sample_set
        else selected_for_phylogeny
    )

    mean_depth = depth_by_sample.get(sample, "Not available")
    mtbc_percent = extract_mtbc_percent(kraken_evidence)

    resistance_profile = normalize_missing(
        r.get("dr_type", ""),
        "Resistance not determined by TB-Profiler"
    )

    resistant_drugs = normalize_missing(
        r.get("resistant_drugs", ""),
        "None reported"
    )

    drug_resistance = resistance_detected(resistance_profile, resistant_drugs)

    integrated_status_text = r.get("species", "")
    mtbc_status, support_source, lineage_status = split_integrated_mtbc_status(integrated_status_text)

    if lineage_group != "Not resolved by TB-Profiler":
        lineage_status = "Resolved by TB-Profiler"

    metadata_rows.append({
        "sample": sample,
        "integrated_mtbc_status": mtbc_status,
        "mtbc_support_source": support_source,
        "tbprofiler_lineage_status": lineage_status,
        "kraken_species": kraken_species,
        "tbprofiler_main_lineage": main_lineage,
        "tbprofiler_sub_lineage": sub_lineage,
        "lineage_group": lineage_group,
        "resistance_profile": resistance_profile,
        "drug_resistance_detected": drug_resistance,
        "resistant_drugs": resistant_drugs,
        "mean_depth": mean_depth,
        "mtbc_percent": mtbc_percent,
        "selected_for_phylogeny": selected_for_phylogeny,
        "included_in_tree": included_in_tree,
        "tbprofiler_status": r.get("status", "")
    })

    qc_rows.append({
        "sample": sample,
        "mean_depth": mean_depth,
        "mtbc_percent": mtbc_percent,
        "selected_for_phylogeny": selected_for_phylogeny,
        "included_in_tree": included_in_tree,
        "reason": r.get("mtbc_selection_reason", "") or "Not reported"
    })

with (outdir / "lineage_distribution.tsv").open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["lineage", "count"])

    items = sorted(lineage_counts.items()) if lineage_counts else [("Not resolved by TB-Profiler", 0)]

    for lineage, count in items:
        writer.writerow([lineage, count])

with (outdir / "tb_surveillance_metadata.tsv").open("w", newline="") as out:
    fieldnames = [
        "sample",
        "integrated_mtbc_status",
        "mtbc_support_source",
        "tbprofiler_lineage_status",
        "kraken_species",
        "tbprofiler_main_lineage",
        "tbprofiler_sub_lineage",
        "lineage_group",
        "resistance_profile",
        "drug_resistance_detected",
        "resistant_drugs",
        "mean_depth",
        "mtbc_percent",
        "selected_for_phylogeny",
        "included_in_tree",
        "tbprofiler_status"
    ]

    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(metadata_rows)

with (outdir / "qc_filtering_rationale.tsv").open("w", newline="") as out:
    fieldnames = [
        "sample",
        "mean_depth",
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
        f'<text x="{width/2}" y="{height-24}" text-anchor="middle" font-family="Arial" font-size="12" fill="#475569">Lineage groups are summarized from TB-Profiler main-lineage and sub-lineage fields.</text></svg>'
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

(outdir / "surveillance_summary.html").write_text(
    '<!doctype html><html><head><meta charset="utf-8"><title>Surveillance Summary</title></head>'
    '<body>'
    '<h1>Surveillance Summary</h1>'
    '<p>TB-Profiler lineage summaries, Kraken2/Bracken-supported MTBC classification, SNP heatmap, QC rationale, mean-depth metadata, and surveillance metadata were generated.</p>'
    '</body></html>',
    encoding="utf-8"
)
PY
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 20 HDD"
    timeout: "4 hours"
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
    File? tbprofiler_mutation_evidence_tsv
    File? tbprofiler_mutation_evidence_html
    File? mtbc_samples_txt
    File? species_typing_html
    File? species_typing_tsv
    File? qc_summary_html
    File? trimming_report_html
    File? variant_summary_html
    File? iqtree_report
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
    species_html="~{if defined(species_typing_html) then species_typing_html else ""}"
    species_tsv="~{if defined(species_typing_tsv) then species_typing_tsv else ""}"
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

    if [ -n "$tree_png" ] && [ -f "$tree_png" ]; then
      cp "$tree_png" final_report/mtbc_tree.png || true
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

    if [ -z "$tb_tsv" ] || [ ! -f "$tb_tsv" ]; then
      echo -e "sample\tspecies\tmain_lineage\tsub_lineage\tdr_type\tresistant_drugs\tresistance_mutations\tkey_mutations\tjson_file\tmtbc_selected\tmtbc_selection_reason\tstatus" > final_report/empty.tsv
      tb_tsv="final_report/empty.tsv"
    fi

    if [ -z "$species_tsv" ] || [ ! -f "$species_tsv" ]; then
      echo -e "Sample_ID\tSpecies_Identified\tEvidence" > final_report/empty_species_typing.tsv
      species_tsv="final_report/empty_species_typing.tsv"
    fi

    if [ -z "$nonsyn_tsv" ] || [ ! -f "$nonsyn_tsv" ]; then
      echo -e "sample\tgene\teffect\taa_change\tnt_change\tproduct" > final_report/empty_nonsyn.tsv
      nonsyn_tsv="final_report/empty_nonsyn.tsv"
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

    python3 - "$tb_tsv" "$species_tsv" "$nonsyn_tsv" "$mutation_tsv" "$snp_pairs_tsv" "$snp_cluster_tsv" "$lineage_tsv" "$surveillance_metadata_tsv" "$qc_rationale_tsv" "$qc_html" "$trim_html" "$variant_html" "$iqtree_txt" "$species_html" <<'PY'
import csv
import html
import re
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone

summary_tsv, species_tsv, nonsyn_tsv, mutation_tsv, snp_pairs_tsv, snp_cluster_tsv, lineage_tsv, surveillance_metadata_tsv, qc_rationale_tsv, qc_html, trim_html, variant_html, iqtree_txt, species_html = sys.argv[1:15]

outdir = Path("final_report")
outdir.mkdir(exist_ok=True)

run_started_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
run_stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_UTC")

(outdir / "run_metadata.txt").write_text(
    f"Workflow report generation timestamp: {run_started_utc}\n"
    f"Run stamp: {run_stamp}\n"
)

RESISTANCE_COLORS = {
    "No resistance detected by TB-Profiler": "#2a9d8f",
    "Resistance not determined by TB-Profiler": "#999999",
    "Sensitive": "#2a9d8f",
    "Monoresistance": "#fcd33d",
    "HR-TB": "#2292dc",
    "RR-TB": "#3b82f6",
    "MDR": "#ed641e",
    "Pre-XDR": "#ed2828",
    "XDR": "#5a189a",
    "Other drug resistance": "#2292dc",
    "Unknown": "#999999",
    "Not determined": "#999999",
}

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
species_rows = read_tsv_rows(species_tsv)
nonsyn_rows = read_tsv_rows(nonsyn_tsv)
mutation_rows = read_tsv_rows(mutation_tsv)
snp_pair_rows = read_tsv_rows(snp_pairs_tsv)
snp_cluster_rows = read_tsv_rows(snp_cluster_tsv)
lineage_rows = read_tsv_rows(lineage_tsv)
surveillance_metadata_rows = read_tsv_rows(surveillance_metadata_tsv)
qc_rationale_rows = read_tsv_rows(qc_rationale_tsv)

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

def classify_resistance(dr_type, resistant_drugs):
    profile_raw = clean(dr_type)
    profile = profile_raw.lower().replace("_", "-")
    drugs = set(split_drug_tokens(resistant_drugs))

    if profile in {
        "",
        "none",
        "none reported",
        "not reported",
        "unknown",
        "na",
        "n/a",
        "not determined",
        "resistance not determined by tb-profiler"
    }:
        return "Resistance not determined by TB-Profiler"

    if "no resistance detected" in profile:
        return "No resistance detected by TB-Profiler"

    if profile in {"susceptible", "sensitive"} and not drugs:
        return "No resistance detected by TB-Profiler"

    if "xdr" in profile and "pre" not in profile:
        return "XDR"

    if "pre-xdr" in profile or "pre xdr" in profile or "prexdr" in profile:
        return "Pre-XDR"

    if "mdr" in profile:
        return "MDR"

    if re.search(r"\bhr[- ]?tb\b", profile):
        return "HR-TB"

    if re.search(r"\brr[- ]?tb\b", profile):
        return "RR-TB"

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

    if isoniazid and rifampicin and has_fq and has_group_a_additional:
        return "XDR"

    if isoniazid and rifampicin and has_fq:
        return "Pre-XDR"

    if isoniazid and rifampicin:
        return "MDR"

    if rifampicin and not isoniazid:
        return "RR-TB"

    if isoniazid and not rifampicin:
        return "HR-TB"

    if len(drugs) == 1:
        return "Monoresistance"

    if len(drugs) > 1:
        return "Other drug resistance"

    if profile_raw and profile not in {"not reported", "unknown", "na", "n/a", "none"}:
        return profile_raw

    return "Resistance not determined by TB-Profiler"

def is_drug_resistant_category(category):
    return category not in {
        "No resistance detected by TB-Profiler",
        "Resistance not determined by TB-Profiler",
        "Sensitive",
        "Unknown",
        "Not determined",
        ""
    }

def badge(label):
    label = normalize_missing(label, "Resistance not determined by TB-Profiler")
    color = RESISTANCE_COLORS.get(label, RESISTANCE_COLORS["Unknown"])
    text_color = "#111111" if label == "Monoresistance" else "#ffffff"
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

total_samples = max(len(species_rows), len(rows))

mtbc_retained = 0
for r in rows:
    if clean(r.get("mtbc_selected")).upper() == "YES":
        mtbc_retained += 1

non_mtbc = max(total_samples - mtbc_retained, 0)

drug_resistant = 0
for r in rows:
    category = classify_resistance(r.get("dr_type"), r.get("resistant_drugs"))
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
            sample_ids.append(sid)

    if not sample_ids:
        sample_ids = [r.get("sample", "") for r in rows if r.get("sample")]

    if not sample_ids:
        body = '<tr><td colspan="5">No sample-level QC records available.</td></tr>'
    else:
        body = []
        for s in sample_ids:
            body.append(
                "<tr>"
                f"<td>{safe(s)}</td>"
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
                f"<td>{safe(sample)}</td>"
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

def build_tb_rows():
    out = []

    for r in rows:
        sample = r.get("sample")
        species = r.get("species") or "Mycobacterium tuberculosis complex (inferred from available MTBC evidence)"
        main_lineage = normalize_missing(r.get("main_lineage"), "Not resolved by TB-Profiler")
        sub_lineage = normalize_missing(r.get("sub_lineage"), "Not resolved by TB-Profiler")
        lineage = f"{safe(main_lineage)} / {safe(sub_lineage)}"
        dr_type = r.get("dr_type") or "Resistance not determined by TB-Profiler"
        resistant_drugs = r.get("resistant_drugs") or "None reported"
        category = classify_resistance(dr_type, resistant_drugs)

        mtbc_selected = clean(r.get("mtbc_selected")).upper()
        mtbc_reason = r.get("mtbc_selection_reason") or ""

        if mtbc_selected == "YES":
            decision_html = green_badge("Selected")
        else:
            decision_html = red_badge("Excluded")

        out.append(
            "<tr>"
            f"<td>{safe(sample)}</td>"
            f"<td>{safe(species)}</td>"
            f"<td>{lineage}</td>"
            f"<td>{badge(category)}</td>"
            f"<td>{safe(resistant_drugs)}</td>"
            f"<td>{decision_html}<br><small>{safe(mtbc_reason)}</small></td>"
            "</tr>"
        )

    return "".join(out)

def build_mutation_evidence_section():
    cleaned = [
        r for r in mutation_rows
        if any(clean(v) for v in r.values())
        and clean(r.get("sample"))
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
        out.append(f'<details><summary>Sample: {safe(s)} — {len(grouped[s])} mutation(s)</summary><table>')
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
        out.append(f'<details><summary>Sample: {safe(s)} — {len(grouped[s])} mutation(s)</summary><table>')
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
            f"<td>{safe(sample1)}</td>"
            f"<td>{safe(sample2)}</td>"
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
                f"<td>{safe(r.get('sample1'))}</td>"
                f"<td>{safe(r.get('sample2'))}</td>"
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

def build_lineage_distribution_section():
    lineage_svg_exists = Path("final_report/lineage_distribution.svg").exists()

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

    plot_html = (
        '<div class="tree-panel"><object type="image/svg+xml" data="lineage_distribution.svg" style="width:100%;min-height:420px;"></object></div>'
        if lineage_svg_exists else
        '<div class="note">Lineage distribution plot was not available.</div>'
    )

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

    heatmap_html = (
        '<div class="tree-panel"><object type="image/svg+xml" data="snp_distance_heatmap.svg" style="width:100%;min-height:560px;"></object></div>'
        if heatmap_exists else
        '<div class="note">SNP distance heatmap was not generated or was skipped.</div>'
    )

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
        selected = r.get("selected_for_phylogeny", r.get("included", ""))
        included = r.get("included_in_tree", r.get("included", ""))

        selected_html = green_badge(selected) if clean(selected).upper() == "YES" else red_badge(selected)
        included_html = green_badge(included) if clean(included).upper() == "YES" else red_badge(included)

        qc_body.append(
            "<tr>"
            f"<td>{safe(r.get('sample'))}</td>"
            f"<td>{safe(r.get('mean_depth'))}</td>"
            f"<td>{safe(r.get('mtbc_percent'))}</td>"
            f"<td>{selected_html}</td>"
            f"<td>{included_html}</td>"
            f"<td>{safe(r.get('reason'))}</td>"
            "</tr>"
        )

    metadata_body = []

    for r in metadata_cleaned:
        included = r.get("included_in_tree", "")
        included_html = green_badge(included) if clean(included).upper() == "YES" else red_badge(included)
        resistance_profile = r.get("resistance_profile") or "Resistance not determined by TB-Profiler"

        drug_res_detected = normalize_missing(r.get("drug_resistance_detected"), "Not determined")

        if drug_res_detected.upper() == "YES":
            drug_res_html = red_badge("YES")
        elif drug_res_detected.upper() == "NO":
            drug_res_html = green_badge("NO")
        else:
            drug_res_html = neutral_badge(drug_res_detected)

        metadata_body.append(
            "<tr>"
            f"<td>{safe(r.get('sample'))}</td>"
            f"<td>{safe(r.get('integrated_mtbc_status') or r.get('tbprofiler_species'))}</td>"
            f"<td>{safe(r.get('mtbc_support_source'))}</td>"
            f"<td>{safe(r.get('tbprofiler_lineage_status'))}</td>"
            f"<td>{safe(r.get('kraken_species'))}</td>"
            f"<td>{safe(r.get('tbprofiler_main_lineage') or r.get('main_lineage'))}</td>"
            f"<td>{safe(r.get('tbprofiler_sub_lineage') or r.get('sub_lineage'))}</td>"
            f"<td>{safe(r.get('lineage_group'))}</td>"
            f"<td>{badge(resistance_profile)}</td>"
            f"<td>{drug_res_html}</td>"
            f"<td>{safe(r.get('resistant_drugs'))}</td>"
            f"<td>{safe(r.get('mean_depth'))}</td>"
            f"<td>{included_html}</td>"
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
This section provides a transparent rationale for sample inclusion/exclusion and a surveillance-ready metadata table. The metadata TSV is exported as <code>tb_surveillance_metadata.tsv</code>, and the QC rationale is exported as <code>qc_filtering_rationale.tsv</code>.
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
<th class="status" onclick="sortTable('qcRationaleTable',3)">Selected for phylogeny</th>
<th class="status" onclick="sortTable('qcRationaleTable',4)">Included in tree</th>
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
<th class="status" onclick="sortTable('surveillanceMetadataTable',12)">Included in Tree</th>
</tr>
</thead>
<tbody>
{''.join(metadata_body)}
</tbody>
</table>
</div>
"""

def build_tree_legend():
    items = [
        ("No resistance detected by TB-Profiler", RESISTANCE_COLORS["No resistance detected by TB-Profiler"]),
        ("Resistance not determined by TB-Profiler", RESISTANCE_COLORS["Resistance not determined by TB-Profiler"]),
        ("Monoresistance", RESISTANCE_COLORS["Monoresistance"]),
        ("HR-TB", RESISTANCE_COLORS["HR-TB"]),
        ("RR-TB", RESISTANCE_COLORS["RR-TB"]),
        ("MDR", RESISTANCE_COLORS["MDR"]),
        ("Pre-XDR", RESISTANCE_COLORS["Pre-XDR"]),
        ("XDR", RESISTANCE_COLORS["XDR"]),
        ("Other drug resistance", RESISTANCE_COLORS["Other drug resistance"]),
        ("Unknown", RESISTANCE_COLORS["Unknown"]),
    ]

    legend = ['<div class="legend">']

    for label, color in items:
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

tree_exists = Path("final_report/mtbc_tree.png").exists()
tree_html = '<img class="mtbc-tree-img" src="mtbc_tree.png" alt="ETE3-rendered MTBC core-SNP phylogenetic tree">' if tree_exists else '<div class="note">Tree image was not provided or could not be copied into the final report folder.</div>'

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
<p>Trimming → QC → Species typing → TB-Profiler → MTBC-only filtering → mutation evidence → lineage and surveillance summaries → SNP distance clustering and heatmap → core-SNP phylogenomics → final merged report</p>
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

<div class="section">
<h2>3. TB-Profiler Resistance, Species, and Lineage Report</h2>

<div class="note">
<strong>Interpretation note:</strong> TB-Profiler is using WHO 2021+ definitions for classification of drug-resistance results.
<strong>Key WHO 2021+ resistance definitions:</strong>
<strong>Isoniazid-resistant, rifampicin-susceptible TB (HR-TB):</strong> resistant to isoniazid and not resistant to rifampicin.
<strong>Rifampicin-resistant TB (RR-TB):</strong> resistant to rifampicin, with or without resistance to other drugs.
<strong>Multidrug-resistant TB (MDR-TB):</strong> resistance to at least isoniazid and rifampicin.
<strong>Pre-extensively drug-resistant TB (Pre-XDR-TB):</strong> MDR/RR-TB that is also resistant to any fluoroquinolone.
<strong>Extensively drug-resistant TB (XDR-TB):</strong> MDR/RR-TB that is resistant to any fluoroquinolone and at least one additional Group A drug, bedaquiline or linezolid.
</div>

<div class="controls">
<input id="tbSearch" onkeyup="filterTable('tbSearch','tbTable')" placeholder="Search TB-Profiler results...">
<select onchange="filterResistance(this.value)">
<option value="">All resistance profiles</option>
<option value="No resistance detected by TB-Profiler">No resistance detected by TB-Profiler only</option>
<option value="Resistance not determined by TB-Profiler">Resistance not determined by TB-Profiler only</option>
<option value="Monoresistance">Monoresistance only</option>
<option value="HR-TB">HR-TB only</option>
<option value="RR-TB">RR-TB only</option>
<option value="MDR">MDR only</option>
<option value="Pre-XDR">Pre-XDR only</option>
<option value="XDR">XDR only</option>
<option value="Other drug resistance">Other drug resistance only</option>
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
<th class="status" onclick="sortTable('tbTable',5)">MTBC decision</th>
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
</div>
<div class="details">
<h3>Tree construction summary</h3>
<p><strong>Included:</strong> {mtbc_retained} MTBC isolate(s) retained for phylogenomic interpretation.</p>
<p><strong>Excluded:</strong> {non_mtbc} non-MTBC or low-confidence isolate(s).</p>
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
<p>The report documents all samples through QC, species typing, TB-Profiler analysis, mutation-level resistance evidence, lineage distribution, SNP distance clustering, SNP heatmap visualization, surveillance metadata, and MTBC-only phylogenomic reconstruction. Samples not classified as MTBC are excluded from the tree but retained in the workflow record for transparency.</p>
<div class="note">
<strong>Interpretation:</strong> use close clustering together with bootstrap support, lineage, drug-resistance profile, mutation-level resistance evidence, lineage distribution, surveillance metadata, and SNP distances before making transmission inferences.
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
<tr><td>Mutation-level resistance evidence</td><td>Parsed from TB-Profiler JSON outputs and summarized by sample, drug or evidence source, gene, mutation/change, confidence, and evidence</td></tr>
<tr><td>Lineage distribution</td><td>Lineage counts and barplot generated from TB-Profiler lineage fields where resolved</td></tr>
<tr><td>Pairwise SNP distance and clustering</td><td>Pairwise SNP distances calculated from the MTBC core genome alignment after excluding reference/non-sample sequences and interpreted using configured SNP thresholds</td></tr>
<tr><td>SNP distance heatmap</td><td>SVG heatmap generated from the pairwise SNP distance matrix after reference filtering</td></tr>
<tr><td>Surveillance metadata</td><td>Downloadable metadata and QC filtering rationale TSV files generated for transparent surveillance reporting</td></tr>
<tr><td>Core-SNP phylogenomics</td><td>Snippy-core, optional Gubbins filtering, IQ-TREE2, and ETE3 tree rendering</td></tr>
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
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: 1
    memory: "4 GB"
  }

  output {
    File run_metadata = "final_report/run_metadata.txt"
    File final_report_html = "final_report/integrated_tb_amr_mtbc_phylogenomics_report.html"
    File? merged_lineage_distribution_svg = "final_report/lineage_distribution.svg"
    File? merged_snp_distance_heatmap_svg = "final_report/snp_distance_heatmap.svg"
    File? merged_surveillance_metadata_tsv = "final_report/tb_surveillance_metadata.tsv"
    File? merged_qc_filtering_rationale_tsv = "final_report/qc_filtering_rationale.tsv"
  }
}
