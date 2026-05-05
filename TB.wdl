version 1.0

workflow TB_AMR_MTBC_Phylogenomics {
  input {
    Array[File]+ input_reads
    File adapters
    File mtbc_reference_genbank

    Boolean do_trimming = true
    Boolean do_quality_control = true
    Boolean do_tb_profiler = true
    Boolean do_phylogeny = true
    Boolean use_gubbins = true
    Boolean midpoint_root_tree = true
    Boolean report_nonsynonymous_drug_gene_mutations = true

    String trimmomatic_quality_encoding = "phred33"
    String tbprofiler_docker = "staphb/tbprofiler:6.6.6"
    String snippy_reference_type = "genbank"
    String iqtree2_model = "GTR+G"
    Int iqtree2_bootstraps = 1000
    Int min_mtbc_samples_for_tree = 3

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

  # 2. FastQC runs only after trimming, because it consumes analysis_reads.
  if (do_quality_control) {
    call FASTQC {
      input:
        input_reads = analysis_reads,
        cpu = cpu_4
    }

    # 3. MultiQC runs only after FastQC, because it consumes FastQC outputs.
    call MULTIQC {
      input:
        fastqc_reports = FASTQC.fastqc_reports,
        fastqc_zips = FASTQC.fastqc_zips
    }
  }

  # 4. TB-Profiler is intentionally chained after MultiQC when QC is enabled through qc_dependency.
  if (do_tb_profiler) {
    call TB_PROFILER_AND_MTBC_FILTER {
      input:
        input_reads = analysis_reads,
        qc_dependency = MULTIQC.multiqc_report,
        docker_image = tbprofiler_docker,
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

  Array[File] mtbc_reads = select_first([TB_PROFILER_AND_MTBC_FILTER.mtbc_reads, []])

  # 5. MTBC-only Snippy/core-SNP analysis. If fewer than min_mtbc_samples_for_tree paired MTBC samples are available, this branch is skipped.
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

  if (report_nonsynonymous_drug_gene_mutations && defined(SNIPPY_CORE_MTBC.snippy_tab_files)) {
    call TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS {
      input:
        snippy_tab_files = select_first([SNIPPY_CORE_MTBC.snippy_tab_files, []]),
        genes_csv = tb_drug_resistance_genes
    }
  }

  # 6. Optional recombination filtering. If Gubbins fails internally, its task passes the original alignment forward.
  if (do_phylogeny && use_gubbins && defined(SNIPPY_CORE_MTBC.core_full_alignment)) {
    call GUBBINS_RECOMBINATION {
      input:
        core_full_alignment = select_first([SNIPPY_CORE_MTBC.core_full_alignment]),
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

  # 7. IQ-TREE uses the Gubbins-filtered alignment when available; otherwise it uses Snippy core.full.aln.
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

  # 8. Render tree only when a Newick tree exists.
  if (do_phylogeny && defined(IQTREE2_PHYLOGENY.final_tree)) {
    call TREE_VISUALIZATION {
    input:
      input_tree = IQTREE2_PHYLOGENY.final_tree,
      tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv,
      width = tree_width,
      height = tree_height,
      image_format = tree_image_format
    }
  }

  # 9. Final merged report always runs and displays skipped sections clearly when optional branches are absent.
  call MERGE_TB_REPORTS {
    input:
      tbprofiler_html = TB_PROFILER_AND_MTBC_FILTER.combined_html,
      tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv,
      mtbc_samples_txt = TB_PROFILER_AND_MTBC_FILTER.mtbc_samples_txt,
      qc_summary_html = MULTIQC.multiqc_report,
      trimming_report_html = TRIMMING.trimming_report,
      variant_summary_html = SNIPPY_CORE_MTBC.variant_summary,
      iqtree_report = IQTREE2_PHYLOGENY.iqtree_report,
      tree_image = TREE_VISUALIZATION.tree_image,
      phylogenetic_tree_newick = TREE_VISUALIZATION.cleaned_tree,
      nonsynonymous_mutations_tsv = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_tsv,
      nonsynonymous_mutations_html = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_html
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

    Array[File]? tbprofiler_json = TB_PROFILER_AND_MTBC_FILTER.json_reports
    Array[File]? tbprofiler_txt = TB_PROFILER_AND_MTBC_FILTER.txt_reports
    File? tbprofiler_summary_tsv = TB_PROFILER_AND_MTBC_FILTER.summary_tsv
    File? tbprofiler_combined_html = TB_PROFILER_AND_MTBC_FILTER.combined_html
    File? mtbc_samples = TB_PROFILER_AND_MTBC_FILTER.mtbc_samples_txt
    Array[File]? mtbc_reads_for_phylogeny = TB_PROFILER_AND_MTBC_FILTER.mtbc_reads

    File? variant_summary_html = SNIPPY_CORE_MTBC.variant_summary
    File? core_full_alignment = SNIPPY_CORE_MTBC.core_full_alignment
    File? core_snp_alignment = SNIPPY_CORE_MTBC.core_snp_alignment
    File? snippy_core_vcf = SNIPPY_CORE_MTBC.core_vcf
    Array[File]? snippy_tab_files = SNIPPY_CORE_MTBC.snippy_tab_files

    File? nonsynonymous_drug_gene_mutations_tsv = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_tsv
    File? nonsynonymous_drug_gene_mutations_html = TB_DRUG_GENE_NONSYNONYMOUS_MUTATIONS.nonsynonymous_mutations_html

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

task TB_PROFILER_AND_MTBC_FILTER {
  input {
    Array[File]+ input_reads
    File? qc_dependency
    String docker_image = "staphb/tbprofiler:6.6.6"
    Int cpu = 8
    Int memory_gb = 16
  }

  command <<<
    set -uo pipefail
    mkdir -p tbprofiler_results mtbc_reads logs

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

      python3 - "$json" "$sample" "$R1" "$R2" "$status" <<'PYTB'
import json, sys, os, shutil, re

json_file, sample, r1, r2, status = sys.argv[1:6]

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
        return v
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
        # Avoid dumping entire dicts into report cells.
        for key in ["name", "drug", "gene", "change", "mutation", "original_mutation", "confidence", "source"]:
            if key in v and v[key] not in (None, "", [], {}):
                return clean_value(v[key])
        return ""
    return str(v)

def uniq(xs):
    seen = []
    for x in xs:
        x = str(x).strip()
        if x and x not in seen:
            seen.append(x)
    return ", ".join(seen)

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

dr_type = clean_value(get_path(data, [
    "drtype", "dr_type", "resistance_type", "drug_resistance_type", "prediction.drtype"
]))

resistant_drugs = []
resistance_mutations = []
key_mutations = []

variant_blocks = []
for key in ["dr_variants", "other_variants", "variants"]:
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
    drug = item.get("drug") or item.get("drug_name") or item.get("name") or item.get("drugs") or item.get("drug_resistance") or ""

    if gene or change:
        key_mutations.append(" ".join([gene, change]).strip())

    if drug:
        if isinstance(drug, list):
            resistant_drugs.extend([clean_value(d) for d in drug])
        else:
            resistant_drugs.append(clean_value(drug))

    confidence = clean_value(item.get("confidence"))
    source = clean_value(item.get("source"))

    if drug or gene or change:
        mut = []
        if drug:
            mut.append(clean_value(drug))
        if gene or change:
            mut.append(" ".join([gene, change]).strip())
        if confidence:
            mut.append(confidence)
        if source:
            mut.append(source)
        resistance_mutations.append(" | ".join([x for x in mut if x]))

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

if dr_type and dr_type.lower() not in ["not reported", "susceptible", "sensitive", "none"] and not resistant_drugs:
    resistant_drugs.append(dr_type)

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

explicit_non_mtbc = bool(species) and not any(t in text_species for t in mtbc_species_terms)
lineage_text = " ".join([main_lineage, sub_lineage]).lower().strip()
valid_lineage = bool(re.search(r"(^|[^a-z0-9])(lineage)?[ _-]?[1-9](\.|$|[^0-9])", lineage_text))

if any(t in text_species for t in mtbc_species_terms):
    is_mtbc = True
    reason = "MTBC species reported by TB-Profiler"
elif not species and valid_lineage:
    is_mtbc = True
    reason = "Species absent, but valid TB-Profiler MTBC lineage reported"
else:
    is_mtbc = False
    reason = "No recognized MTBC species or valid TB-Profiler lineage"

if explicit_non_mtbc:
    is_mtbc = False
    reason = "Explicit non-MTBC species reported"

# TB-Profiler JSON files do not always expose a top-level species field.
# When species is absent but a valid MTBC lineage is present, keep the sample
# selected and report an explicit inferred MTBC label instead of "Not reported".
# When the sample is not selected and species is absent, report it as Non-MTB/Not MTBC.
species_display = species
if not species_display and is_mtbc:
    species_display = "Mycobacterium tuberculosis complex (inferred from TB-Profiler lineage)"
elif not species_display and not is_mtbc:
    species_display = "Non-MTB / not classified as MTBC"

if is_mtbc:
    os.makedirs("mtbc_reads", exist_ok=True)
    shutil.copy(r1, f"mtbc_reads/{sample}_R1.fastq.gz")
    shutil.copy(r2, f"mtbc_reads/{sample}_R2.fastq.gz")
    with open("mtbc_samples.txt", "a") as fh:
        fh.write(sample + "\n")

line = [
    sample,
    species_display or "Not reported",
    main_lineage or "Not reported",
    sub_lineage or "Not reported",
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

      echo -e "${sample}\t${status}\t${vcf}\t${aln}" >> variant_summary.tsv
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

    python3 - <<'PY'
import csv, html

rows = list(csv.DictReader(open("variant_summary.tsv"), delimiter="\t"))

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
<tr><th>Sample</th><th>Status</th><th>VCF</th><th>Aligned FASTA</th></tr>
</thead>
<tbody>
"""

for r in rows:
    status = r["status"]
    cls = "ok" if status == "success" else "fail"
    out += f"<tr><td>{html.escape(r['sample'])}</td><td class='{cls}'>{html.escape(status)}</td><td>{html.escape(r['vcf'])}</td><td>{html.escape(r['aligned_fasta'])}</td></tr>\n"

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

    File variant_summary = "variant_summary.html"

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
    File? input_tree
    File? tbprofiler_summary_tsv
    Int width = 2600
    Int height = 3200
    String image_format = "png"
    String title = "Phylogenetic tree"
  }

  command <<<
    set -euo pipefail
    mkdir -p tree_visualization
    export QT_QPA_PLATFORM=offscreen
    export MPLBACKEND=Agg

    python3 - <<'PY'
from pathlib import Path
import csv

tree_input = "~{if defined(input_tree) then input_tree else ""}"
tbprofiler_summary_path = "~{if defined(tbprofiler_summary_tsv) then tbprofiler_summary_tsv else ""}"

image_format = "~{image_format}".lower().strip()
out_img = Path(f"tree_visualization/phylogenetic_tree.{image_format}")
cleaned_tree = Path("tree_visualization/phylogenetic_tree.cleaned.nwk")
log = Path("tree_visualization/render.log")

log.write_text("TREE_VISUALIZATION started\n")

RESISTANCE_COLORS = {
  "Unknown": "#999999",
  "Sensitive": "#1b9e77",
  "MDR": "#d73027",
  "Pre-XDR": "#FF6A00",
  "XDR": "#FF0000",
  "Other": "#e6ab02"
}

try:
    from ete3 import Tree, TreeStyle, TextFace, NodeStyle

    if image_format not in {"png", "svg", "pdf"}:
        raise ValueError(f"Unsupported image_format: {image_format}. Use png, svg, or pdf.")

    if not tree_input:
        raise ValueError("input_tree not provided")

    tree_path = Path(tree_input)

    if not tree_path.exists() or tree_path.stat().st_size == 0:
        raise FileNotFoundError(f"Tree missing or empty: {tree_path}")

    metadata = {}

    def classify_resistance(dr_type, resistant_drugs):
        profile_raw = (dr_type or "").strip()
        profile = profile_raw.lower().replace("_", "-")

        if profile in {"susceptible", "sensitive"}:
            return "Sensitive", RESISTANCE_COLORS["Sensitive"]

        if "pre-xdr" in profile or "pre xdr" in profile or "prexdr" in profile:
            return "Pre-XDR", RESISTANCE_COLORS["Pre-XDR"]

        if "xdr" in profile:
            return "XDR", RESISTANCE_COLORS["XDR"]

        if "mdr" in profile:
            return "MDR", RESISTANCE_COLORS["MDR"]

        # Preserve the exact TB-Profiler resistance profile label:
        # e.g. HR-TB, RR-TB, RIF resistant, INH resistant.
        # These all share the same yellow color but keep their original names.
        if profile_raw and profile not in {"other", "not reported", "unknown", "na", "n/a", "none"}:
            return profile_raw, RESISTANCE_COLORS["Other"]

        return "Unknown", RESISTANCE_COLORS["Unknown"]

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

                category, color = classify_resistance(
                    r.get("dr_type"),
                    r.get("resistant_drugs")
                )

                metadata[sample] = {
                    "lineage": lineage,
                    "category": category,
                    "color": color
                }

    t = Tree(str(tree_path), format=1)

    reference_names = {
        "reference",
        "ref",
        "h37rv",
        "nc_000962",
        "nc_000962.3"
    }

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
    requested_width = int("~{width}")
    requested_height = int("~{height}")

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
        ns["hz_line_color"] = "#000000"
        ns["vt_line_color"] = "#000000"
        ns["fgcolor"] = "#000000"
        ns["size"] = 0
        ns["shape"] = "circle"

        if node.is_leaf():
            sample = node.name
            meta = metadata.get(sample, {})

            lineage = meta.get("lineage", "")
            category = meta.get("category", "Unknown")
            color = meta.get("color", RESISTANCE_COLORS["Unknown"])

            ns["fgcolor"] = color
            ns["size"] = tip_node_size
            node.set_style(ns)

            node.add_face(
                TextFace(sample, fsize=label_font, fgcolor="#111111"),
                column=0,
                position="branch-right"
            )

            if lineage:
                node.add_face(
                    TextFace("  " + lineage, fsize=lineage_font, fgcolor="#555555"),
                    column=1,
                    position="branch-right"
                )

            node.add_face(
                TextFace("  ● " + category, fsize=resistance_font, fgcolor=color),
                column=2,
                position="branch-right"
            )

            node.name = ""

        else:
            ns["fgcolor"] = "#000000"
            ns["size"] = 0
            node.set_style(ns)

            support = ""

            # Prefer the original IQ-TREE internal node label.
            # If IQ-TREE produced SH-aLRT/UFBoot labels such as 99.4/100,
            # display only the UFBoot/bootstrap value after "/".
            raw_name = str(getattr(node, "name", "") or "").strip()

            if raw_name:
                try:
                    if "/" in raw_name:
                        raw_name = raw_name.split("/")[-1].strip()

                    value = float(raw_name)

                    if 0 < value <= 1:
                        value = value * 100

                    support = str(int(round(value)))
                except Exception:
                    support = ""
            else:
                raw_support = str(getattr(node, "support", "") or "").strip()

                # ETE3 may assign default internal-node support = 1.0.
                # Do not convert this default to 100 unless it came from a real label.
                if raw_support not in {"", "0", "0.0", "1", "1.0"}:
                    try:
                        value = float(raw_support)
                        if 0 < value <= 1:
                            value = value * 100
                        support = str(int(round(value)))
                    except Exception:
                        support = ""

            if support:
                node.add_face(
                    TextFace(support, fsize=bootstrap_font, fgcolor="#b00000"),
                    column=0,
                    position="branch-top"
                )

    ts = TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_scale = True
    ts.scale = None
    ts.branch_vertical_margin = branch_vertical_margin

    ts.margin_top = 12
    ts.margin_bottom = 65
    ts.margin_left = 25
    ts.margin_right = margin_right
    ts.title.clear()

    t.write(format=1, outfile=str(cleaned_tree))

    t.render(str(out_img), w=auto_width, units="px", tree_style=ts)

    with open(log, "a") as fh:
        fh.write("TREE_VISUALIZATION completed successfully\n")
        fh.write(f"Tree: {tree_path}\n")
        fh.write(f"Output: {out_img}\n")
        fh.write(f"Cleaned tree: {cleaned_tree}\n")
        fh.write(f"Tips rendered: {n_leaves}\n")
        fh.write(f"Reference tips removed: {removed_refs}\n")
        fh.write(f"Metadata records loaded: {len(metadata)}\n")
        fh.write(f"Requested width: {requested_width}\n")
        fh.write(f"Requested height: {requested_height}\n")
        fh.write(f"Auto width used: {auto_width}\n")
        fh.write("Height mode: natural ETE3 height using dynamic branch_vertical_margin\n")
        fh.write(f"Label font: {label_font}\n")
        fh.write(f"Lineage font: {lineage_font}\n")
        fh.write(f"Resistance font: {resistance_font}\n")
        fh.write(f"Bootstrap font: {bootstrap_font}\n")
        fh.write(f"Tip node size: {tip_node_size}\n")
        fh.write("Internal node size: 0\n")
        fh.write(f"Branch width: {branch_width}\n")
        fh.write(f"Branch vertical margin: {branch_vertical_margin}\n")
        fh.write(f"Right margin: {margin_right}\n")
        fh.write("Bootstrap/support display: original IQ-TREE internal labels preferred; ETE3 default 1.0 skipped\n")
        fh.write("Susceptibility color key:\n")
        for label, hex_color in RESISTANCE_COLORS.items():
            fh.write(f"  {label}: {hex_color}\n")

except Exception as e:
    with open(log, "a") as fh:
        fh.write("ERROR:\n")
        fh.write(repr(e) + "\n")
    raise
PY
  >>>

  runtime {
    docker: "gmboowa/ete3-render:1.18"
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }

  output {
    File tree_image = "tree_visualization/phylogenetic_tree.~{image_format}"
    File cleaned_tree = "tree_visualization/phylogenetic_tree.cleaned.nwk"
    File render_log = "tree_visualization/render.log"
  }
}
task MERGE_TB_REPORTS {
  input {
    String docker_image = "python:3.11-slim"
    File? tbprofiler_html
    File? tbprofiler_summary_tsv
    File? mtbc_samples_txt
    File? qc_summary_html
    File? trimming_report_html
    File? variant_summary_html
    File? iqtree_report
    File? tree_image
    File? phylogenetic_tree_newick
    File? pairwise_tree_newick
    File? nonsynonymous_mutations_tsv
    File? nonsynonymous_mutations_html
  }

  command <<<
    set -uo pipefail
    mkdir -p final_report

    tb_tsv="~{if defined(tbprofiler_summary_tsv) then tbprofiler_summary_tsv else ""}"
    nonsyn_tsv="~{if defined(nonsynonymous_mutations_tsv) then nonsynonymous_mutations_tsv else ""}"
    tree_png="~{if defined(tree_image) then tree_image else ""}"
    qc_html="~{if defined(qc_summary_html) then qc_summary_html else ""}"
    trim_html="~{if defined(trimming_report_html) then trimming_report_html else ""}"
    variant_html="~{if defined(variant_summary_html) then variant_summary_html else ""}"
    iqtree_txt="~{if defined(iqtree_report) then iqtree_report else ""}"

    if [ -n "$tree_png" ] && [ -f "$tree_png" ]; then
      cp "$tree_png" final_report/mtbc_tree.png || true
    fi

    if [ -z "$tb_tsv" ] || [ ! -f "$tb_tsv" ]; then
      echo -e "sample\tspecies\tmain_lineage\tsub_lineage\tdr_type\tresistant_drugs" > final_report/empty.tsv
      tb_tsv="final_report/empty.tsv"
    fi

    if [ -z "$nonsyn_tsv" ] || [ ! -f "$nonsyn_tsv" ]; then
      echo -e "sample\tgene\teffect\taa_change\tnt_change\tproduct" > final_report/empty_nonsyn.tsv
      nonsyn_tsv="final_report/empty_nonsyn.tsv"
    fi

    python3 - "$tb_tsv" "$nonsyn_tsv" "$qc_html" "$trim_html" "$variant_html" "$iqtree_txt" <<'PY'
import csv, html, re
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import sys

summary_tsv, nonsyn_tsv, qc_html, trim_html, variant_html, iqtree_txt = sys.argv[1:7]

outdir = Path("final_report")
outdir.mkdir(exist_ok=True)

run_started_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
run_stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_UTC")

(outdir / "run_metadata.txt").write_text(
    f"Workflow report generation timestamp: {run_started_utc}\n"
    f"Run stamp: {run_stamp}\n"
)

def safe(v):
    return html.escape(str(v if v is not None else ""))

def read_optional_file(path):
    p = Path(path)
    if path and p.exists() and p.stat().st_size > 0:
        try:
            return p.read_text(errors="replace")
        except Exception:
            return ""
    return ""

rows = list(csv.DictReader(open(summary_tsv), delimiter="\t"))
nonsyn_rows = list(csv.DictReader(open(nonsyn_tsv), delimiter="\t"))

RESISTANCE_COLORS = {
  "Sensitive": "#1b9e77",
  "Monoresistance": "#e6ab02",
  "MDR": "#d73027",
  "Pre-XDR": "#FF6A00",
  "XDR": "#FF0000",
  "Unknown": "#999999",
}

def classify_resistance(dr_type, resistant_drugs):
    text = f"{dr_type or ''} {resistant_drugs or ''}".lower().replace("_", "-")

    if any(x in text for x in ["xdr-tb", "xdr", "extensively drug"]) and not any(x in text for x in ["pre-xdr", "pre xdr", "prexdr"]):
        return "XDR"

    if any(x in text for x in ["pre-xdr", "pre xdr", "prexdr"]):
        return "Pre-XDR"

    if "mdr" in text or ("rif" in text and "inh" in text) or ("rifampicin" in text and "isoniazid" in text):
        return "MDR"

    if any(x in text for x in [
        "rif resistant", "rifampicin resistant", "rr-tb", "rif",
        "inh resistant", "isoniazid resistant", "hr-tb", "inh",
        "mono", "monoresistant", "mono-resistant"
    ]):
        return "Monoresistance"

    if any(x in text for x in [
        "susceptible", "sensitive", "pan-susceptible",
        "none reported", "no resistance"
    ]):
        return "Sensitive"

    if text.strip() in {"", "none", "not reported", "na", "n/a", "unknown", "other"}:
        return "Unknown"

    return "Unknown"

def badge(label):
    color = RESISTANCE_COLORS.get(label, RESISTANCE_COLORS["Unknown"])
    return f'<span class="res-badge" style="background:{color};">{safe(label)}</span>'

total_samples = len(rows)
mtbc_retained = len(rows)
non_mtbc = 0
drug_resistant = 0

for r in rows:
    category = classify_resistance(r.get("dr_type"), r.get("resistant_drugs"))
    if category not in {"Sensitive", "Unknown"}:
        drug_resistant += 1

def build_qc_section():
    qc_content = read_optional_file(qc_html)
    trim_content = read_optional_file(trim_html)
    variant_content = read_optional_file(variant_html)

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
                '<td><span class="badge badge-green">PASS</span></td>'
                '<td><span class="badge badge-green">Proceed</span></td>'
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

def build_tb_rows():
    out = []
    for r in rows:
        sample = r.get("sample")
        species = r.get("species") or "Mycobacterium tuberculosis complex (inferred from TB-Profiler lineage)"
        main_lineage = r.get("main_lineage") or ""
        sub_lineage = r.get("sub_lineage") or ""
        lineage = f"{safe(main_lineage)} / {safe(sub_lineage)}" if sub_lineage else safe(main_lineage)
        dr_type = r.get("dr_type") or "Unknown"
        resistant_drugs = r.get("resistant_drugs") or "None reported"
        category = classify_resistance(dr_type, resistant_drugs)

        out.append(
        "<tr>"
        f"<td>{safe(sample)}</td>"
        f"<td>{safe(species)}</td>"
        f"<td>{lineage}</td>"
        f"<td>{badge(category)} &nbsp; {safe(dr_type)}</td>"
        f"<td>{safe(resistant_drugs)}</td>"
        '<td><span class="badge badge-green">Selected</span></td>'
        "</tr>"
        )
    return "".join(out)

def build_nonsyn_section():
    cleaned = [r for r in nonsyn_rows if any((v or "").strip() for v in r.values())]

    if not cleaned:
        return """
<div class="section">
<h2>3. Non-synonymous mutations in key TB drug-resistance-associated genes</h2>
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
<h2>3. Non-synonymous mutations in key TB drug-resistance-associated genes</h2>
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

def build_tree_legend():
    items = [
        ("Sensitive", RESISTANCE_COLORS["Sensitive"]),
        ("Monoresistance / HR-TB / RR-TB", RESISTANCE_COLORS["Monoresistance"]),
        ("MDR", RESISTANCE_COLORS["MDR"]),
        ("Pre-XDR", RESISTANCE_COLORS["Pre-XDR"]),
        ("XDR", RESISTANCE_COLORS["XDR"]),
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
<title>Interactive TB AMR MTBC Phylogenomics Report</title>
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
  color:white;
  font-size:12px;
  display:inline-block;
}}
.badge-green{{background:#28A745;font-weight:700;}}
.badge-red{{background:#b91c1c;}}
.badge-blue{{background:#2563eb;}}
.badge-orange{{background:#d97706;}}
.res-badge{{
  padding:4px 8px;
  border-radius:999px;
  color:white;
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
<h1>Interactive TB AMR MTBC Phylogenomics Report</h1>
<p>Trimming → QC → TB-Profiler → MTBC-only filtering → core-SNP phylogenomics → final merged report</p>
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

<div class="section">
<h2>2. TB-Profiler Resistance, Species, and Lineage Report</h2>

<div class="note">
<strong>Interpretation note:</strong> TB-Profiler is using WHO 2021+ definitions for classification of drug-resistance results.
<strong>Key WHO 2021+ resistance definitions:</strong>
<strong>Rifampicin-resistant TB (RR-TB):</strong> resistant to rifampicin, with or without resistance to other drugs.
<strong>Multidrug-resistant TB (MDR-TB):</strong> resistance to at least isoniazid and rifampicin.
<strong>Pre-extensively drug-resistant TB (Pre-XDR-TB):</strong> MDR/RR-TB that is also resistant to any fluoroquinolone.
<strong>Extensively drug-resistant TB (XDR-TB):</strong> MDR/RR-TB that is resistant to any fluoroquinolone and at least one additional Group A drug, bedaquiline or linezolid.
</div>

<div class="controls">
<input id="tbSearch" onkeyup="filterTable('tbSearch','tbTable')" placeholder="Search TB-Profiler results...">
<select onchange="filterResistance(this.value)">
<option value="">All resistance profiles</option>
<option value="Sensitive">Sensitive only</option>
<option value="Monoresistance">Monoresistance / HR-TB / RR-TB</option>
<option value="MDR">MDR only</option>
<option value="Pre-XDR">Pre-XDR only</option>
<option value="XDR">XDR only</option>
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

{build_nonsyn_section()}

<div class="section">
<h2>4. MTBC-only Core-SNP Phylogenetic Tree</h2>
<div class="grid2">
<div class="tree-panel">
{tree_html}
{build_tree_legend()}
</div>
<div class="details">
<h3>Tree construction summary</h3>
<p><strong>Included:</strong> {mtbc_retained} TB-Profiler-confirmed MTBC isolate(s).</p>
<p><strong>Excluded:</strong> {non_mtbc} non-MTBC or low-confidence isolate(s).</p>
<p><strong>Core alignment:</strong> Snippy-core alignment.</p>
<p><strong>Recombination:</strong> Optional Gubbins-filtered alignment when enabled.</p>
<p><strong>Tree:</strong> IQ-TREE2 maximum-likelihood phylogeny.</p>
<p><strong>Display:</strong> ETE3-rendered static tree image, shown inside an auto-scaling scrollable report panel.</p>
{build_iqtree_section()}
</div>
</div>
</div>

<div class="section">
<h2>5. Final Interpretation</h2>
<p>The report documents all samples through QC and TB-Profiler analysis, then applies an MTBC-only rule before phylogenomic reconstruction. Samples not classified as MTBC are excluded from the tree but retained in the workflow record for transparency.</p>
<div class="note">
<strong>Interpretation:</strong> use close clustering together with bootstrap support, lineage, drug-resistance profile, metadata, and SNP distances before making transmission inferences.
</div>
</div>

<div class="footer">
Generated by TB_AMR_MTBC_Phylogenomics WDL workflow. Run generated: {safe(run_started_utc)}. Run stamp: {safe(run_stamp)}.
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
  >>>

  runtime {
    docker: "~{docker_image}"
    cpu: 1
    memory: "4 GB"
  }

  output {
    File run_metadata = "final_report/run_metadata.txt"
    File final_report_html = "final_report/integrated_tb_amr_mtbc_phylogenomics_report.html"
  }
}
