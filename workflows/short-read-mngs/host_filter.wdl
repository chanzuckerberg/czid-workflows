version 1.0

# CZ ID short-read-mngs pipeline stage 1 (2022 version):
# - input validation & QC
# - host & human filtering
# - deduplication
# - subsampling
workflow czid_host_filter {
  input {
    File fastqs_0
    File? fastqs_1
    String nucleotide_type = "DNA"

    File adapter_fasta

    String host_genome
    File bowtie2_index_tar
    File hisat2_index_tar
    File kallisto_idx

    File human_bowtie2_index_tar
    File human_hisat2_index_tar

    Int max_input_fragments
    Int max_subsample_fragments

    Int cpu = 16
    String docker_image_id

    # legacy idseq-dag inputs:
    String file_ext = "fastq"
    String s3_wd_uri
  }

  # Validate input reads (and truncate if very large)
  call RunValidateInput {
    input:
    reads1_fastq = fastqs_0,
    reads2_fastq = fastqs_1,
    file_ext = file_ext,
    max_input_fragments = max_input_fragments,
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri
  }

  # Adapter trimming and QC filtering
  call fastp_qc {
    input:
    reads1_fastq = RunValidateInput.valid_input1_fastq,
    reads2_fastq = RunValidateInput.valid_input2_fastq,
    adapter_fasta = adapter_fasta,
    docker_image_id = docker_image_id,
    cpu = cpu
  }

  # If RNAseq, quantify host transcripts and ERCC
  if (nucleotide_type == "RNA") {
    call kallisto {
      input:
      reads1_fastq = fastp_qc.fastp1_fastq,
      reads2_fastq = fastp_qc.fastp2_fastq,
      kallisto_idx = kallisto_idx,
      docker_image_id = docker_image_id,
      cpu = cpu
    }
  }

  # Filter out host reads.
  # Two stages: bowtie2 --very-sensitive-local, followed by splice-aware HISAT2.
  call bowtie2_filter {
    input:
    reads1_fastq = fastp_qc.fastp1_fastq,
    reads2_fastq = fastp_qc.fastp2_fastq,
    index_tar = bowtie2_index_tar,
    docker_image_id = docker_image_id,
    cpu = cpu
  }

  call hisat2_filter {
    input:
    reads1_fastq = bowtie2_filter.filtered1_fastq,
    reads2_fastq = bowtie2_filter.filtered2_fastq,
    index_tar = hisat2_index_tar,
    docker_image_id = docker_image_id,
    cpu = cpu
  }

  # If paired-end, collect insert size metrics from unfiltered, host-aligned bowtie2 BAM.
  if (defined(fastqs_1)) {
    call collect_insert_size_metrics {
      input:
      bam = bowtie2_filter.bam,
      docker_image_id = docker_image_id
    }
  }

  # Additionally filter out human reads, if the host is non-human.
  if (host_genome != "human") {
    call bowtie2_filter as bowtie2_human_filter {
      input:
      reads1_fastq = hisat2_filter.filtered1_fastq,
      reads2_fastq = hisat2_filter.filtered2_fastq,
      filter_type = "human",
      index_tar = human_bowtie2_index_tar,
      docker_image_id = docker_image_id,
      cpu = cpu
    }
    call hisat2_filter as hisat2_human_filter {
      input:
      reads1_fastq = bowtie2_human_filter.filtered1_fastq,
      reads2_fastq = bowtie2_human_filter.filtered2_fastq,
      filter_type = "human",
      index_tar = human_hisat2_index_tar,
      docker_image_id = docker_image_id,
      cpu = cpu
    }
  }

  # Collect effectively filtered reads from the previous conditional
  File filtered1_fastq = select_first([hisat2_human_filter.filtered1_fastq, hisat2_filter.filtered1_fastq])
  File? filtered2_fastq = if defined(hisat2_human_filter.filtered2_fastq) then hisat2_human_filter.filtered2_fastq
                                                                          else hisat2_filter.filtered2_fastq

  # Deduplicate filtered reads using custom czid-dedup tool.
  # It retains one exemplar [pair] from each duplicate cluster, and produces mapping from exemplar
  # read name to cluster size.
  call RunCZIDDedup {
    input:
    reads1_fastq = filtered1_fastq,
    reads2_fastq = filtered2_fastq,
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri
  }

  # Subsample remaining reads.
  call RunSubsample {
    input:
    reads1_fastq = RunCZIDDedup.dedup1_fastq,
    reads2_fastq = RunCZIDDedup.dedup2_fastq,
    duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
    max_subsample_fragments = max_subsample_fragments,
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri
  }

  output {
    File reads_in_count = RunValidateInput.reads_in_count
    File validate_input_out_validate_input_summary_json = RunValidateInput.validate_input_summary_json
    File validate_input_out_count = RunValidateInput.reads_out_count

    File fastp_out_fastp1_fastq = fastp_qc.fastp1_fastq
    File? fastp_out_fastp2_fastq = fastp_qc.fastp2_fastq
    File fastp_out_count = fastp_qc.reads_out_count
    File fastp_html = fastp_qc.fastp_html
    File fastp_json = fastp_qc.fastp_json

    File? kallisto_abundance_tsv = kallisto.abundance_tsv

    File bowtie2_host_filtered1_fastq = bowtie2_filter.filtered1_fastq
    File? bowtie2_host_filtered2_fastq = bowtie2_filter.filtered2_fastq
    File bowtie2_host_filtered_out_count = bowtie2_filter.reads_out_count
    File hisat2_host_filtered1_fastq = hisat2_filter.filtered1_fastq
    File? hisat2_host_filtered2_fastq = hisat2_filter.filtered2_fastq
    File hisat2_host_filtered_out_count = hisat2_filter.reads_out_count

    File? insert_size_metrics = collect_insert_size_metrics.insert_size_metrics
    File? insert_size_histogram = collect_insert_size_metrics.insert_size_histogram

    File? bowtie2_human_filtered1_fastq = bowtie2_human_filter.filtered1_fastq
    File? bowtie2_human_filtered2_fastq = bowtie2_human_filter.filtered2_fastq
    File? bowtie2_human_filtered_out_count = bowtie2_human_filter.reads_out_count
    File? hisat2_human_filtered1_fastq = hisat2_human_filter.filtered1_fastq
    File? hisat2_human_filtered2_fastq = hisat2_human_filter.filtered2_fastq
    File? hisat2_human_filtered_out_count = hisat2_human_filter.reads_out_count

    File czid_dedup_out_dedup1_fastq = RunCZIDDedup.dedup1_fastq
    File? czid_dedup_out_dedup2_fastq = RunCZIDDedup.dedup2_fastq
    File czid_dedup_out_duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv
    File czid_dedup_out_duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv
    File czid_dedup_out_count = RunCZIDDedup.reads_out_count

    File subsampled_out_subsampled_1_fa = RunSubsample.subsampled_1_fa
    File? subsampled_out_subsampled_2_fa = RunSubsample.subsampled_2_fa
    File? subsampled_out_subsampled_merged_fa = RunSubsample.subsampled_merged_fa
    File subsampled_out_count = RunSubsample.reads_out_count
  }
}

task RunValidateInput {
  input {
    File reads1_fastq
    File? reads2_fastq
    String file_ext

    Int max_input_fragments

    String docker_image_id
    String s3_wd_uri
  }
  Boolean paired = defined(reads2_fastq)
  command<<<
    set -euxo pipefail
    idseq-dag-run-step --workflow-name host_filter \
      --step-module idseq_dag.steps.run_validate_input \
      --step-class PipelineStepRunValidateInput \
      --step-name validate_input_out \
      --input-files '[["~{sep='","' select_all([reads1_fastq, reads2_fastq])}"]]' \
      --output-files '["validate_input_summary.json", ~{if paired then '"valid_input1.fastq", "valid_input2.fastq"' else '"valid_input1.fastq"'}]' \
      --output-dir-s3 '~{s3_wd_uri}' \
      --additional-files '{}' \
      --additional-attributes '{"truncate_fragments_to": ~{max_input_fragments}, "file_ext": "~{file_ext}"}'
  >>>
  output {
    String step_description_md = read_string("validate_input_out.description.md")
    File validate_input_summary_json = "validate_input_summary.json"
    File valid_input1_fastq = "valid_input1.fastq"
    File? valid_input2_fastq = "valid_input2.fastq"
    File reads_out_count = "validate_input_out.count"
    File reads_in_count = "fastqs.count"
  }
  runtime {
    docker: docker_image_id
    cpu: 4
    memory: "8G"
  }
}

task fastp_qc {
    # fastp all-in-one for
    # - adapter trimming
    # - quality filtering
    # - complexity filtering
  input {
    File reads1_fastq
    File? reads2_fastq
    File adapter_fasta

    # These default QC thresholds are loosely based on the pre-2022 pipeline using PriceSeq & LZW
    String fastp_options = "--dont_eval_duplication --length_required 35" +
                           " --qualified_quality_phred 17 --unqualified_percent_limit 15 --n_base_limit 15" +
                           " --sdust_complexity_filter --complexity_threshold 60"

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(reads2_fastq)
  command<<<
    set -euxo pipefail
    fastp \
        -i ~{reads1_fastq} ~{"-I " + reads2_fastq} \
        -o fastp1.fastq ~{if (paired) then "-O fastp2.fastq" else ""} \
        -w ~{cpu} ~{fastp_options} \
        --adapter_fasta ~{adapter_fasta} ~{if (paired) then "--detect_adapter_for_pe" else ""}
    count="$(jq .read1_after_filtering.total_reads fastp.json)"
    if [ '~{paired}' == 'true' ]; then
        count=$((2 * count))
    fi
    jq --null-input --arg count "$count" '{"fastp_out":$count}' > fastp_out.count
    # TODO: extract insert size metrics from JSON, also render histogram?

    python3 - <<EOF
    import textwrap
    with open("fastp.description.md", "w") as outfile:
      print(textwrap.dedent("""
      # fastp read trimming & filtering

      *PLACEHOLDER TEXT*
      """).strip(), file=outfile)
    EOF
  >>>
  output {
    String step_description_md = read_string("fastp.description.md")
    File fastp1_fastq = "fastp1.fastq"
    File? fastp2_fastq = "fastp2.fastq"
    File fastp_html = "fastp.html"
    File fastp_json = "fastp.json"
    File reads_out_count = "fastp_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu}G"
  }
}

task kallisto {
  input {
    File reads1_fastq
    File? reads2_fastq
    File kallisto_idx
    String kallisto_options = ""

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(reads2_fastq)

  command <<<
    set -euxo pipefail

    single=""
    if [[ '~{paired}' != 'true' ]]; then
      # TODO: input fragment length parameters (l = average, s = std dev)
      single="--single -l 200 -s 20"
    fi
    # shellcheck disable=SC2086
    /kallisto/kallisto quant -i '~{kallisto_idx}' -o "$(pwd)" --plaintext $single ~{kallisto_options} -t ~{cpu} \
      ~{sep=' ' select_all([reads1_fastq, reads2_fastq])}
    >&2 jq . run_info.json

    python3 - <<EOF
    import textwrap
    with open("kallisto.description.md", "w") as outfile:
      print(textwrap.dedent("""
      # kallisto RNA quantification

      *PLACEHOLDER TEXT*
      """).strip(), file=outfile)
    EOF
  >>>

  output {
    String step_description_md = read_string("kallisto.description.md")
    File abundance_tsv = "abundance.tsv"
  }

  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu}G"
  }
}

task bowtie2_filter {
  # Remove reads [pairs] with bowtie2 hits to the given index
  input {
    File reads1_fastq
    File? reads2_fastq
    String filter_type = "host" # or human

    # GENOME_NAME.bowtie2.tar should contain GENOME_NAME/GENOME_NAME.*.bt*
    File index_tar
    String bowtie2_options = "--very-sensitive-local"

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(reads2_fastq)
  command <<<
    set -euxo pipefail
    TMPDIR="${TMPDIR:-/tmp}"

    genome_name="$(basename '~{index_tar}' .bowtie2.tar)"
    tar xf '~{index_tar}' -C "$TMPDIR"

    if [[ '~{paired}' == 'true' ]]; then
        bowtie2 -x "$TMPDIR/$genome_name/$genome_name" ~{bowtie2_options} -p ~{cpu} \
            -q -1 '~{reads1_fastq}' -2 '~{reads2_fastq}' \
            -S "$TMPDIR/bowtie2.sam"
    else
        bowtie2 -x "$TMPDIR/$genome_name/$genome_name" ~{bowtie2_options} -p ~{cpu} \
            -q -U '~{reads1_fastq}' \
            -S "$TMPDIR/bowtie2.sam"
    fi

    # generate sort & compressed BAM file for archival
    samtools sort -o "bowtie2_~{filter_type}.bam" -@ 4 -T "$TMPDIR" "$TMPDIR/bowtie2.sam" & samtools_pid=$!

    # Extract reads [pairs] that did NOT map to the index
    if [[ '~{paired}' == 'true' ]]; then
        #    1 (read paired)
        #    4 (read unmapped)
        # +  8 (mate unmapped)
        # ----
        #   13
        samtools fastq -f 13 -1 'bowtie2_~{filter_type}_filtered1.fastq' -2 'bowtie2_~{filter_type}_filtered2.fastq' -0 /dev/null -s /dev/null "$TMPDIR/bowtie2.sam"
    else
        samtools fastq -f 4 "$TMPDIR/bowtie2.sam" > 'bowtie2_~{filter_type}_filtered1.fastq'
    fi

    count="$(cat bowtie2_~{filter_type}_filtered{1,2}.fastq | wc -l)"
    count=$((count / 4))
    jq --null-input --arg count "$count" '{"bowtie2_~{filter_type}_filtered_out":$count}' > 'bowtie2_~{filter_type}_filtered_out.count'

    python3 - <<EOF
    import textwrap
    with open("bowtie2.description.md", "w") as outfile:
      print(textwrap.dedent("""
      # bowtie2 ~{filter_type} filtering

      *PLACEHOLDER TEXT*
      """).strip(), file=outfile)
    EOF

    wait $samtools_pid
  >>>

  output {
    String step_description_md = read_string("bowtie2.description.md")
    File filtered1_fastq = glob("bowtie2_*_filtered1.fastq")[0]
    File? filtered2_fastq = if paired then glob("bowtie2_*_filtered2.fastq")[0] else reads2_fastq
    File reads_out_count = "bowtie2_~{filter_type}_filtered_out.count"
    File bam = "bowtie2_~{filter_type}.bam"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu*2}G"
  }
}

task hisat2_filter {
  # Remove reads [pairs] with HISAT2 hits to the given index
  input {
    File reads1_fastq
    File? reads2_fastq
    String filter_type = "host" # or human

    # GENOME_NAME.hisat2.tar should contain GENOME_NAME/GENOME_NAME.*.ht2
    File index_tar
    String hisat2_options = ""

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(reads2_fastq)
  command <<<
    set -euxo pipefail
    TMPDIR="${TMPDIR:-/tmp}"

    genome_name="$(basename '~{index_tar}' .hisat2.tar)"
    tar xf '~{index_tar}' -C "$TMPDIR"

    if [[ '~{paired}' == 'true' ]]; then
        /hisat2/hisat2 -x "$TMPDIR/$genome_name/$genome_name" ~{hisat2_options} -p ~{cpu} \
            -q -1 '~{reads1_fastq}' -2 '~{reads2_fastq}' \
            -S "$TMPDIR/hisat2.sam"
    else
        /hisat2/hisat2 -x "$TMPDIR/$genome_name/$genome_name" ~{hisat2_options} -p ~{cpu} \
            -q -U '~{reads1_fastq}' \
            -S "$TMPDIR/hisat2.sam"
    fi

    # Extract reads [pairs] that did NOT map to the index
    if [[ '~{paired}' == 'true' ]]; then
        #    1 (read paired)
        #    4 (read unmapped)
        # +  8 (mate unmapped)
        # ----
        #   13
        samtools fastq -f 13 -1 'hisat2_~{filter_type}_filtered1.fastq' -2 'hisat2_~{filter_type}_filtered2.fastq' -0 /dev/null -s /dev/null "$TMPDIR/hisat2.sam"
    else
        samtools fastq -f 4 "$TMPDIR/hisat2.sam" > 'hisat2_~{filter_type}_filtered1.fastq'
    fi

    count="$(cat hisat2_~{filter_type}_filtered{1,2}.fastq | wc -l)"
    count=$((count / 4))
    jq --null-input --arg count "$count" '{"hisat2_~{filter_type}_filtered_out":$count}' > 'hisat2_~{filter_type}_filtered_out.count'

    python3 - <<EOF
    import textwrap
    with open("hisat2.description.md", "w") as outfile:
      print(textwrap.dedent("""
      # HISAT2 ~{filter_type} filtering

      *PLACEHOLDER TEXT*
      """).strip(), file=outfile)
    EOF
  >>>

  output {
    String step_description_md = read_string("hisat2.description.md")
    File filtered1_fastq = glob("hisat2_*_filtered1.fastq")[0]
    File? filtered2_fastq = if paired then glob("hisat2_*_filtered2.fastq")[0] else reads2_fastq
    File reads_out_count = "hisat2_~{filter_type}_filtered_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu*2}G"
  }
}

task collect_insert_size_metrics {
  input {
    File bam
    String docker_image_id
  }

  command <<<
    picard CollectInsertSizeMetrics 'I=~{bam}' O=insert_size_metrics.txt H=insert_size_histogram.pdf
    python3 - <<EOF
    import textwrap
    with open("collect_insert_size_metrics.description.md", "w") as outfile:
      print(textwrap.dedent("""
      # Picard CollectInsertSizeMetrics

      *PLACEHOLDER TEXT*
      """).strip(), file=outfile)
    EOF
  >>>

  output {
    String step_description_md = read_string("collect_insert_size_metrics.description.md")
    # If no reads mapped to the host, then picard exits "successfully" without creating these files.
    File? insert_size_metrics = "insert_size_metrics.txt"
    File? insert_size_histogram = "insert_size_histogram.pdf"
  }

  runtime {
    docker: docker_image_id
    cpu: 1
    memory: "4G"
  }
}

task RunCZIDDedup {
  input {
    File reads1_fastq
    File? reads2_fastq
    String docker_image_id
    String s3_wd_uri
  }
  Boolean paired = defined(reads2_fastq)
  command<<<
    set -euxo pipefail

    >&2 idseq-dag-run-step --workflow-name host_filter \
      --step-module idseq_dag.steps.run_czid_dedup \
      --step-class PipelineStepRunCZIDDedup \
      --step-name czid_dedup_out \
      --input-files '[["~{sep='","' select_all([reads1_fastq, reads2_fastq])}"]]' \
      --output-files '[~{if paired then '"dedup1.fastq","dedup2.fastq"' else '"dedup1.fastq"'}, "clusters.csv", "duplicate_cluster_sizes.tsv"]' \
      --output-dir-s3 '~{s3_wd_uri}' \
      --additional-files '{}' \
      --additional-attributes '{}'
  >>>
  output {
    String step_description_md = read_string("czid_dedup_out.description.md")
    File dedup1_fastq = "dedup1.fastq"
    File? dedup2_fastq = "dedup2.fastq"
    File duplicate_clusters_csv = "clusters.csv"
    File duplicate_cluster_sizes_tsv = "duplicate_cluster_sizes.tsv"
    File reads_out_count = "czid_dedup_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: 4
    memory: "16G"
  }
}

task RunSubsample {
  input {
    File reads1_fastq
    File? reads2_fastq
    File duplicate_cluster_sizes_tsv
    Int max_subsample_fragments

    String docker_image_id
    String s3_wd_uri
  }
  Boolean paired = defined(reads2_fastq)
  command<<<
  set -euxo pipefail
  TMPDIR="${TMPDIR:-/tmp}"

  # Convert FASTQs to FASTAs: the idseq-dag subsampling tool inputs and outputs FASTAs, and
  # downstream pipeline stages consume the FASTAs.
  seqtk seq -a '~{reads1_fastq}' > "$TMPDIR/reads1.fasta" & pid=$!
  fastas="\"$TMPDIR/reads1.fasta\""
  if [[ '~{paired}' == 'true' ]]; then
    seqtk seq -a '~{reads2_fastq}' > "$TMPDIR/reads2.fasta"
    wait $pid
    # also generate merged FASTA. `seqtk mergepe` interleaves the reads but doesn't append /1 /2 to
    # the names, so we add an awk kludge to do that.
    seqtk mergepe "$TMPDIR/reads1.fasta" "$TMPDIR/reads2.fasta" | awk '
        BEGIN {
          name = "";
        }
        /^>.*/ {
          if ($0 != name) {
            name = $0;
            printf("%s/1\n", $0);
          } else {
            printf("%s/2\n", $0);
          }
        }
        ! /^>.*/ { print; }
      ' > "$TMPDIR/reads_merged.fasta"
    fastas="\"$TMPDIR/reads1.fasta\",\"$TMPDIR/reads2.fasta\",\"$TMPDIR/reads_merged.fasta\""
  else
    wait $pid
  fi

  # subsample FASTAs
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_subsample \
    --step-class PipelineStepRunSubsample \
    --step-name subsampled_out \
    --input-files '[['"$fastas"'], ["~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if paired then '"subsampled_1.fa", "subsampled_2.fa", "subsampled_merged.fa"' else '"subsampled_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"max_fragments": ~{max_subsample_fragments}}'
  >>>
  output {
    String step_description_md = read_string("subsampled_out.description.md")
    File subsampled_1_fa = "subsampled_1.fa"
    File? subsampled_2_fa = "subsampled_2.fa"
    File? subsampled_merged_fa = "subsampled_merged.fa"
    File reads_out_count = "subsampled_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunStarDownstream {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] subsampled_fa
    File validate_input_summary_json
    Array[File] valid_input_fastq
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    File human_star_genome
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_star_downstream \
    --step-class PipelineStepRunStarDownstream \
    --step-name star_human_out \
    --input-files '[["~{sep='","' subsampled_fa}"], ["~{validate_input_summary_json}", "~{sep='","' valid_input_fastq}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"unmapped_human_1.fa", "unmapped_human_2.fa"' else '"unmapped_human_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"star_genome": "~{human_star_genome}"}' \
    --additional-attributes '{}'
  STAR --version > star_human_version.txt
  >>>
  output {
    String step_description_md = read_string("star_human_out.description.md")
    File unmapped_human_1_fa = "unmapped_human_1.fa"
    File? unmapped_human_2_fa = "unmapped_human_2.fa"
    File? output_read_count = "star_human_out.count"
    File? version = "star_human_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunBowtie2_bowtie2_human_out {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] unmapped_human_fa
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    File human_bowtie2_genome
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_bowtie2 \
    --step-class PipelineStepRunBowtie2 \
    --step-name bowtie2_human_out \
    --input-files '[["~{sep='","' unmapped_human_fa}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"bowtie2_human_1.fa", "bowtie2_human_2.fa", "bowtie2_human_merged.fa"' else '"bowtie2_human_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"bowtie2_genome": "~{human_bowtie2_genome}"}' \
    --additional-attributes '{"output_sam_file": "bowtie2_human.sam"}'
  bowtie2 --version > bowtie2_human_version.txt
  >>>
  output {
    String step_description_md = read_string("bowtie2_human_out.description.md")
    File bowtie2_human_1_fa = "bowtie2_human_1.fa"
    File? bowtie2_human_2_fa = "bowtie2_human_2.fa"
    File? bowtie2_human_merged_fa = "bowtie2_human_merged.fa"
    File? output_read_count = "bowtie2_human_out.count"
    File? version = "bowtie2_human_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunGsnapFilter {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] subsampled_fa
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    File gsnap_genome
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_gsnap_filter \
    --step-class PipelineStepRunGsnapFilter \
    --step-name gsnap_filter_out \
    --input-files '[["~{sep='","' subsampled_fa}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"gsnap_filter_1.fa", "gsnap_filter_2.fa", "gsnap_filter_merged.fa"' else '"gsnap_filter_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"gsnap_genome": "~{gsnap_genome}"}' \
    --additional-attributes '{"output_sam_file": "gsnap_filter.sam"}'
  gsnap --version > gsnap_filter_version.txt
  >>>
  output {
    String step_description_md = read_string("gsnap_filter_out.description.md")
    File gsnap_filter_1_fa = "gsnap_filter_1.fa"
    File? gsnap_filter_2_fa = "gsnap_filter_2.fa"
    File? gsnap_filter_merged_fa = "gsnap_filter_merged.fa"
    File? output_read_count = "gsnap_filter_out.count"
    File? version = "gsnap_filter_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}


workflow czid_host_filter {
  input {
    String docker_image_id
    String s3_wd_uri
    File fastqs_0
    File? fastqs_1
    String file_ext
    String nucleotide_type
    String host_genome
    File adapter_fasta
    File star_genome
    File bowtie2_genome
    File gsnap_genome = "s3://czid-public-references/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_pantro5_k16.tar"
    String human_star_genome
    String human_bowtie2_genome
    Int max_input_fragments
    Int max_subsample_fragments
  }

  call RunValidateInput {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      fastqs = select_all([fastqs_0, fastqs_1]),
      file_ext = file_ext,
      max_input_fragments = max_input_fragments
  }

  call RunStar {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      validate_input_summary_json = RunValidateInput.validate_input_summary_json,
      valid_input_fastq = select_all([RunValidateInput.valid_input1_fastq, RunValidateInput.valid_input2_fastq]),
      star_genome = star_genome,
      nucleotide_type = nucleotide_type,
      host_genome = host_genome
  }

  call RunTrimmomatic {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      unmapped_fastq = select_all([RunStar.unmapped1_fastq, RunStar.unmapped2_fastq]),
      adapter_fasta = adapter_fasta
  }

  call RunPriceSeq {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      trimmomatic_fastq = select_all([RunTrimmomatic.trimmomatic1_fastq, RunTrimmomatic.trimmomatic2_fastq])
  }

  call RunCZIDDedup {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      priceseq_fa = select_all([RunPriceSeq.priceseq1_fa, RunPriceSeq.priceseq2_fa])
  }

  call RunLZW {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv
  }

  call RunBowtie2_bowtie2_out {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      lzw_fa = select_all([RunLZW.lzw1_fa, RunLZW.lzw2_fa]),
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
      bowtie2_genome = bowtie2_genome
  }

  call RunSubsample {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      bowtie2_fa = select_all([RunBowtie2_bowtie2_out.bowtie2_1_fa, RunBowtie2_bowtie2_out.bowtie2_2_fa, RunBowtie2_bowtie2_out.bowtie2_merged_fa]),
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
      max_subsample_fragments = max_subsample_fragments
  }

  if (host_genome != "human") {
    call RunStarDownstream {
      input:
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri,
        subsampled_fa = select_all([RunSubsample.subsampled_1_fa, RunSubsample.subsampled_2_fa, RunSubsample.subsampled_merged_fa]),
        validate_input_summary_json = RunValidateInput.validate_input_summary_json,
        valid_input_fastq = select_all([RunValidateInput.valid_input1_fastq, RunValidateInput.valid_input2_fastq]),
        dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
        duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
        duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
        human_star_genome = human_star_genome
    }

    call RunBowtie2_bowtie2_human_out {
      input:
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri,
        unmapped_human_fa = select_all([RunStarDownstream.unmapped_human_1_fa, RunStarDownstream.unmapped_human_2_fa]),
        dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
        duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
        duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
        human_bowtie2_genome = human_bowtie2_genome
    }
  }

  Array[File] gsnap_filter_input = if (host_genome == "human")
    then select_all([RunSubsample.subsampled_1_fa, RunSubsample.subsampled_2_fa, RunSubsample.subsampled_merged_fa])
    else select_all([RunBowtie2_bowtie2_human_out.bowtie2_human_1_fa, RunBowtie2_bowtie2_human_out.bowtie2_human_2_fa, RunBowtie2_bowtie2_human_out.bowtie2_human_merged_fa])

  call RunGsnapFilter {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      subsampled_fa = gsnap_filter_input,
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
      gsnap_genome = gsnap_genome
  }

  output {
    File validate_input_out_validate_input_summary_json = RunValidateInput.validate_input_summary_json
    File? validate_input_out_count = RunValidateInput.output_read_count
    File star_out_unmapped1_fastq = RunStar.unmapped1_fastq
    File? star_out_unmapped2_fastq = RunStar.unmapped2_fastq
    File? star_out_log_file = RunStar.output_log_file
    File? star_out_count = RunStar.output_read_count
    File? star_version = RunStar.version
    File trimmomatic_out_trimmomatic1_fastq = RunTrimmomatic.trimmomatic1_fastq
    File? trimmomatic_out_trimmomatic2_fastq = RunTrimmomatic.trimmomatic2_fastq
    File? trimmomatic_out_count = RunTrimmomatic.output_read_count
    File? trimmomatic_version = RunTrimmomatic.version
    File priceseq_out_priceseq1_fa = RunPriceSeq.priceseq1_fa
    File? priceseq_out_priceseq2_fa = RunPriceSeq.priceseq2_fa
    File? priceseq_out_count = RunPriceSeq.output_read_count
    File? priceseq_version = RunPriceSeq.version
    File czid_dedup_out_dedup1_fa = RunCZIDDedup.dedup1_fa
    File? czid_dedup_out_dedup2_fa = RunCZIDDedup.dedup2_fa
    File czid_dedup_out_duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv
    File czid_dedup_out_duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv
    File? czid_dedup_out_count = RunCZIDDedup.output_read_count
    File? czid_dedup_version = RunCZIDDedup.version
    File lzw_out_lzw1_fa = RunLZW.lzw1_fa
    File? lzw_out_lzw2_fa = RunLZW.lzw2_fa
    File? lzw_out_count = RunLZW.output_read_count
    File bowtie2_out_bowtie2_1_fa = RunBowtie2_bowtie2_out.bowtie2_1_fa
    File? bowtie2_out_bowtie2_2_fa = RunBowtie2_bowtie2_out.bowtie2_2_fa
    File? bowtie2_out_bowtie2_merged_fa = RunBowtie2_bowtie2_out.bowtie2_merged_fa
    File? bowtie2_out_count = RunBowtie2_bowtie2_out.output_read_count
    File? bowtie2_version = RunBowtie2_bowtie2_out.version
    File subsampled_out_subsampled_1_fa = RunSubsample.subsampled_1_fa
    File? subsampled_out_subsampled_2_fa = RunSubsample.subsampled_2_fa
    File? subsampled_out_subsampled_merged_fa = RunSubsample.subsampled_merged_fa
    File? subsampled_out_count = RunSubsample.output_read_count
    File? star_human_out_unmapped_human_1_fa = RunStarDownstream.unmapped_human_1_fa
    File? star_human_out_unmapped_human_2_fa = RunStarDownstream.unmapped_human_2_fa
    File? star_human_out_count = RunStarDownstream.output_read_count
    File? star_human_version = RunStarDownstream.version
    File? bowtie2_human_out_bowtie2_human_1_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_1_fa
    File? bowtie2_human_out_bowtie2_human_2_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_2_fa
    File? bowtie2_human_out_bowtie2_human_merged_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_merged_fa
    File? bowtie2_human_out_count = RunBowtie2_bowtie2_human_out.output_read_count
    File? bowtie2_human_version = RunBowtie2_bowtie2_human_out.version
    File gsnap_filter_out_gsnap_filter_1_fa = RunGsnapFilter.gsnap_filter_1_fa
    File? gsnap_filter_out_gsnap_filter_2_fa = RunGsnapFilter.gsnap_filter_2_fa
    File? gsnap_filter_out_gsnap_filter_merged_fa = RunGsnapFilter.gsnap_filter_merged_fa
    File? gsnap_filter_out_count = RunGsnapFilter.output_read_count
    File? gsnap_filter_version = RunGsnapFilter.version
    File? input_read_count = RunValidateInput.input_read_count
    File? output_gene_file = RunStar.output_gene_file
    File? output_metrics_file = RunStar.output_metrics_file
    File? output_histogram_file = RunStar.output_histogram_file
  }
}
