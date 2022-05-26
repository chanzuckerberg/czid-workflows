version 1.0

task RunValidateInput {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] fastqs
    Int max_input_fragments
    String file_ext
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_validate_input \
    --step-class PipelineStepRunValidateInput \
    --step-name validate_input_out \
    --input-files '[["~{sep='","' fastqs}"]]' \
    --output-files '["validate_input_summary.json", ~{if length(fastqs) == 2 then '"valid_input1.fastq", "valid_input2.fastq"' else '"valid_input1.fastq"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"truncate_fragments_to": ~{max_input_fragments}, "file_ext": "~{file_ext}"}'
  >>>
  output {
    String step_description_md = read_string("validate_input_out.description.md")
    File validate_input_summary_json = "validate_input_summary.json"
    File valid_input1_fastq = "valid_input1.fastq"
    File? valid_input2_fastq = "valid_input2.fastq"
    File? output_read_count = "validate_input_out.count"
    File? input_read_count = "fastqs.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task fastp_qc {
    # fastp all-in-one for
    # - adapter trimming
    # - quality filtering
    # - complexity filtering
  input {
    String docker_image_id
    File reads1_fastq
    File? reads2_fastq
    File adapter_fasta

    Int cpu = 16
  }
  Boolean paired = defined(reads2_fastq)
  command<<<
    set -euxo pipefail
    fastp \
        -i ~{reads1_fastq} ~{"-I " + reads2_fastq} \
        -o fastp1.fastq ~{if (paired) then "-O fastp2.fastq" else ""} \
        -w ~{cpu} \
        --dont_eval_duplication --length_required 35 \
        --qualified_quality_phred 17 --unqualified_percent_limit 15 --n_base_limit 15 \
        --low_complexity_filter --complexity_threshold 30 \
        --adapter_fasta ~{adapter_fasta} ~{if (paired) then "--detect_adapter_for_pe" else ""}
    if [ '~{paired}' == 'true' ]; then
        expr 2 \* "$(jq .read1_after_filtering.total_reads fastp.json)" > fastp_out.count
    else
        jq .read1_after_filtering.total_reads fastp.json > fastp_out.count
    fi
  >>>
  output {
    #String step_description_md = read_string("fastp_out.description.md")
    File fastp1_fastq = "fastp1.fastq"
    File? fastp2_fastq = "fastp2.fastq"
    File fastp_html = "fastp.html"
    File fastp_json = "fastp.json"
    File output_read_count = "fastp_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task bowtie2_filter {
  # Remove reads [pairs] with bowtie2 hits to the given index
  input {
    File reads1_fastq
    File? reads2_fastq

    # GENOME_NAME.tar should contain GENOME_NAME/GENOME_NAME.*.bt*
    File index_tar
    String bowtie2_options = "--very-sensitive-local"

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(reads2_fastq)
  command <<<
    set -euxo pipefail
    TMPDIR="${TMPDIR:-/tmp}"

    genome_name="$(basename '~{index_tar}' .tar)"
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

    if [[ '~{paired}' == 'true' ]]; then
        samtools fastq -f 13 -1 bowtie2_filtered1.fastq -2 bowtie2_filtered2.fastq -0 /dev/null -s /dev/null "$TMPDIR/bowtie2.sam"
    else
        samtools fastq -f 4 "$TMPDIR/bowtie2.sam" > bowtie2_filtered1.fastq
    fi
  >>>

  output {
    #String step_description_md = read_string("fastp_out.description.md")
    File filtered1_fastq = "bowtie2_filtered1.fastq"
    File? filtered2_fastq = "bowtie2_filtered2.fastq"
    #File output_read_count = "fastp_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task hisat2_filter {
  # Remove reads [pairs] with HISAT2 hits to the given index
  input {
    File reads1_fastq
    File? reads2_fastq

    # GENOME_NAME.tar should contain GENOME_NAME/GENOME_NAME.*.ht2
    File index_tar
    String hisat2_options = ""

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(reads2_fastq)
  command <<<
    set -euxo pipefail
    TMPDIR="${TMPDIR:-/tmp}"

    genome_name="$(basename '~{index_tar}' .tar)"
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

    if [[ '~{paired}' == 'true' ]]; then
        samtools fastq -f 13 -1 hisat2_filtered1.fastq -2 hisat2_filtered2.fastq -0 /dev/null -s /dev/null "$TMPDIR/hisat2.sam"
    else
        samtools fastq -f 4 "$TMPDIR/hisat2.sam" > hisat2_filtered1.fastq
    fi
  >>>

  output {
    #String step_description_md = read_string("fastp_out.description.md")
    File filtered1_fastq = "hisat2_filtered1.fastq"
    File? filtered2_fastq = "hisat2_filtered2.fastq"
    #File output_read_count = "fastp_out.count"
  }
  runtime {
    docker: docker_image_id
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
    TMPDIR="${TMPDIR:-/tmp}"

    seqtk seq -a '~{reads1_fastq}' > "$TMPDIR/reads1.fa"
    if [[ '~{paired}' == 'true' ]]; then
        seqtk seq -a '~{reads2_fastq}' > "$TMPDIR/reads2.fa"
    fi

    idseq-dag-run-step --workflow-name host_filter \
      --step-module idseq_dag.steps.run_czid_dedup \
      --step-class PipelineStepRunCZIDDedup \
      --step-name czid_dedup_out \
      --input-files '[["~{sep='","' select_all([reads1_fastq, reads2_fastq])}"]]' \
      --output-files '[~{if paired then '"dedup1.fa", "dedup2.fa"' else '"dedup1.fa"'}, "clusters.csv", "duplicate_cluster_sizes.tsv"]' \
      --output-dir-s3 '~{s3_wd_uri}' \
      --additional-files '{}' \
      --additional-attributes '{}'
  >>>
  output {
    String step_description_md = read_string("czid_dedup_out.description.md")
    File dedup1_fa = "dedup1.fa"
    File? dedup2_fa = "dedup2.fa"
    File duplicate_clusters_csv = "clusters.csv"
    File duplicate_cluster_sizes_tsv = "duplicate_cluster_sizes.tsv"
    File? output_read_count = "czid_dedup_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunSubsample {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] bowtie2_fa
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    Int max_subsample_fragments
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_subsample \
    --step-class PipelineStepRunSubsample \
    --step-name subsampled_out \
    --input-files '[["~{sep='","' bowtie2_fa}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"subsampled_1.fa", "subsampled_2.fa", "subsampled_merged.fa"' else '"subsampled_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"max_fragments": ~{max_subsample_fragments}}'
  >>>
  output {
    String step_description_md = read_string("subsampled_out.description.md")
    File subsampled_1_fa = "subsampled_1.fa"
    File? subsampled_2_fa = "subsampled_2.fa"
    File? subsampled_merged_fa = "subsampled_merged.fa"
    File? output_read_count = "subsampled_out.count"
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
    #String nucleotide_type
    #String host_genome
    File adapter_fasta
    #File star_genome
    File bowtie2_index_tar
    File hisat2_index_tar
    #File gsnap_genome = "s3://czid-public-references/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_pantro5_k16.tar"
    #String human_star_genome
    #String human_bowtie2_genome
    Int max_input_fragments
    #Int max_subsample_fragments

    Int cpu = 16
  }

  call RunValidateInput {
    input:
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri,
    fastqs = select_all([fastqs_0, fastqs_1]),
    file_ext = file_ext,
    max_input_fragments = max_input_fragments
  }

  call fastp_qc {
    input:
    docker_image_id = docker_image_id,
    reads1_fastq = RunValidateInput.valid_input1_fastq,
    reads2_fastq = RunValidateInput.valid_input2_fastq,
    adapter_fasta = adapter_fasta,
    cpu = cpu
  }

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

  call RunCZIDDedup {
    input:
    reads1_fastq = hisat2_filter.filtered1_fastq,
    reads2_fastq = hisat2_filter.filtered2_fastq,
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri
  }

  output {
    File validate_input_out_validate_input_summary_json = RunValidateInput.validate_input_summary_json
    File? validate_input_out_count = RunValidateInput.output_read_count
    File? input_read_count = RunValidateInput.input_read_count
    File fastp_out_fastp1_fastq = fastp_qc.fastp1_fastq
    File? fastp_out_fastp2_fastq = fastp_qc.fastp2_fastq
    File fastp_out_count = fastp_qc.output_read_count
    File bowtie2_filtered1_fastq = bowtie2_filter.filtered1_fastq
    File? bowtie2_filtered2_fastq = bowtie2_filter.filtered2_fastq
    File hisat2_filtered1_fastq = hisat2_filter.filtered1_fastq
    File? hisat2_filtered2_fastq = hisat2_filter.filtered2_fastq
    File czid_dedup_out_dedup1_fa = RunCZIDDedup.dedup1_fa
    File? czid_dedup_out_dedup2_fa = RunCZIDDedup.dedup2_fa
    File czid_dedup_out_duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv
    File czid_dedup_out_duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv
    File? czid_dedup_out_count = RunCZIDDedup.output_read_count
  }
}