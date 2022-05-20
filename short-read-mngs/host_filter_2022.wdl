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

task fastp {
    # fastp all-in-one for
    # - adapter trimming
    # - quality filtering
    # - complexity filtering
  input {
    String docker_image_id
    File reads1_fastq
    File? reads2_fastq
    File adapter_fasta

    Int cpu = 8
  }
  command<<<
    set -euxo pipefail
    fastp \
        -i ~{reads1_fastq} ~{"-I " + reads2_fastq} \
        -o fastp1.fastq ~{if (defined(reads2_fastq)) then "-O fastp2.fastq" else ""} \
        -w ~{cpu} \
        --dont_eval_duplication --length_required 35 \
        --qualified_quality_phred 17 --unqualified_percent_limit 15 --n_base_limit 15 \
        --low_complexity_filter --complexity_threshold 30 \
        --adapter_fasta ~{adapter_fasta} ~{if (defined(reads2_fastq)) then "--detect_adapter_for_pe" else ""}
    # TODO: double the count & make sure this works right for non-paired
    jq .read1_after_filtering.total_reads fastp.json > fastp_out.count
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

task RunCZIDDedup {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] priceseq_fa
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_czid_dedup \
    --step-class PipelineStepRunCZIDDedup \
    --step-name czid_dedup_out \
    --input-files '[["~{sep='","' priceseq_fa}"]]' \
    --output-files '[~{if length(priceseq_fa) == 2 then '"dedup1.fa", "dedup2.fa"' else '"dedup1.fa"'}, "clusters.csv", "duplicate_cluster_sizes.tsv"]' \
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

  call fastp {
      input:
          docker_image_id = docker_image_id,
          reads1_fastq = RunValidateInput.valid_input1_fastq,
          reads2_fastq = RunValidateInput.valid_input2_fastq,
          adapter_fasta = adapter_fasta
  }

  output {
    File validate_input_out_validate_input_summary_json = RunValidateInput.validate_input_summary_json
    File? validate_input_out_count = RunValidateInput.output_read_count
    File? input_read_count = RunValidateInput.input_read_count
    File fastp_out_fastp1_fastq = fastp.fastp1_fastq
    File? fastp_out_fastp2_fastq = fastp.fastp2_fastq
    File fastp_out_count = fastp.output_read_count
  }
}