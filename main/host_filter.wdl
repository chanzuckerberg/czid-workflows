version 1.0

task RunValidateInput {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File fastqs_0
    File? fastqs_1
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_validate_input \
    --step-class PipelineStepRunValidateInput \
    --step-name validate_input_out \
    --input-files '[["~{fastqs_0}", "~{fastqs_1}"]]' \
    --output-files '["validate_input_summary.json", "valid_input1.fastq", "valid_input2.fastq"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"truncate_fragments_to": 75000000, "file_ext": "fastq"}'
  >>>
  output {
    File validate_input_summary_json = "validate_input_summary.json"
    File valid_input1_fastq = "valid_input1.fastq"
    File valid_input2_fastq = "valid_input2.fastq"
    File? output_read_count = "validate_input_out.count"
    File? input_read_count = "fastqs.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunStar {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File validate_input_summary_json
    File valid_input1_fastq
    File valid_input2_fastq
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_star \
    --step-class PipelineStepRunStar \
    --step-name star_out \
    --input-files '[["~{validate_input_summary_json}", "~{valid_input1_fastq}", "~{valid_input2_fastq}"]]' \
    --output-files '["unmapped1.fastq", "unmapped2.fastq"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"star_genome": "s3://idseq-database/host_filter/mosquitos/2019-10-02/mosquito_STAR_genome.tar"}' \
    --additional-attributes '{"output_gene_file": "reads_per_gene.star.tab", "nucleotide_type": "DNA", "host_genome": "mosquito", "output_metrics_file": "picard_insert_metrics.txt", "output_histogram_file": "insert_size_histogram.pdf"}'
  >>>
  output {
    File unmapped1_fastq = "unmapped1.fastq"
    File unmapped2_fastq = "unmapped2.fastq"
    File? output_read_count = "star_out.count"
    File? output_gene_file = "reads_per_gene.star.tab"
    File? output_metrics_file = "picard_insert_metrics.txt"
    File? output_histogram_file = "insert_size_histogram.pdf"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunTrimmomatic {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File unmapped1_fastq
    File unmapped2_fastq
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_trimmomatic \
    --step-class PipelineStepRunTrimmomatic \
    --step-name trimmomatic_out \
    --input-files '[["~{unmapped1_fastq}", "~{unmapped2_fastq}"]]' \
    --output-files '["trimmomatic1.fastq", "trimmomatic2.fastq"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"adapter_fasta": "s3://idseq-database/adapter_sequences/illumina_TruSeq3-PE-2_NexteraPE-PE.fasta"}' \
    --additional-attributes '{}'
  >>>
  output {
    File trimmomatic1_fastq = "trimmomatic1.fastq"
    File trimmomatic2_fastq = "trimmomatic2.fastq"
    File? output_read_count = "trimmomatic_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunPriceSeq {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File trimmomatic1_fastq
    File trimmomatic2_fastq
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_priceseq \
    --step-class PipelineStepRunPriceSeq \
    --step-name priceseq_out \
    --input-files '[["~{trimmomatic1_fastq}", "~{trimmomatic2_fastq}"]]' \
    --output-files '["priceseq1.fa", "priceseq2.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File priceseq1_fa = "priceseq1.fa"
    File priceseq2_fa = "priceseq2.fa"
    File? output_read_count = "priceseq_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunCDHitDup {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File priceseq1_fa
    File priceseq2_fa
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_cdhitdup \
    --step-class PipelineStepRunCDHitDup \
    --step-name cdhitdup_out \
    --input-files '[["~{priceseq1_fa}", "~{priceseq2_fa}"]]' \
    --output-files '["dedup1.fa", "dedup2.fa", "dedup1.fa.clstr", "cdhitdup_cluster_sizes.tsv"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File dedup1_fa = "dedup1.fa"
    File dedup2_fa = "dedup2.fa"
    File dedup1_fa_clstr = "dedup1.fa.clstr"
    File cdhitdup_cluster_sizes_tsv = "cdhitdup_cluster_sizes.tsv"
    File? output_read_count = "cdhitdup_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunLZW {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File dedup1_fa
    File dedup2_fa
    File dedup1_fa_clstr
    File cdhitdup_cluster_sizes_tsv
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_lzw \
    --step-class PipelineStepRunLZW \
    --step-name lzw_out \
    --input-files '[["~{dedup1_fa}", "~{dedup2_fa}", "~{dedup1_fa_clstr}", "~{cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["lzw1.fa", "lzw2.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"thresholds": [0.45, 0.42], "threshold_readlength": 150}'
  >>>
  output {
    File lzw1_fa = "lzw1.fa"
    File lzw2_fa = "lzw2.fa"
    File? output_read_count = "lzw_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunBowtie2_bowtie2_out {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File lzw1_fa
    File lzw2_fa
    File dedup1_fa
    File dedup2_fa
    File dedup1_fa_clstr
    File cdhitdup_cluster_sizes_tsv
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_bowtie2 \
    --step-class PipelineStepRunBowtie2 \
    --step-name bowtie2_out \
    --input-files '[["~{lzw1_fa}", "~{lzw2_fa}"], ["~{dedup1_fa}", "~{dedup2_fa}", "~{dedup1_fa_clstr}", "~{cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["bowtie2_1.fa", "bowtie2_2.fa", "bowtie2_merged.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"bowtie2_genome": "s3://idseq-database/host_filter/mosquitos/2019-10-02/mosquito_bowtie2_genome.tar"}' \
    --additional-attributes '{"output_sam_file": "bowtie2.sam"}'
  >>>
  output {
    File bowtie2_1_fa = "bowtie2_1.fa"
    File bowtie2_2_fa = "bowtie2_2.fa"
    File bowtie2_merged_fa = "bowtie2_merged.fa"
    File? output_read_count = "bowtie2_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunSubsample {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File bowtie2_1_fa
    File bowtie2_2_fa
    File bowtie2_merged_fa
    File dedup1_fa
    File dedup2_fa
    File dedup1_fa_clstr
    File cdhitdup_cluster_sizes_tsv
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_subsample \
    --step-class PipelineStepRunSubsample \
    --step-name subsampled_out \
    --input-files '[["~{bowtie2_1_fa}", "~{bowtie2_2_fa}", "~{bowtie2_merged_fa}"], ["~{dedup1_fa}", "~{dedup2_fa}", "~{dedup1_fa_clstr}", "~{cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["subsampled_1.fa", "subsampled_2.fa", "subsampled_merged.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"max_fragments": 1000000}'
  >>>
  output {
    File subsampled_1_fa = "subsampled_1.fa"
    File subsampled_2_fa = "subsampled_2.fa"
    File subsampled_merged_fa = "subsampled_merged.fa"
    File? output_read_count = "subsampled_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunStarDownstream {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File subsampled_1_fa
    File subsampled_2_fa
    File subsampled_merged_fa
    File validate_input_summary_json
    File valid_input1_fastq
    File valid_input2_fastq
    File dedup1_fa
    File dedup2_fa
    File dedup1_fa_clstr
    File cdhitdup_cluster_sizes_tsv
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_star_downstream \
    --step-class PipelineStepRunStarDownstream \
    --step-name star_human_out \
    --input-files '[["~{subsampled_1_fa}", "~{subsampled_2_fa}", "~{subsampled_merged_fa}"], ["~{validate_input_summary_json}", "~{valid_input1_fastq}", "~{valid_input2_fastq}"], ["~{dedup1_fa}", "~{dedup2_fa}", "~{dedup1_fa_clstr}", "~{cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["unmapped_human_1.fa", "unmapped_human_2.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"star_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/human_STAR_genome.tar"}' \
    --additional-attributes '{}'
  >>>
  output {
    File unmapped_human_1_fa = "unmapped_human_1.fa"
    File unmapped_human_2_fa = "unmapped_human_2.fa"
    File? output_read_count = "star_human_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunBowtie2_bowtie2_human_out {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File unmapped_human_1_fa
    File unmapped_human_2_fa
    File dedup1_fa
    File dedup2_fa
    File dedup1_fa_clstr
    File cdhitdup_cluster_sizes_tsv
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_bowtie2 \
    --step-class PipelineStepRunBowtie2 \
    --step-name bowtie2_human_out \
    --input-files '[["~{unmapped_human_1_fa}", "~{unmapped_human_2_fa}"], ["~{dedup1_fa}", "~{dedup2_fa}", "~{dedup1_fa_clstr}", "~{cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["bowtie2_human_1.fa", "bowtie2_human_2.fa", "bowtie2_human_merged.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"bowtie2_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/human_bowtie2_genome.tar"}' \
    --additional-attributes '{"output_sam_file": "bowtie2_human.sam"}'
  >>>
  output {
    File bowtie2_human_1_fa = "bowtie2_human_1.fa"
    File bowtie2_human_2_fa = "bowtie2_human_2.fa"
    File bowtie2_human_merged_fa = "bowtie2_human_merged.fa"
    File? output_read_count = "bowtie2_human_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunGsnapFilter {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File bowtie2_human_1_fa
    File bowtie2_human_2_fa
    File bowtie2_human_merged_fa
    File dedup1_fa
    File dedup2_fa
    File dedup1_fa_clstr
    File cdhitdup_cluster_sizes_tsv
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_gsnap_filter \
    --step-class PipelineStepRunGsnapFilter \
    --step-name gsnap_filter_out \
    --input-files '[["~{bowtie2_human_1_fa}", "~{bowtie2_human_2_fa}", "~{bowtie2_human_merged_fa}"], ["~{dedup1_fa}", "~{dedup2_fa}", "~{dedup1_fa_clstr}", "~{cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["gsnap_filter_1.fa", "gsnap_filter_2.fa", "gsnap_filter_merged.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"gsnap_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_pantro5_k16.tar"}' \
    --additional-attributes '{"output_sam_file": "gsnap_filter.sam"}'
  >>>
  output {
    File gsnap_filter_1_fa = "gsnap_filter_1.fa"
    File gsnap_filter_2_fa = "gsnap_filter_2.fa"
    File gsnap_filter_merged_fa = "gsnap_filter_merged.fa"
    File? output_read_count = "gsnap_filter_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

workflow idseq_host_filter {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File fastqs_0
    File? fastqs_1
  }

  call RunValidateInput {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      fastqs_0 = fastqs_0,
      fastqs_1 = fastqs_1
  }

  call RunStar {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      validate_input_summary_json = RunValidateInput.validate_input_summary_json,
      valid_input1_fastq = RunValidateInput.valid_input1_fastq,
      valid_input2_fastq = RunValidateInput.valid_input2_fastq
  }

  call RunTrimmomatic {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      unmapped1_fastq = RunStar.unmapped1_fastq,
      unmapped2_fastq = RunStar.unmapped2_fastq
  }

  call RunPriceSeq {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      trimmomatic1_fastq = RunTrimmomatic.trimmomatic1_fastq,
      trimmomatic2_fastq = RunTrimmomatic.trimmomatic2_fastq
  }

  call RunCDHitDup {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      priceseq1_fa = RunPriceSeq.priceseq1_fa,
      priceseq2_fa = RunPriceSeq.priceseq2_fa
  }

  call RunLZW {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      dedup1_fa = RunCDHitDup.dedup1_fa,
      dedup2_fa = RunCDHitDup.dedup2_fa,
      dedup1_fa_clstr = RunCDHitDup.dedup1_fa_clstr,
      cdhitdup_cluster_sizes_tsv = RunCDHitDup.cdhitdup_cluster_sizes_tsv
  }

  call RunBowtie2_bowtie2_out {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      lzw1_fa = RunLZW.lzw1_fa,
      lzw2_fa = RunLZW.lzw2_fa,
      dedup1_fa = RunCDHitDup.dedup1_fa,
      dedup2_fa = RunCDHitDup.dedup2_fa,
      dedup1_fa_clstr = RunCDHitDup.dedup1_fa_clstr,
      cdhitdup_cluster_sizes_tsv = RunCDHitDup.cdhitdup_cluster_sizes_tsv
  }

  call RunSubsample {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      bowtie2_1_fa = RunBowtie2_bowtie2_out.bowtie2_1_fa,
      bowtie2_2_fa = RunBowtie2_bowtie2_out.bowtie2_2_fa,
      bowtie2_merged_fa = RunBowtie2_bowtie2_out.bowtie2_merged_fa,
      dedup1_fa = RunCDHitDup.dedup1_fa,
      dedup2_fa = RunCDHitDup.dedup2_fa,
      dedup1_fa_clstr = RunCDHitDup.dedup1_fa_clstr,
      cdhitdup_cluster_sizes_tsv = RunCDHitDup.cdhitdup_cluster_sizes_tsv
  }

  call RunStarDownstream {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      subsampled_1_fa = RunSubsample.subsampled_1_fa,
      subsampled_2_fa = RunSubsample.subsampled_2_fa,
      subsampled_merged_fa = RunSubsample.subsampled_merged_fa,
      validate_input_summary_json = RunValidateInput.validate_input_summary_json,
      valid_input1_fastq = RunValidateInput.valid_input1_fastq,
      valid_input2_fastq = RunValidateInput.valid_input2_fastq,
      dedup1_fa = RunCDHitDup.dedup1_fa,
      dedup2_fa = RunCDHitDup.dedup2_fa,
      dedup1_fa_clstr = RunCDHitDup.dedup1_fa_clstr,
      cdhitdup_cluster_sizes_tsv = RunCDHitDup.cdhitdup_cluster_sizes_tsv
  }

  call RunBowtie2_bowtie2_human_out {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      unmapped_human_1_fa = RunStarDownstream.unmapped_human_1_fa,
      unmapped_human_2_fa = RunStarDownstream.unmapped_human_2_fa,
      dedup1_fa = RunCDHitDup.dedup1_fa,
      dedup2_fa = RunCDHitDup.dedup2_fa,
      dedup1_fa_clstr = RunCDHitDup.dedup1_fa_clstr,
      cdhitdup_cluster_sizes_tsv = RunCDHitDup.cdhitdup_cluster_sizes_tsv
  }

  call RunGsnapFilter {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      bowtie2_human_1_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_1_fa,
      bowtie2_human_2_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_2_fa,
      bowtie2_human_merged_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_merged_fa,
      dedup1_fa = RunCDHitDup.dedup1_fa,
      dedup2_fa = RunCDHitDup.dedup2_fa,
      dedup1_fa_clstr = RunCDHitDup.dedup1_fa_clstr,
      cdhitdup_cluster_sizes_tsv = RunCDHitDup.cdhitdup_cluster_sizes_tsv
  }

  output {
    File validate_input_out_validate_input_summary_json = RunValidateInput.validate_input_summary_json
    File validate_input_out_valid_input1_fastq = RunValidateInput.valid_input1_fastq
    File validate_input_out_valid_input2_fastq = RunValidateInput.valid_input2_fastq
    File? validate_input_out_count = RunValidateInput.output_read_count
    File star_out_unmapped1_fastq = RunStar.unmapped1_fastq
    File star_out_unmapped2_fastq = RunStar.unmapped2_fastq
    File? star_out_count = RunStar.output_read_count
    File trimmomatic_out_trimmomatic1_fastq = RunTrimmomatic.trimmomatic1_fastq
    File trimmomatic_out_trimmomatic2_fastq = RunTrimmomatic.trimmomatic2_fastq
    File? trimmomatic_out_count = RunTrimmomatic.output_read_count
    File priceseq_out_priceseq1_fa = RunPriceSeq.priceseq1_fa
    File priceseq_out_priceseq2_fa = RunPriceSeq.priceseq2_fa
    File? priceseq_out_count = RunPriceSeq.output_read_count
    File cdhitdup_out_dedup1_fa = RunCDHitDup.dedup1_fa
    File cdhitdup_out_dedup2_fa = RunCDHitDup.dedup2_fa
    File cdhitdup_out_dedup1_fa_clstr = RunCDHitDup.dedup1_fa_clstr
    File cdhitdup_out_cdhitdup_cluster_sizes_tsv = RunCDHitDup.cdhitdup_cluster_sizes_tsv
    File? cdhitdup_out_count = RunCDHitDup.output_read_count
    File lzw_out_lzw1_fa = RunLZW.lzw1_fa
    File lzw_out_lzw2_fa = RunLZW.lzw2_fa
    File? lzw_out_count = RunLZW.output_read_count
    File bowtie2_out_bowtie2_1_fa = RunBowtie2_bowtie2_out.bowtie2_1_fa
    File bowtie2_out_bowtie2_2_fa = RunBowtie2_bowtie2_out.bowtie2_2_fa
    File bowtie2_out_bowtie2_merged_fa = RunBowtie2_bowtie2_out.bowtie2_merged_fa
    File? bowtie2_out_count = RunBowtie2_bowtie2_out.output_read_count
    File subsampled_out_subsampled_1_fa = RunSubsample.subsampled_1_fa
    File subsampled_out_subsampled_2_fa = RunSubsample.subsampled_2_fa
    File subsampled_out_subsampled_merged_fa = RunSubsample.subsampled_merged_fa
    File? subsampled_out_count = RunSubsample.output_read_count
    File star_human_out_unmapped_human_1_fa = RunStarDownstream.unmapped_human_1_fa
    File star_human_out_unmapped_human_2_fa = RunStarDownstream.unmapped_human_2_fa
    File? star_human_out_count = RunStarDownstream.output_read_count
    File bowtie2_human_out_bowtie2_human_1_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_1_fa
    File bowtie2_human_out_bowtie2_human_2_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_2_fa
    File bowtie2_human_out_bowtie2_human_merged_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_merged_fa
    File? bowtie2_human_out_count = RunBowtie2_bowtie2_human_out.output_read_count
    File gsnap_filter_out_gsnap_filter_1_fa = RunGsnapFilter.gsnap_filter_1_fa
    File gsnap_filter_out_gsnap_filter_2_fa = RunGsnapFilter.gsnap_filter_2_fa
    File gsnap_filter_out_gsnap_filter_merged_fa = RunGsnapFilter.gsnap_filter_merged_fa
    File? gsnap_filter_out_count = RunGsnapFilter.output_read_count
    File? input_read_count = RunValidateInput.input_read_count
    File? output_gene_file = RunStar.output_gene_file
    File? output_metrics_file = RunStar.output_metrics_file
    File? output_histogram_file = RunStar.output_histogram_file
  }
}
