version 1.0

task RunAlignmentRemotely_gsnap_out {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    Array[File] host_filter_out_gsnap_filter_fa
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
    String lineage_db
    String accession2taxid_db
    String taxon_blacklist
    String? deuterostome_db
    String index_dir_suffix
    Boolean use_taxon_whitelist
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.run_alignment_remotely \
    --step-class PipelineStepRunAlignmentRemotely \
    --step-name gsnap_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["gsnap.m8", "gsnap.deduped.m8", "gsnap.hitsummary.tab", "gsnap_counts_with_dcr.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "accession2taxid_db": "~{accession2taxid_db}", "taxon_blacklist": "~{taxon_blacklist}", "deuterostome_db": "~{deuterostome_db}"}' \
    --additional-attributes '{"alignment_algorithm": "gsnap", "index_dir_suffix": "~{index_dir_suffix}", "use_taxon_whitelist": "~{use_taxon_whitelist}"}'
  >>>
  output {
    File gsnap_m8 = "gsnap.m8"
    File gsnap_deduped_m8 = "gsnap.deduped.m8"
    File gsnap_hitsummary_tab = "gsnap.hitsummary.tab"
    File gsnap_counts_with_dcr_json = "gsnap_counts_with_dcr.json"
    File? output_read_count = "gsnap_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunAlignmentRemotely_rapsearch2_out {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    Array[File] host_filter_out_gsnap_filter_fa
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
    String lineage_db
    String accession2taxid_db
    String taxon_blacklist
    String index_dir_suffix
    Boolean use_taxon_whitelist
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.run_alignment_remotely \
    --step-class PipelineStepRunAlignmentRemotely \
    --step-name rapsearch2_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["rapsearch2.m8", "rapsearch2.deduped.m8", "rapsearch2.hitsummary.tab", "rapsearch2_counts_with_dcr.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "accession2taxid_db": "~{accession2taxid_db}", "taxon_blacklist": "~{taxon_blacklist}"}' \
    --additional-attributes '{"alignment_algorithm": "rapsearch2", "index_dir_suffix": "~{index_dir_suffix}", "use_taxon_whitelist": "~{use_taxon_whitelist}"}'
  >>>
  output {
    File rapsearch2_m8 = "rapsearch2.m8"
    File rapsearch2_deduped_m8 = "rapsearch2.deduped.m8"
    File rapsearch2_hitsummary_tab = "rapsearch2.hitsummary.tab"
    File rapsearch2_counts_with_dcr_json = "rapsearch2_counts_with_dcr.json"
    File? output_read_count = "rapsearch2_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task CombineTaxonCounts {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File gsnap_m8
    File gsnap_deduped_m8
    File gsnap_hitsummary_tab
    File gsnap_counts_with_dcr_json
    File rapsearch2_m8
    File rapsearch2_deduped_m8
    File rapsearch2_hitsummary_tab
    File rapsearch2_counts_with_dcr_json
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.combine_taxon_counts \
    --step-class PipelineStepCombineTaxonCounts \
    --step-name taxon_count_out \
    --input-files '[["~{gsnap_m8}", "~{gsnap_deduped_m8}", "~{gsnap_hitsummary_tab}", "~{gsnap_counts_with_dcr_json}"], ["~{rapsearch2_m8}", "~{rapsearch2_deduped_m8}", "~{rapsearch2_hitsummary_tab}", "~{rapsearch2_counts_with_dcr_json}"]]' \
    --output-files '["taxon_counts_with_dcr.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File taxon_counts_with_dcr_json = "taxon_counts_with_dcr.json"
    File? output_read_count = "taxon_count_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateAnnotatedFasta {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    Array[File] host_filter_out_gsnap_filter_fa
    File gsnap_m8
    File gsnap_deduped_m8
    File gsnap_hitsummary_tab
    File gsnap_counts_with_dcr_json
    File rapsearch2_m8
    File rapsearch2_deduped_m8
    File rapsearch2_hitsummary_tab
    File rapsearch2_counts_with_dcr_json
    File cdhitdup_out_dedup1_fa_clstr
    File cdhitdup_out_dedup1_fa
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.generate_annotated_fasta \
    --step-class PipelineStepGenerateAnnotatedFasta \
    --step-name annotated_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{gsnap_m8}", "~{gsnap_deduped_m8}", "~{gsnap_hitsummary_tab}", "~{gsnap_counts_with_dcr_json}"], ["~{rapsearch2_m8}", "~{rapsearch2_deduped_m8}", "~{rapsearch2_hitsummary_tab}", "~{rapsearch2_counts_with_dcr_json}"], ["~{cdhitdup_out_dedup1_fa_clstr}", "~{cdhitdup_out_dedup1_fa}"], ["~{cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["annotated_merged.fa", "unidentified.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File annotated_merged_fa = "annotated_merged.fa"
    File unidentified_fa = "unidentified.fa"
    File? output_read_count = "annotated_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

workflow idseq_non_host_alignment {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File host_filter_out_gsnap_filter_1_fa
    File? host_filter_out_gsnap_filter_2_fa
    File? host_filter_out_gsnap_filter_merged_fa
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
    File cdhitdup_out_dedup1_fa_clstr
    File cdhitdup_out_dedup1_fa
    String lineage_db
    String accession2taxid_db
    String taxon_blacklist
    String index_dir_suffix
    String? deuterostome_db
    Boolean use_taxon_whitelist
  }

  call RunAlignmentRemotely_gsnap_out {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv,
      lineage_db = lineage_db,
      accession2taxid_db = accession2taxid_db,
      taxon_blacklist = taxon_blacklist,
      deuterostome_db = deuterostome_db,
      index_dir_suffix = index_dir_suffix,
      use_taxon_whitelist = use_taxon_whitelist
  }

  call RunAlignmentRemotely_rapsearch2_out {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv,
      lineage_db = lineage_db,
      accession2taxid_db = accession2taxid_db,
      taxon_blacklist = taxon_blacklist,
      index_dir_suffix = index_dir_suffix,
      use_taxon_whitelist = use_taxon_whitelist
  }

  call CombineTaxonCounts {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      gsnap_m8 = RunAlignmentRemotely_gsnap_out.gsnap_m8,
      gsnap_deduped_m8 = RunAlignmentRemotely_gsnap_out.gsnap_deduped_m8,
      gsnap_hitsummary_tab = RunAlignmentRemotely_gsnap_out.gsnap_hitsummary_tab,
      gsnap_counts_with_dcr_json = RunAlignmentRemotely_gsnap_out.gsnap_counts_with_dcr_json,
      rapsearch2_m8 = RunAlignmentRemotely_rapsearch2_out.rapsearch2_m8,
      rapsearch2_deduped_m8 = RunAlignmentRemotely_rapsearch2_out.rapsearch2_deduped_m8,
      rapsearch2_hitsummary_tab = RunAlignmentRemotely_rapsearch2_out.rapsearch2_hitsummary_tab,
      rapsearch2_counts_with_dcr_json = RunAlignmentRemotely_rapsearch2_out.rapsearch2_counts_with_dcr_json
  }

  call GenerateAnnotatedFasta {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      gsnap_m8 = RunAlignmentRemotely_gsnap_out.gsnap_m8,
      gsnap_deduped_m8 = RunAlignmentRemotely_gsnap_out.gsnap_deduped_m8,
      gsnap_hitsummary_tab = RunAlignmentRemotely_gsnap_out.gsnap_hitsummary_tab,
      gsnap_counts_with_dcr_json = RunAlignmentRemotely_gsnap_out.gsnap_counts_with_dcr_json,
      rapsearch2_m8 = RunAlignmentRemotely_rapsearch2_out.rapsearch2_m8,
      rapsearch2_deduped_m8 = RunAlignmentRemotely_rapsearch2_out.rapsearch2_deduped_m8,
      rapsearch2_hitsummary_tab = RunAlignmentRemotely_rapsearch2_out.rapsearch2_hitsummary_tab,
      rapsearch2_counts_with_dcr_json = RunAlignmentRemotely_rapsearch2_out.rapsearch2_counts_with_dcr_json,
      cdhitdup_out_dedup1_fa_clstr = cdhitdup_out_dedup1_fa_clstr,
      cdhitdup_out_dedup1_fa = cdhitdup_out_dedup1_fa,
      cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
  }

  output {
    File gsnap_out_gsnap_m8 = RunAlignmentRemotely_gsnap_out.gsnap_m8
    File gsnap_out_gsnap_deduped_m8 = RunAlignmentRemotely_gsnap_out.gsnap_deduped_m8
    File gsnap_out_gsnap_hitsummary_tab = RunAlignmentRemotely_gsnap_out.gsnap_hitsummary_tab
    File gsnap_out_gsnap_counts_with_dcr_json = RunAlignmentRemotely_gsnap_out.gsnap_counts_with_dcr_json
    File? gsnap_out_count = RunAlignmentRemotely_gsnap_out.output_read_count
    File rapsearch2_out_rapsearch2_m8 = RunAlignmentRemotely_rapsearch2_out.rapsearch2_m8
    File rapsearch2_out_rapsearch2_deduped_m8 = RunAlignmentRemotely_rapsearch2_out.rapsearch2_deduped_m8
    File rapsearch2_out_rapsearch2_hitsummary_tab = RunAlignmentRemotely_rapsearch2_out.rapsearch2_hitsummary_tab
    File rapsearch2_out_rapsearch2_counts_with_dcr_json = RunAlignmentRemotely_rapsearch2_out.rapsearch2_counts_with_dcr_json
    File? rapsearch2_out_count = RunAlignmentRemotely_rapsearch2_out.output_read_count
    File taxon_count_out_taxon_counts_with_dcr_json = CombineTaxonCounts.taxon_counts_with_dcr_json
    File? taxon_count_out_count = CombineTaxonCounts.output_read_count
    File annotated_out_annotated_merged_fa = GenerateAnnotatedFasta.annotated_merged_fa
    File annotated_out_unidentified_fa = GenerateAnnotatedFasta.unidentified_fa
    File? annotated_out_count = GenerateAnnotatedFasta.output_read_count
  }
}
