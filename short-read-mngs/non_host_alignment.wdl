version 1.0

task RunAlignment_gsnap_out {
  input {
    String docker_image_id
    String s3_wd_uri
    String? genome_name
    Array[File] host_filter_out_gsnap_filter_fa
    File duplicate_cluster_sizes_tsv
    File? index
    File lineage_db
    File accession2taxid_db
    File taxon_blacklist
    File deuterostome_db
    String index_dir_suffix
    Boolean use_deuterostome_filter
    Boolean use_taxon_whitelist
    Boolean? run_locally = false
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.run_alignment \
    --step-class PipelineStepRunAlignment \
    --step-name gsnap_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '["gsnap.m8", "gsnap.deduped.m8", "gsnap.hitsummary.tab", "gsnap_counts_with_dcr.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "accession2taxid_db": "~{accession2taxid_db}", "taxon_blacklist": "~{taxon_blacklist}", "deuterostome_db": "~{if use_deuterostome_filter then '~{deuterostome_db}' else ''}" ~{if defined(index) then ', "index": "~{index}"' else ''} }' \
    --additional-attributes '{"alignment_algorithm": "gsnap", "index_dir_suffix": "~{index_dir_suffix}", "use_taxon_whitelist": ~{use_taxon_whitelist}, "run_locally": ~{run_locally} ~{if defined(genome_name) then ' , "genome_name": "~{genome_name}"' else ''} }'
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

task RunAlignment_rapsearch2_out {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] host_filter_out_gsnap_filter_fa
    File duplicate_cluster_sizes_tsv
    File lineage_db
    File accession2taxid_db
    File taxon_blacklist
    File? index
    String index_dir_suffix
    Boolean use_taxon_whitelist
	  Boolean? run_locally = false
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.run_alignment \
    --step-class PipelineStepRunAlignment \
    --step-name rapsearch2_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '["rapsearch2.m8", "rapsearch2.deduped.m8", "rapsearch2.hitsummary.tab", "rapsearch2_counts_with_dcr.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "accession2taxid_db": "~{accession2taxid_db}", "taxon_blacklist": "~{taxon_blacklist}" ~{if defined(index) then ', "index": "~{index}"' else ''} }' \
    --additional-attributes '{"alignment_algorithm": "rapsearch2", "index_dir_suffix": "~{index_dir_suffix}", "use_taxon_whitelist": ~{use_taxon_whitelist}, "run_locally": ~{run_locally} }'
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
  set -euxo pipefail
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
    File idseq_dedup_out_duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.generate_annotated_fasta \
    --step-class PipelineStepGenerateAnnotatedFasta \
    --step-name annotated_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{gsnap_m8}", "~{gsnap_deduped_m8}", "~{gsnap_hitsummary_tab}", "~{gsnap_counts_with_dcr_json}"], ["~{rapsearch2_m8}", "~{rapsearch2_deduped_m8}", "~{rapsearch2_hitsummary_tab}", "~{rapsearch2_counts_with_dcr_json}"], ["~{idseq_dedup_out_duplicate_clusters_csv}"], ["~{duplicate_cluster_sizes_tsv}"]]' \
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
    String s3_wd_uri
    File host_filter_out_gsnap_filter_1_fa
    File? host_filter_out_gsnap_filter_2_fa
    File? host_filter_out_gsnap_filter_merged_fa
    File duplicate_cluster_sizes_tsv
    File idseq_dedup_out_duplicate_clusters_csv
    String index_version = "2020-04-20"
    File lineage_db = "s3://idseq-public-references/taxonomy/2020-04-20/taxid-lineages.db"
    File accession2taxid_db = "s3://idseq-public-references/alignment_data/2020-04-20/accession2taxid.db"
    File taxon_blacklist = "s3://idseq-public-references/taxonomy/2020-04-20/taxon_blacklist.txt"
    String index_dir_suffix = index_version
    File deuterostome_db = "s3://idseq-public-references/taxonomy/2020-04-20/deuterostome_taxids.txt"
    Boolean use_deuterostome_filter = true
    Boolean use_taxon_whitelist = false
    File? local_gsnap_index
    String? local_gsnap_genome_name
    File? local_rapsearch2_index
  }

  call RunAlignment_gsnap_out {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      duplicate_cluster_sizes_tsv = duplicate_cluster_sizes_tsv,
      lineage_db = lineage_db,
      accession2taxid_db = accession2taxid_db,
      taxon_blacklist = taxon_blacklist,
      deuterostome_db = deuterostome_db,
      index_dir_suffix = index_dir_suffix,
      use_deuterostome_filter = use_deuterostome_filter,
      use_taxon_whitelist = use_taxon_whitelist,
      run_locally = defined(local_gsnap_index),
      index = local_gsnap_index,
      genome_name = local_gsnap_genome_name
  }

  call RunAlignment_rapsearch2_out {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      duplicate_cluster_sizes_tsv = duplicate_cluster_sizes_tsv,
      lineage_db = lineage_db,
      accession2taxid_db = accession2taxid_db,
      taxon_blacklist = taxon_blacklist,
      index_dir_suffix = index_dir_suffix,
      use_taxon_whitelist = use_taxon_whitelist,
      run_locally = defined(local_rapsearch2_index),
      index = local_rapsearch2_index
  }

  call CombineTaxonCounts {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      gsnap_m8 = RunAlignment_gsnap_out.gsnap_m8,
      gsnap_deduped_m8 = RunAlignment_gsnap_out.gsnap_deduped_m8,
      gsnap_hitsummary_tab = RunAlignment_gsnap_out.gsnap_hitsummary_tab,
      gsnap_counts_with_dcr_json = RunAlignment_gsnap_out.gsnap_counts_with_dcr_json,
      rapsearch2_m8 = RunAlignment_rapsearch2_out.rapsearch2_m8,
      rapsearch2_deduped_m8 = RunAlignment_rapsearch2_out.rapsearch2_deduped_m8,
      rapsearch2_hitsummary_tab = RunAlignment_rapsearch2_out.rapsearch2_hitsummary_tab,
      rapsearch2_counts_with_dcr_json = RunAlignment_rapsearch2_out.rapsearch2_counts_with_dcr_json
  }

  call GenerateAnnotatedFasta {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      gsnap_m8 = RunAlignment_gsnap_out.gsnap_m8,
      gsnap_deduped_m8 = RunAlignment_gsnap_out.gsnap_deduped_m8,
      gsnap_hitsummary_tab = RunAlignment_gsnap_out.gsnap_hitsummary_tab,
      gsnap_counts_with_dcr_json = RunAlignment_gsnap_out.gsnap_counts_with_dcr_json,
      rapsearch2_m8 = RunAlignment_rapsearch2_out.rapsearch2_m8,
      rapsearch2_deduped_m8 = RunAlignment_rapsearch2_out.rapsearch2_deduped_m8,
      rapsearch2_hitsummary_tab = RunAlignment_rapsearch2_out.rapsearch2_hitsummary_tab,
      rapsearch2_counts_with_dcr_json = RunAlignment_rapsearch2_out.rapsearch2_counts_with_dcr_json,
      idseq_dedup_out_duplicate_clusters_csv = idseq_dedup_out_duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = duplicate_cluster_sizes_tsv
  }

  output {
    File gsnap_out_gsnap_m8 = RunAlignment_gsnap_out.gsnap_m8
    File gsnap_out_gsnap_deduped_m8 = RunAlignment_gsnap_out.gsnap_deduped_m8
    File gsnap_out_gsnap_hitsummary_tab = RunAlignment_gsnap_out.gsnap_hitsummary_tab
    File gsnap_out_gsnap_counts_with_dcr_json = RunAlignment_gsnap_out.gsnap_counts_with_dcr_json
    File? gsnap_out_count = RunAlignment_gsnap_out.output_read_count
    File rapsearch2_out_rapsearch2_m8 = RunAlignment_rapsearch2_out.rapsearch2_m8
    File rapsearch2_out_rapsearch2_deduped_m8 = RunAlignment_rapsearch2_out.rapsearch2_deduped_m8
    File rapsearch2_out_rapsearch2_hitsummary_tab = RunAlignment_rapsearch2_out.rapsearch2_hitsummary_tab
    File rapsearch2_out_rapsearch2_counts_with_dcr_json = RunAlignment_rapsearch2_out.rapsearch2_counts_with_dcr_json
    File? rapsearch2_out_count = RunAlignment_rapsearch2_out.output_read_count
    File taxon_count_out_taxon_counts_with_dcr_json = CombineTaxonCounts.taxon_counts_with_dcr_json
    File? taxon_count_out_count = CombineTaxonCounts.output_read_count
    File annotated_out_annotated_merged_fa = GenerateAnnotatedFasta.annotated_merged_fa
    File annotated_out_unidentified_fa = GenerateAnnotatedFasta.unidentified_fa
    File? annotated_out_count = GenerateAnnotatedFasta.output_read_count
  }
}
