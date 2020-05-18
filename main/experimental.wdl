version 1.0

task GenerateTaxidFasta {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File taxid_fasta_in_annotated_merged_fa
    File taxid_fasta_in_gsnap_hitsummary_tab
    File taxid_fasta_in_rapsearch2_hitsummary_tab
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.generate_taxid_fasta \
    --step-class PipelineStepGenerateTaxidFasta \
    --step-name taxid_fasta_out \
    --input-files '[["~{taxid_fasta_in_annotated_merged_fa}", "~{taxid_fasta_in_gsnap_hitsummary_tab}", "~{taxid_fasta_in_rapsearch2_hitsummary_tab}"]]' \
    --output-files '["taxid_annot.fasta"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "s3://idseq-database/taxonomy/2020-02-10/taxid-lineages.db"}' \
    --additional-attributes '{}'
  >>>
  output {
    File taxid_annot_fasta = "taxid_annot.fasta"
    File? output_read_count = "taxid_fasta_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateTaxidLocator {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File taxid_annot_fasta
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.generate_taxid_locator \
    --step-class PipelineStepGenerateTaxidLocator \
    --step-name taxid_locator_out \
    --input-files '[["~{taxid_annot_fasta}"]]' \
    --output-files '["taxid_annot_sorted_nt.fasta", "taxid_locations_nt.json", "taxid_annot_sorted_nr.fasta", "taxid_locations_nr.json", "taxid_annot_sorted_genus_nt.fasta", "taxid_locations_genus_nt.json", "taxid_annot_sorted_genus_nr.fasta", "taxid_locations_genus_nr.json", "taxid_annot_sorted_family_nt.fasta", "taxid_locations_family_nt.json", "taxid_annot_sorted_family_nr.fasta", "taxid_locations_family_nr.json", "taxid_locations_combined.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File taxid_annot_sorted_nt_fasta = "taxid_annot_sorted_nt.fasta"
    File taxid_locations_nt_json = "taxid_locations_nt.json"
    File taxid_annot_sorted_nr_fasta = "taxid_annot_sorted_nr.fasta"
    File taxid_locations_nr_json = "taxid_locations_nr.json"
    File taxid_annot_sorted_genus_nt_fasta = "taxid_annot_sorted_genus_nt.fasta"
    File taxid_locations_genus_nt_json = "taxid_locations_genus_nt.json"
    File taxid_annot_sorted_genus_nr_fasta = "taxid_annot_sorted_genus_nr.fasta"
    File taxid_locations_genus_nr_json = "taxid_locations_genus_nr.json"
    File taxid_annot_sorted_family_nt_fasta = "taxid_annot_sorted_family_nt.fasta"
    File taxid_locations_family_nt_json = "taxid_locations_family_nt.json"
    File taxid_annot_sorted_family_nr_fasta = "taxid_annot_sorted_family_nr.fasta"
    File taxid_locations_family_nr_json = "taxid_locations_family_nr.json"
    File taxid_locations_combined_json = "taxid_locations_combined.json"
    File? output_read_count = "taxid_locator_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateAlignmentViz {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File gsnap_m8_gsnap_deduped_m8
    File taxid_annot_sorted_nt_fasta
    File taxid_locations_nt_json
    File taxid_annot_sorted_nr_fasta
    File taxid_locations_nr_json
    File taxid_annot_sorted_genus_nt_fasta
    File taxid_locations_genus_nt_json
    File taxid_annot_sorted_genus_nr_fasta
    File taxid_locations_genus_nr_json
    File taxid_annot_sorted_family_nt_fasta
    File taxid_locations_family_nt_json
    File taxid_annot_sorted_family_nr_fasta
    File taxid_locations_family_nr_json
    File taxid_locations_combined_json
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.generate_alignment_viz \
    --step-class PipelineStepGenerateAlignmentViz \
    --step-name alignment_viz_out \
    --input-files '[["~{gsnap_m8_gsnap_deduped_m8}"], ["~{taxid_annot_sorted_nt_fasta}", "~{taxid_locations_nt_json}", "~{taxid_annot_sorted_nr_fasta}", "~{taxid_locations_nr_json}", "~{taxid_annot_sorted_genus_nt_fasta}", "~{taxid_locations_genus_nt_json}", "~{taxid_annot_sorted_genus_nr_fasta}", "~{taxid_locations_genus_nr_json}", "~{taxid_annot_sorted_family_nt_fasta}", "~{taxid_locations_family_nt_json}", "~{taxid_annot_sorted_family_nr_fasta}", "~{taxid_locations_family_nr_json}", "~{taxid_locations_combined_json}"]]' \
    --output-files '["align_viz.summary"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"nt_loc_db": "s3://idseq-database/alignment_data/2020-02-10/nt_loc.db", "nt_db": "s3://idseq-database/alignment_data/2020-02-10/nt"}' \
    --additional-attributes '{"nt_db": "s3://idseq-database/alignment_data/2020-02-10/nt"}'
  >>>
  output {
    File align_viz_summary = "align_viz.summary"
    File? output_read_count = "alignment_viz_out.count"
    Array[File] align_viz = glob("align_viz/*.align_viz.json")
  }
  runtime {
    docker: docker_image_id
  }
}

task RunSRST2 {
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
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.run_srst2 \
    --step-class PipelineStepRunSRST2 \
    --step-name srst2_out \
    --input-files '[["~{fastqs_0}", "~{fastqs_1}"]]' \
    --output-files '["out.log", "out__genes__ARGannot_r2__results.txt", "out__fullgenes__ARGannot_r2__results.txt", "amr_processed_results.csv", "amr_summary_results.csv", "output__.ARGannot_r2.sorted.bam"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"resist_gene_db": "s3://idseq-database/amr/ARGannot_r2.fasta", "resist_genome_bed": "s3://idseq-database/amr/argannot_genome.bed"}' \
    --additional-attributes '{"min_cov": 0, "n_threads": 16, "file_ext": "fastq"}'
  >>>
  output {
    File out_log = "out.log"
    File out__genes__ARGannot_r2__results_txt = "out__genes__ARGannot_r2__results.txt"
    File out__fullgenes__ARGannot_r2__results_txt = "out__fullgenes__ARGannot_r2__results.txt"
    File amr_processed_results_csv = "amr_processed_results.csv"
    File amr_summary_results_csv = "amr_summary_results.csv"
    File output___ARGannot_r2_sorted_bam = "output__.ARGannot_r2.sorted.bam"
    File? output_read_count = "srst2_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateCoverageViz {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File refined_gsnap_in_gsnap_reassigned_m8
    File refined_gsnap_in_gsnap_hitsummary2_tab
    File refined_gsnap_in_gsnap_blast_top_m8
    File contig_in_contig_coverage_json
    File contig_in_contig_stats_json
    File contig_in_contigs_fasta
    File gsnap_m8_gsnap_deduped_m8
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.generate_coverage_viz \
    --step-class PipelineStepGenerateCoverageViz \
    --step-name coverage_viz_out \
    --input-files '[["~{refined_gsnap_in_gsnap_reassigned_m8}", "~{refined_gsnap_in_gsnap_hitsummary2_tab}", "~{refined_gsnap_in_gsnap_blast_top_m8}"], ["~{contig_in_contig_coverage_json}", "~{contig_in_contig_stats_json}", "~{contig_in_contigs_fasta}"], ["~{gsnap_m8_gsnap_deduped_m8}"]]' \
    --output-files '["coverage_viz_summary.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"info_db": "s3://idseq-database/alignment_data/2020-02-10/nt_info.db"}' \
    --additional-attributes '{}'
  >>>
  output {
    File coverage_viz_summary_json = "coverage_viz_summary.json"
    File? output_read_count = "coverage_viz_out.count"
    Array[File] coverage_viz = glob("coverage_viz/*_coverage_viz.json")
  }
  runtime {
    docker: docker_image_id
  }
}

task NonhostFastq {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File fastqs_0
    File? fastqs_1
    File nonhost_fasta_refined_taxid_annot_fasta
    File cdhitdup_clusters_dedup1_fa_clstr
    File deduped_fasta_dedup1_fa
  }
  command<<<
  export AWS_DEFAULT_REGION=~{aws_region} DEPLOYMENT_ENVIRONMENT=~{deployment_env}
  if ! [[ -z "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  set -x
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.nonhost_fastq \
    --step-class PipelineStepNonhostFastq \
    --step-name nonhost_fastq_out \
    --input-files '[["~{fastqs_0}", "~{fastqs_1}"], ["~{nonhost_fasta_refined_taxid_annot_fasta}"], ["~{cdhitdup_clusters_dedup1_fa_clstr}"], ["~{deduped_fasta_dedup1_fa}"]]' \
    --output-files '["nonhost_R1.fastq", "nonhost_R2.fastq"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"use_taxon_whitelist": false}'
  >>>
  output {
    File nonhost_R1_fastq = "nonhost_R1.fastq"
    File nonhost_R2_fastq = "nonhost_R2.fastq"
    File? output_read_count = "nonhost_fastq_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

workflow idseq_experimental {
  input {
    String docker_image_id
    String aws_region
    String deployment_env
    String dag_branch
    String s3_wd_uri
    File taxid_fasta_in_annotated_merged_fa
    File taxid_fasta_in_gsnap_hitsummary_tab
    File taxid_fasta_in_rapsearch2_hitsummary_tab
    File gsnap_m8_gsnap_deduped_m8
    File refined_gsnap_in_gsnap_reassigned_m8
    File refined_gsnap_in_gsnap_hitsummary2_tab
    File refined_gsnap_in_gsnap_blast_top_m8
    File contig_in_contig_coverage_json
    File contig_in_contig_stats_json
    File contig_in_contigs_fasta
    File fastqs_0
    File? fastqs_1
    File nonhost_fasta_refined_taxid_annot_fasta
    File cdhitdup_clusters_dedup1_fa_clstr
    File deduped_fasta_dedup1_fa
  }

  call GenerateTaxidFasta {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      taxid_fasta_in_annotated_merged_fa = taxid_fasta_in_annotated_merged_fa,
      taxid_fasta_in_gsnap_hitsummary_tab = taxid_fasta_in_gsnap_hitsummary_tab,
      taxid_fasta_in_rapsearch2_hitsummary_tab = taxid_fasta_in_rapsearch2_hitsummary_tab
  }

  call GenerateTaxidLocator {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      taxid_annot_fasta = GenerateTaxidFasta.taxid_annot_fasta
  }

  call GenerateAlignmentViz {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      gsnap_m8_gsnap_deduped_m8 = gsnap_m8_gsnap_deduped_m8,
      taxid_annot_sorted_nt_fasta = GenerateTaxidLocator.taxid_annot_sorted_nt_fasta,
      taxid_locations_nt_json = GenerateTaxidLocator.taxid_locations_nt_json,
      taxid_annot_sorted_nr_fasta = GenerateTaxidLocator.taxid_annot_sorted_nr_fasta,
      taxid_locations_nr_json = GenerateTaxidLocator.taxid_locations_nr_json,
      taxid_annot_sorted_genus_nt_fasta = GenerateTaxidLocator.taxid_annot_sorted_genus_nt_fasta,
      taxid_locations_genus_nt_json = GenerateTaxidLocator.taxid_locations_genus_nt_json,
      taxid_annot_sorted_genus_nr_fasta = GenerateTaxidLocator.taxid_annot_sorted_genus_nr_fasta,
      taxid_locations_genus_nr_json = GenerateTaxidLocator.taxid_locations_genus_nr_json,
      taxid_annot_sorted_family_nt_fasta = GenerateTaxidLocator.taxid_annot_sorted_family_nt_fasta,
      taxid_locations_family_nt_json = GenerateTaxidLocator.taxid_locations_family_nt_json,
      taxid_annot_sorted_family_nr_fasta = GenerateTaxidLocator.taxid_annot_sorted_family_nr_fasta,
      taxid_locations_family_nr_json = GenerateTaxidLocator.taxid_locations_family_nr_json,
      taxid_locations_combined_json = GenerateTaxidLocator.taxid_locations_combined_json
  }

  call RunSRST2 {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      fastqs_0 = fastqs_0,
      fastqs_1 = fastqs_1
  }

  call GenerateCoverageViz {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      refined_gsnap_in_gsnap_reassigned_m8 = refined_gsnap_in_gsnap_reassigned_m8,
      refined_gsnap_in_gsnap_hitsummary2_tab = refined_gsnap_in_gsnap_hitsummary2_tab,
      refined_gsnap_in_gsnap_blast_top_m8 = refined_gsnap_in_gsnap_blast_top_m8,
      contig_in_contig_coverage_json = contig_in_contig_coverage_json,
      contig_in_contig_stats_json = contig_in_contig_stats_json,
      contig_in_contigs_fasta = contig_in_contigs_fasta,
      gsnap_m8_gsnap_deduped_m8 = gsnap_m8_gsnap_deduped_m8
  }

  call NonhostFastq {
    input:
      docker_image_id = docker_image_id,
      aws_region = aws_region,
      deployment_env = deployment_env,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      fastqs_0 = fastqs_0,
      fastqs_1 = fastqs_1,
      nonhost_fasta_refined_taxid_annot_fasta = nonhost_fasta_refined_taxid_annot_fasta,
      cdhitdup_clusters_dedup1_fa_clstr = cdhitdup_clusters_dedup1_fa_clstr,
      deduped_fasta_dedup1_fa = deduped_fasta_dedup1_fa
  }

  output {
    File taxid_fasta_out_taxid_annot_fasta = GenerateTaxidFasta.taxid_annot_fasta
    File? taxid_fasta_out_count = GenerateTaxidFasta.output_read_count
    File taxid_locator_out_taxid_annot_sorted_nt_fasta = GenerateTaxidLocator.taxid_annot_sorted_nt_fasta
    File taxid_locator_out_taxid_locations_nt_json = GenerateTaxidLocator.taxid_locations_nt_json
    File taxid_locator_out_taxid_annot_sorted_nr_fasta = GenerateTaxidLocator.taxid_annot_sorted_nr_fasta
    File taxid_locator_out_taxid_locations_nr_json = GenerateTaxidLocator.taxid_locations_nr_json
    File taxid_locator_out_taxid_annot_sorted_genus_nt_fasta = GenerateTaxidLocator.taxid_annot_sorted_genus_nt_fasta
    File taxid_locator_out_taxid_locations_genus_nt_json = GenerateTaxidLocator.taxid_locations_genus_nt_json
    File taxid_locator_out_taxid_annot_sorted_genus_nr_fasta = GenerateTaxidLocator.taxid_annot_sorted_genus_nr_fasta
    File taxid_locator_out_taxid_locations_genus_nr_json = GenerateTaxidLocator.taxid_locations_genus_nr_json
    File taxid_locator_out_taxid_annot_sorted_family_nt_fasta = GenerateTaxidLocator.taxid_annot_sorted_family_nt_fasta
    File taxid_locator_out_taxid_locations_family_nt_json = GenerateTaxidLocator.taxid_locations_family_nt_json
    File taxid_locator_out_taxid_annot_sorted_family_nr_fasta = GenerateTaxidLocator.taxid_annot_sorted_family_nr_fasta
    File taxid_locator_out_taxid_locations_family_nr_json = GenerateTaxidLocator.taxid_locations_family_nr_json
    File taxid_locator_out_taxid_locations_combined_json = GenerateTaxidLocator.taxid_locations_combined_json
    File? taxid_locator_out_count = GenerateTaxidLocator.output_read_count
    File alignment_viz_out_align_viz_summary = GenerateAlignmentViz.align_viz_summary
    File? alignment_viz_out_count = GenerateAlignmentViz.output_read_count
    File srst2_out_out_log = RunSRST2.out_log
    File srst2_out_out__genes__ARGannot_r2__results_txt = RunSRST2.out__genes__ARGannot_r2__results_txt
    File srst2_out_out__fullgenes__ARGannot_r2__results_txt = RunSRST2.out__fullgenes__ARGannot_r2__results_txt
    File srst2_out_amr_processed_results_csv = RunSRST2.amr_processed_results_csv
    File srst2_out_amr_summary_results_csv = RunSRST2.amr_summary_results_csv
    File srst2_out_output___ARGannot_r2_sorted_bam = RunSRST2.output___ARGannot_r2_sorted_bam
    File? srst2_out_count = RunSRST2.output_read_count
    File coverage_viz_out_coverage_viz_summary_json = GenerateCoverageViz.coverage_viz_summary_json
    File? coverage_viz_out_count = GenerateCoverageViz.output_read_count
    File nonhost_fastq_out_nonhost_R1_fastq = NonhostFastq.nonhost_R1_fastq
    File nonhost_fastq_out_nonhost_R2_fastq = NonhostFastq.nonhost_R2_fastq
    File? nonhost_fastq_out_count = NonhostFastq.output_read_count
    Array[File] align_viz = GenerateAlignmentViz.align_viz
    Array[File] coverage_viz = GenerateCoverageViz.coverage_viz
  }
}
