version 1.0

task GenerateTaxidFasta {
  input {
    String docker_image_id
    String s3_wd_uri
    File taxid_fasta_in_annotated_merged_fa
    File taxid_fasta_in_gsnap_hitsummary_tab
    File taxid_fasta_in_rapsearch2_hitsummary_tab
    File lineage_db
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.generate_taxid_fasta \
    --step-class PipelineStepGenerateTaxidFasta \
    --step-name taxid_fasta_out \
    --input-files '[["~{taxid_fasta_in_annotated_merged_fa}", "~{taxid_fasta_in_gsnap_hitsummary_tab}", "~{taxid_fasta_in_rapsearch2_hitsummary_tab}"]]' \
    --output-files '["taxid_annot.fasta"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}"}' \
    --additional-attributes '{}'
  >>>
  output {
    String step_description_md = read_string("taxid_fasta_out.description.md")
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
    String s3_wd_uri
    File taxid_annot_fasta
  }
  command<<<
  set -euxo pipefail
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
    String step_description_md = read_string("taxid_locator_out.description.md")
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

task GenerateCoverageViz {
  input {
    String docker_image_id
    String s3_wd_uri
    File refined_gsnap_in_gsnap_reassigned_m8
    File refined_gsnap_in_gsnap_hitsummary2_tab
    File refined_gsnap_in_gsnap_blast_top_m8
    File contig_in_contig_coverage_json
    File contig_in_contig_stats_json
    File contig_in_contigs_fasta
    File gsnap_m8_gsnap_deduped_m8
    File nt_info_db
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.generate_coverage_viz \
    --step-class PipelineStepGenerateCoverageViz \
    --step-name coverage_viz_out \
    --input-files '[["~{refined_gsnap_in_gsnap_reassigned_m8}", "~{refined_gsnap_in_gsnap_hitsummary2_tab}", "~{refined_gsnap_in_gsnap_blast_top_m8}"], ["~{contig_in_contig_coverage_json}", "~{contig_in_contig_stats_json}", "~{contig_in_contigs_fasta}"], ["~{gsnap_m8_gsnap_deduped_m8}"]]' \
    --output-files '["coverage_viz_summary.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"info_db": "~{nt_info_db}"}' \
    --additional-attributes '{}'
  >>>
  output {
    String step_description_md = read_string("coverage_viz_out.description.md")
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
    String s3_wd_uri
    Array[File] fastqs
    File nonhost_fasta_refined_taxid_annot_fasta
    File duplicate_clusters_csv
    Boolean use_taxon_whitelist
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name experimental \
    --step-module idseq_dag.steps.nonhost_fastq \
    --step-class PipelineStepNonhostFastq \
    --step-name nonhost_fastq_out \
    --input-files '[["~{sep='","' fastqs}"], ["~{nonhost_fasta_refined_taxid_annot_fasta}"], ["~{duplicate_clusters_csv}"]]' \
    --output-files '[~{if length(fastqs) == 2 then '"nonhost_R1.fastq", "nonhost_R2.fastq"' else '"nonhost_R1.fastq"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"use_taxon_whitelist": ~{use_taxon_whitelist}}'
  >>>
  output {
    String step_description_md = read_string("nonhost_fastq_out.description.md")
    File nonhost_R1_fastq = "nonhost_R1.fastq"
    File? nonhost_R2_fastq = "nonhost_R2.fastq"
    File? output_read_count = "nonhost_fastq_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

workflow czid_experimental {
  input {
    String docker_image_id
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
    File duplicate_clusters_csv
    String file_ext = "fastq"
    String nt_db = "s3://czid-public-references/ncbi-sources/2021-01-22/nt"
    File nt_loc_db = "s3://czid-public-references/alignment_data/2021-01-22/nt_loc.db"
    File nt_info_db = "s3://czid-public-references/alignment_data/2021-01-22/nt_info.db"
    File lineage_db = "s3://czid-public-references/taxonomy/2021-01-22/taxid-lineages.db"
    File resist_genome_db = "s3://czid-public-references/amr/ARGannot_r2.fasta"
    File resist_genome_bed = "s3://czid-public-references/amr/argannot_genome.bed"
    Boolean use_taxon_whitelist = false
  }

  call GenerateTaxidFasta {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      taxid_fasta_in_annotated_merged_fa = taxid_fasta_in_annotated_merged_fa,
      taxid_fasta_in_gsnap_hitsummary_tab = taxid_fasta_in_gsnap_hitsummary_tab,
      taxid_fasta_in_rapsearch2_hitsummary_tab = taxid_fasta_in_rapsearch2_hitsummary_tab,
      lineage_db = lineage_db
  }

  call GenerateTaxidLocator {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      taxid_annot_fasta = GenerateTaxidFasta.taxid_annot_fasta
  }

  call GenerateCoverageViz {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      refined_gsnap_in_gsnap_reassigned_m8 = refined_gsnap_in_gsnap_reassigned_m8,
      refined_gsnap_in_gsnap_hitsummary2_tab = refined_gsnap_in_gsnap_hitsummary2_tab,
      refined_gsnap_in_gsnap_blast_top_m8 = refined_gsnap_in_gsnap_blast_top_m8,
      contig_in_contig_coverage_json = contig_in_contig_coverage_json,
      contig_in_contig_stats_json = contig_in_contig_stats_json,
      contig_in_contigs_fasta = contig_in_contigs_fasta,
      gsnap_m8_gsnap_deduped_m8 = gsnap_m8_gsnap_deduped_m8,
      nt_info_db = nt_info_db
  }

  call NonhostFastq {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      fastqs = select_all([fastqs_0, fastqs_1]),
      nonhost_fasta_refined_taxid_annot_fasta = nonhost_fasta_refined_taxid_annot_fasta,
      duplicate_clusters_csv = duplicate_clusters_csv,
      use_taxon_whitelist = use_taxon_whitelist
  }

  output {
    File taxid_fasta_out_taxid_annot_fasta = GenerateTaxidFasta.taxid_annot_fasta
    File? taxid_fasta_out_count = GenerateTaxidFasta.output_read_count
    File? taxid_locator_out_count = GenerateTaxidLocator.output_read_count
    File coverage_viz_out_coverage_viz_summary_json = GenerateCoverageViz.coverage_viz_summary_json
    File? coverage_viz_out_count = GenerateCoverageViz.output_read_count
    File nonhost_fastq_out_nonhost_R1_fastq = NonhostFastq.nonhost_R1_fastq
    File? nonhost_fastq_out_nonhost_R2_fastq = NonhostFastq.nonhost_R2_fastq
    File? nonhost_fastq_out_count = NonhostFastq.output_read_count
    Array[File] coverage_viz = GenerateCoverageViz.coverage_viz
  }
}
