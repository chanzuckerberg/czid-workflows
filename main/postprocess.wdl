version 1.0

task RunAssembly {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    Array[File] host_filter_out_gsnap_filter_fa
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.run_assembly \
    --step-class PipelineStepRunAssembly \
    --step-name assembly_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["assembly/contigs.fasta", "assembly/scaffolds.fasta", "assembly/read-contig.sam", "assembly/contig_stats.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"memory": 200}'
  >>>
  output {
    File assembly_contigs_fasta = "assembly/contigs.fasta"
    File assembly_scaffolds_fasta = "assembly/scaffolds.fasta"
    File assembly_read_contig_sam = "assembly/read-contig.sam"
    File assembly_contig_stats_json = "assembly/contig_stats.json"
    File? output_read_count = "assembly_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateCoverageStats {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File assembly_contigs_fasta
    File assembly_scaffolds_fasta
    File assembly_read_contig_sam
    File assembly_contig_stats_json
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.generate_coverage_stats \
    --step-class PipelineStepGenerateCoverageStats \
    --step-name coverage_out \
    --input-files '[["~{assembly_contigs_fasta}", "~{assembly_scaffolds_fasta}", "~{assembly_read_contig_sam}", "~{assembly_contig_stats_json}"]]' \
    --output-files '["assembly/contig_coverage.json", "assembly/contig_coverage_summary.csv"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File assembly_contig_coverage_json = "assembly/contig_coverage.json"
    File assembly_contig_coverage_summary_csv = "assembly/contig_coverage_summary.csv"
    File? output_read_count = "coverage_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task DownloadAccessions_gsnap_accessions_out {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File gsnap_out_gsnap_m8
    File gsnap_out_gsnap_deduped_m8
    File gsnap_out_gsnap_hitsummary_tab
    File gsnap_out_gsnap_counts_with_dcr_json
    String nt_db
    String nt_loc_db
    String lineage_db
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.download_accessions \
    --step-class PipelineStepDownloadAccessions \
    --step-name gsnap_accessions_out \
    --input-files '[["~{gsnap_out_gsnap_m8}", "~{gsnap_out_gsnap_deduped_m8}", "~{gsnap_out_gsnap_hitsummary_tab}", "~{gsnap_out_gsnap_counts_with_dcr_json}"]]' \
    --output-files '["assembly/nt.refseq.fasta"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "loc_db": "~{nt_loc_db}", "db": "~{nt_db}"}' \
    --additional-attributes '{"db": "~{nt_db}", "db_type": "nt"}'
  >>>
  output {
    File assembly_nt_refseq_fasta = "assembly/nt.refseq.fasta"
    File? output_read_count = "gsnap_accessions_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task DownloadAccessions_rapsearch2_accessions_out {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File rapsearch2_out_rapsearch2_m8
    File rapsearch2_out_rapsearch2_deduped_m8
    File rapsearch2_out_rapsearch2_hitsummary_tab
    File rapsearch2_out_rapsearch2_counts_with_dcr_json
    String lineage_db
    String nr_loc_db
    String nr_db
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.download_accessions \
    --step-class PipelineStepDownloadAccessions \
    --step-name rapsearch2_accessions_out \
    --input-files '[["~{rapsearch2_out_rapsearch2_m8}", "~{rapsearch2_out_rapsearch2_deduped_m8}", "~{rapsearch2_out_rapsearch2_hitsummary_tab}", "~{rapsearch2_out_rapsearch2_counts_with_dcr_json}"]]' \
    --output-files '["assembly/nr.refseq.fasta"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "loc_db": "~{nr_loc_db}", "db": "~{nr_db}"}' \
    --additional-attributes '{"db": "~{nr_db}", "db_type": "nr"}'
  >>>
  output {
    File assembly_nr_refseq_fasta = "assembly/nr.refseq.fasta"
    File? output_read_count = "rapsearch2_accessions_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task BlastContigs_refined_gsnap_out {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File gsnap_out_gsnap_m8
    File gsnap_out_gsnap_deduped_m8
    File gsnap_out_gsnap_hitsummary_tab
    File gsnap_out_gsnap_counts_with_dcr_json
    File assembly_contigs_fasta
    File assembly_scaffolds_fasta
    File assembly_read_contig_sam
    File assembly_contig_stats_json
    File assembly_nt_refseq_fasta
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
    String lineage_db
    String taxon_blacklist
    String? deuterostome_db
    Boolean use_deuterostome_filter
    Boolean use_taxon_whitelist
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.blast_contigs \
    --step-class PipelineStepBlastContigs \
    --step-name refined_gsnap_out \
    --input-files '[["~{gsnap_out_gsnap_m8}", "~{gsnap_out_gsnap_deduped_m8}", "~{gsnap_out_gsnap_hitsummary_tab}", "~{gsnap_out_gsnap_counts_with_dcr_json}"], ["~{assembly_contigs_fasta}", "~{assembly_scaffolds_fasta}", "~{assembly_read_contig_sam}", "~{assembly_contig_stats_json}"], ["~{assembly_nt_refseq_fasta}"], ["~{cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["assembly/gsnap.blast.m8", "assembly/gsnap.reassigned.m8", "assembly/gsnap.hitsummary2.tab", "assembly/refined_gsnap_counts_with_dcr.json", "assembly/gsnap_contig_summary.json", "assembly/gsnap.blast.top.m8"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "taxon_blacklist": "~{taxon_blacklist}", "deuterostome_db": "~{if use_deuterostome_filter then '~{deuterostome_db}' else ''}"}' \
    --additional-attributes '{"db_type": "nt", "use_taxon_whitelist": ~{use_taxon_whitelist}}'
  >>>
  output {
    File assembly_gsnap_blast_m8 = "assembly/gsnap.blast.m8"
    File assembly_gsnap_reassigned_m8 = "assembly/gsnap.reassigned.m8"
    File assembly_gsnap_hitsummary2_tab = "assembly/gsnap.hitsummary2.tab"
    File assembly_refined_gsnap_counts_with_dcr_json = "assembly/refined_gsnap_counts_with_dcr.json"
    File assembly_gsnap_contig_summary_json = "assembly/gsnap_contig_summary.json"
    File assembly_gsnap_blast_top_m8 = "assembly/gsnap.blast.top.m8"
    File? output_read_count = "refined_gsnap_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task BlastContigs_refined_rapsearch2_out {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File rapsearch2_out_rapsearch2_m8
    File rapsearch2_out_rapsearch2_deduped_m8
    File rapsearch2_out_rapsearch2_hitsummary_tab
    File rapsearch2_out_rapsearch2_counts_with_dcr_json
    File assembly_contigs_fasta
    File assembly_scaffolds_fasta
    File assembly_read_contig_sam
    File assembly_contig_stats_json
    File assembly_nr_refseq_fasta
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
    String lineage_db
    String taxon_blacklist
    Boolean use_taxon_whitelist
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.blast_contigs \
    --step-class PipelineStepBlastContigs \
    --step-name refined_rapsearch2_out \
    --input-files '[["~{rapsearch2_out_rapsearch2_m8}", "~{rapsearch2_out_rapsearch2_deduped_m8}", "~{rapsearch2_out_rapsearch2_hitsummary_tab}", "~{rapsearch2_out_rapsearch2_counts_with_dcr_json}"], ["~{assembly_contigs_fasta}", "~{assembly_scaffolds_fasta}", "~{assembly_read_contig_sam}", "~{assembly_contig_stats_json}"], ["~{assembly_nr_refseq_fasta}"], ["~{cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["assembly/rapsearch2.blast.m8", "assembly/rapsearch2.reassigned.m8", "assembly/rapsearch2.hitsummary2.tab", "assembly/refined_rapsearch2_counts_with_dcr.json", "assembly/rapsearch2_contig_summary.json", "assembly/rapsearch2.blast.top.m8"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}", "taxon_blacklist": "~{taxon_blacklist}"}' \
    --additional-attributes '{"db_type": "nr", "use_taxon_whitelist": ~{use_taxon_whitelist}}'
  >>>
  output {
    File assembly_rapsearch2_blast_m8 = "assembly/rapsearch2.blast.m8"
    File assembly_rapsearch2_reassigned_m8 = "assembly/rapsearch2.reassigned.m8"
    File assembly_rapsearch2_hitsummary2_tab = "assembly/rapsearch2.hitsummary2.tab"
    File assembly_refined_rapsearch2_counts_with_dcr_json = "assembly/refined_rapsearch2_counts_with_dcr.json"
    File assembly_rapsearch2_contig_summary_json = "assembly/rapsearch2_contig_summary.json"
    File assembly_rapsearch2_blast_top_m8 = "assembly/rapsearch2.blast.top.m8"
    File? output_read_count = "refined_rapsearch2_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task CombineTaxonCounts {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File assembly_gsnap_blast_m8
    File assembly_gsnap_reassigned_m8
    File assembly_gsnap_hitsummary2_tab
    File assembly_refined_gsnap_counts_with_dcr_json
    File assembly_gsnap_contig_summary_json
    File assembly_gsnap_blast_top_m8
    File assembly_rapsearch2_blast_m8
    File assembly_rapsearch2_reassigned_m8
    File assembly_rapsearch2_hitsummary2_tab
    File assembly_refined_rapsearch2_counts_with_dcr_json
    File assembly_rapsearch2_contig_summary_json
    File assembly_rapsearch2_blast_top_m8
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.combine_taxon_counts \
    --step-class PipelineStepCombineTaxonCounts \
    --step-name refined_taxon_count_out \
    --input-files '[["~{assembly_gsnap_blast_m8}", "~{assembly_gsnap_reassigned_m8}", "~{assembly_gsnap_hitsummary2_tab}", "~{assembly_refined_gsnap_counts_with_dcr_json}", "~{assembly_gsnap_contig_summary_json}", "~{assembly_gsnap_blast_top_m8}"], ["~{assembly_rapsearch2_blast_m8}", "~{assembly_rapsearch2_reassigned_m8}", "~{assembly_rapsearch2_hitsummary2_tab}", "~{assembly_refined_rapsearch2_counts_with_dcr_json}", "~{assembly_rapsearch2_contig_summary_json}", "~{assembly_rapsearch2_blast_top_m8}"]]' \
    --output-files '["assembly/refined_taxon_counts_with_dcr.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File assembly_refined_taxon_counts_with_dcr_json = "assembly/refined_taxon_counts_with_dcr.json"
    File? output_read_count = "refined_taxon_count_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task CombineJson {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File assembly_gsnap_blast_m8
    File assembly_gsnap_reassigned_m8
    File assembly_gsnap_hitsummary2_tab
    File assembly_refined_gsnap_counts_with_dcr_json
    File assembly_gsnap_contig_summary_json
    File assembly_gsnap_blast_top_m8
    File assembly_rapsearch2_blast_m8
    File assembly_rapsearch2_reassigned_m8
    File assembly_rapsearch2_hitsummary2_tab
    File assembly_refined_rapsearch2_counts_with_dcr_json
    File assembly_rapsearch2_contig_summary_json
    File assembly_rapsearch2_blast_top_m8
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.combine_json \
    --step-class PipelineStepCombineJson \
    --step-name contig_summary_out \
    --input-files '[["~{assembly_gsnap_blast_m8}", "~{assembly_gsnap_reassigned_m8}", "~{assembly_gsnap_hitsummary2_tab}", "~{assembly_refined_gsnap_counts_with_dcr_json}", "~{assembly_gsnap_contig_summary_json}", "~{assembly_gsnap_blast_top_m8}"], ["~{assembly_rapsearch2_blast_m8}", "~{assembly_rapsearch2_reassigned_m8}", "~{assembly_rapsearch2_hitsummary2_tab}", "~{assembly_refined_rapsearch2_counts_with_dcr_json}", "~{assembly_rapsearch2_contig_summary_json}", "~{assembly_rapsearch2_blast_top_m8}"]]' \
    --output-files '["assembly/combined_contig_summary.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"field_idx": 4}'
  >>>
  output {
    File assembly_combined_contig_summary_json = "assembly/combined_contig_summary.json"
    File? output_read_count = "contig_summary_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateAnnotatedFasta {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    Array[File] host_filter_out_gsnap_filter_fa
    File assembly_gsnap_blast_m8
    File assembly_gsnap_reassigned_m8
    File assembly_gsnap_hitsummary2_tab
    File assembly_refined_gsnap_counts_with_dcr_json
    File assembly_gsnap_contig_summary_json
    File assembly_gsnap_blast_top_m8
    File assembly_rapsearch2_blast_m8
    File assembly_rapsearch2_reassigned_m8
    File assembly_rapsearch2_hitsummary2_tab
    File assembly_refined_rapsearch2_counts_with_dcr_json
    File assembly_rapsearch2_contig_summary_json
    File assembly_rapsearch2_blast_top_m8
    File cdhitdup_out_dedup1_fa_clstr
    File cdhitdup_out_dedup1_fa
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.generate_annotated_fasta \
    --step-class PipelineStepGenerateAnnotatedFasta \
    --step-name refined_annotated_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{assembly_gsnap_blast_m8}", "~{assembly_gsnap_reassigned_m8}", "~{assembly_gsnap_hitsummary2_tab}", "~{assembly_refined_gsnap_counts_with_dcr_json}", "~{assembly_gsnap_contig_summary_json}", "~{assembly_gsnap_blast_top_m8}"], ["~{assembly_rapsearch2_blast_m8}", "~{assembly_rapsearch2_reassigned_m8}", "~{assembly_rapsearch2_hitsummary2_tab}", "~{assembly_refined_rapsearch2_counts_with_dcr_json}", "~{assembly_rapsearch2_contig_summary_json}", "~{assembly_rapsearch2_blast_top_m8}"], ["~{cdhitdup_out_dedup1_fa_clstr}", "~{cdhitdup_out_dedup1_fa}"], ["~{cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv}"]]' \
    --output-files '["assembly/refined_annotated_merged.fa", "assembly/refined_unidentified.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File assembly_refined_annotated_merged_fa = "assembly/refined_annotated_merged.fa"
    File assembly_refined_unidentified_fa = "assembly/refined_unidentified.fa"
    File? output_read_count = "refined_annotated_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateTaxidFasta {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File assembly_refined_annotated_merged_fa
    File assembly_refined_unidentified_fa
    File assembly_gsnap_blast_m8
    File assembly_gsnap_reassigned_m8
    File assembly_gsnap_hitsummary2_tab
    File assembly_refined_gsnap_counts_with_dcr_json
    File assembly_gsnap_contig_summary_json
    File assembly_gsnap_blast_top_m8
    File assembly_rapsearch2_blast_m8
    File assembly_rapsearch2_reassigned_m8
    File assembly_rapsearch2_hitsummary2_tab
    File assembly_refined_rapsearch2_counts_with_dcr_json
    File assembly_rapsearch2_contig_summary_json
    File assembly_rapsearch2_blast_top_m8
    String lineage_db
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.generate_taxid_fasta \
    --step-class PipelineStepGenerateTaxidFasta \
    --step-name refined_taxid_fasta_out \
    --input-files '[["~{assembly_refined_annotated_merged_fa}", "~{assembly_refined_unidentified_fa}"], ["~{assembly_gsnap_blast_m8}", "~{assembly_gsnap_reassigned_m8}", "~{assembly_gsnap_hitsummary2_tab}", "~{assembly_refined_gsnap_counts_with_dcr_json}", "~{assembly_gsnap_contig_summary_json}", "~{assembly_gsnap_blast_top_m8}"], ["~{assembly_rapsearch2_blast_m8}", "~{assembly_rapsearch2_reassigned_m8}", "~{assembly_rapsearch2_hitsummary2_tab}", "~{assembly_refined_rapsearch2_counts_with_dcr_json}", "~{assembly_rapsearch2_contig_summary_json}", "~{assembly_rapsearch2_blast_top_m8}"]]' \
    --output-files '["assembly/refined_taxid_annot.fasta"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"lineage_db": "~{lineage_db}"}' \
    --additional-attributes '{}'
  >>>
  output {
    File assembly_refined_taxid_annot_fasta = "assembly/refined_taxid_annot.fasta"
    File? output_read_count = "refined_taxid_fasta_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task GenerateTaxidLocator {
  input {
    String docker_image_id
    String dag_branch
    String s3_wd_uri
    File assembly_refined_taxid_annot_fasta
  }
  command<<<
  set -euxo pipefail
  if [[ -n "~{dag_branch}" ]]; then
    pip3 install --upgrade https://github.com/chanzuckerberg/idseq-dag/archive/~{dag_branch}.tar.gz
  fi
  idseq-dag-run-step --workflow-name postprocess \
    --step-module idseq_dag.steps.generate_taxid_locator \
    --step-class PipelineStepGenerateTaxidLocator \
    --step-name refined_taxid_locator_out \
    --input-files '[["~{assembly_refined_taxid_annot_fasta}"]]' \
    --output-files '["assembly/refined_taxid_annot_sorted_nt.fasta", "assembly/refined_taxid_locations_nt.json", "assembly/refined_taxid_annot_sorted_nr.fasta", "assembly/refined_taxid_locations_nr.json", "assembly/refined_taxid_annot_sorted_genus_nt.fasta", "assembly/refined_taxid_locations_genus_nt.json", "assembly/refined_taxid_annot_sorted_genus_nr.fasta", "assembly/refined_taxid_locations_genus_nr.json", "assembly/refined_taxid_annot_sorted_family_nt.fasta", "assembly/refined_taxid_locations_family_nt.json", "assembly/refined_taxid_annot_sorted_family_nr.fasta", "assembly/refined_taxid_locations_family_nr.json", "assembly/refined_taxid_locations_combined.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    File assembly_refined_taxid_annot_sorted_nt_fasta = "assembly/refined_taxid_annot_sorted_nt.fasta"
    File assembly_refined_taxid_locations_nt_json = "assembly/refined_taxid_locations_nt.json"
    File assembly_refined_taxid_annot_sorted_nr_fasta = "assembly/refined_taxid_annot_sorted_nr.fasta"
    File assembly_refined_taxid_locations_nr_json = "assembly/refined_taxid_locations_nr.json"
    File assembly_refined_taxid_annot_sorted_genus_nt_fasta = "assembly/refined_taxid_annot_sorted_genus_nt.fasta"
    File assembly_refined_taxid_locations_genus_nt_json = "assembly/refined_taxid_locations_genus_nt.json"
    File assembly_refined_taxid_annot_sorted_genus_nr_fasta = "assembly/refined_taxid_annot_sorted_genus_nr.fasta"
    File assembly_refined_taxid_locations_genus_nr_json = "assembly/refined_taxid_locations_genus_nr.json"
    File assembly_refined_taxid_annot_sorted_family_nt_fasta = "assembly/refined_taxid_annot_sorted_family_nt.fasta"
    File assembly_refined_taxid_locations_family_nt_json = "assembly/refined_taxid_locations_family_nt.json"
    File assembly_refined_taxid_annot_sorted_family_nr_fasta = "assembly/refined_taxid_annot_sorted_family_nr.fasta"
    File assembly_refined_taxid_locations_family_nr_json = "assembly/refined_taxid_locations_family_nr.json"
    File assembly_refined_taxid_locations_combined_json = "assembly/refined_taxid_locations_combined.json"
    File? output_read_count = "refined_taxid_locator_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

workflow idseq_postprocess {
  input {
    String docker_image_id
    String dag_branch = ""
    String s3_wd_uri
    File host_filter_out_gsnap_filter_1_fa
    File? host_filter_out_gsnap_filter_2_fa
    File? host_filter_out_gsnap_filter_merged_fa
    File gsnap_out_gsnap_m8
    File gsnap_out_gsnap_deduped_m8
    File gsnap_out_gsnap_hitsummary_tab
    File gsnap_out_gsnap_counts_with_dcr_json
    File rapsearch2_out_rapsearch2_m8
    File rapsearch2_out_rapsearch2_deduped_m8
    File rapsearch2_out_rapsearch2_hitsummary_tab
    File rapsearch2_out_rapsearch2_counts_with_dcr_json
    File cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
    File cdhitdup_out_dedup1_fa_clstr
    File cdhitdup_out_dedup1_fa
    String idseq_db_bucket = "idseq-database"
    String index_version = "2020-04-20"
    String nt_db = "s3://~{idseq_db_bucket}/ncbi-sources/~{index_version}/nt"
    String nt_loc_db = "s3://~{idseq_db_bucket}/ncbi-sources/~{index_version}/nt_loc.db"
    String nr_db = "s3://~{idseq_db_bucket}/ncbi-sources/~{index_version}/nr"
    String nr_loc_db = "s3://~{idseq_db_bucket}/ncbi-sources/~{index_version}/nr_loc.db"
    String lineage_db = "s3://~{idseq_db_bucket}/taxonomy/~{index_version}/taxid-lineages.db"
    String taxon_blacklist = "s3://~{idseq_db_bucket}/taxonomy/~{index_version}/taxon_blacklist.txt"
    String deuterostome_db = "s3://~{idseq_db_bucket}/taxonomy/~{index_version}/deuterostome_taxids.txt"
    Boolean use_deuterostome_filter = true
    Boolean use_taxon_whitelist = false
  }

  call RunAssembly {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
  }

  call GenerateCoverageStats {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      assembly_contigs_fasta = RunAssembly.assembly_contigs_fasta,
      assembly_scaffolds_fasta = RunAssembly.assembly_scaffolds_fasta,
      assembly_read_contig_sam = RunAssembly.assembly_read_contig_sam,
      assembly_contig_stats_json = RunAssembly.assembly_contig_stats_json
  }

  call DownloadAccessions_gsnap_accessions_out {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      gsnap_out_gsnap_m8 = gsnap_out_gsnap_m8,
      gsnap_out_gsnap_deduped_m8 = gsnap_out_gsnap_deduped_m8,
      gsnap_out_gsnap_hitsummary_tab = gsnap_out_gsnap_hitsummary_tab,
      gsnap_out_gsnap_counts_with_dcr_json = gsnap_out_gsnap_counts_with_dcr_json,
      nt_db = nt_db,
      nt_loc_db = nt_loc_db,
      lineage_db = lineage_db
  }

  call DownloadAccessions_rapsearch2_accessions_out {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      rapsearch2_out_rapsearch2_m8 = rapsearch2_out_rapsearch2_m8,
      rapsearch2_out_rapsearch2_deduped_m8 = rapsearch2_out_rapsearch2_deduped_m8,
      rapsearch2_out_rapsearch2_hitsummary_tab = rapsearch2_out_rapsearch2_hitsummary_tab,
      rapsearch2_out_rapsearch2_counts_with_dcr_json = rapsearch2_out_rapsearch2_counts_with_dcr_json,
      nr_db = nr_db,
      nr_loc_db = nr_loc_db,
      lineage_db = lineage_db
  }

  call BlastContigs_refined_gsnap_out {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      gsnap_out_gsnap_m8 = gsnap_out_gsnap_m8,
      gsnap_out_gsnap_deduped_m8 = gsnap_out_gsnap_deduped_m8,
      gsnap_out_gsnap_hitsummary_tab = gsnap_out_gsnap_hitsummary_tab,
      gsnap_out_gsnap_counts_with_dcr_json = gsnap_out_gsnap_counts_with_dcr_json,
      assembly_contigs_fasta = RunAssembly.assembly_contigs_fasta,
      assembly_scaffolds_fasta = RunAssembly.assembly_scaffolds_fasta,
      assembly_read_contig_sam = RunAssembly.assembly_read_contig_sam,
      assembly_contig_stats_json = RunAssembly.assembly_contig_stats_json,
      assembly_nt_refseq_fasta = DownloadAccessions_gsnap_accessions_out.assembly_nt_refseq_fasta,
      cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv,
      lineage_db = lineage_db,
      taxon_blacklist = taxon_blacklist,
      deuterostome_db = deuterostome_db,
      use_deuterostome_filter = use_deuterostome_filter,
      use_taxon_whitelist = use_taxon_whitelist
  }

  call BlastContigs_refined_rapsearch2_out {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      rapsearch2_out_rapsearch2_m8 = rapsearch2_out_rapsearch2_m8,
      rapsearch2_out_rapsearch2_deduped_m8 = rapsearch2_out_rapsearch2_deduped_m8,
      rapsearch2_out_rapsearch2_hitsummary_tab = rapsearch2_out_rapsearch2_hitsummary_tab,
      rapsearch2_out_rapsearch2_counts_with_dcr_json = rapsearch2_out_rapsearch2_counts_with_dcr_json,
      assembly_contigs_fasta = RunAssembly.assembly_contigs_fasta,
      assembly_scaffolds_fasta = RunAssembly.assembly_scaffolds_fasta,
      assembly_read_contig_sam = RunAssembly.assembly_read_contig_sam,
      assembly_contig_stats_json = RunAssembly.assembly_contig_stats_json,
      assembly_nr_refseq_fasta = DownloadAccessions_rapsearch2_accessions_out.assembly_nr_refseq_fasta,
      cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv,
      lineage_db = lineage_db,
      taxon_blacklist = taxon_blacklist,
      use_taxon_whitelist = use_taxon_whitelist
  }

  call CombineTaxonCounts {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      assembly_gsnap_blast_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_m8,
      assembly_gsnap_reassigned_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_reassigned_m8,
      assembly_gsnap_hitsummary2_tab = BlastContigs_refined_gsnap_out.assembly_gsnap_hitsummary2_tab,
      assembly_refined_gsnap_counts_with_dcr_json = BlastContigs_refined_gsnap_out.assembly_refined_gsnap_counts_with_dcr_json,
      assembly_gsnap_contig_summary_json = BlastContigs_refined_gsnap_out.assembly_gsnap_contig_summary_json,
      assembly_gsnap_blast_top_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_top_m8,
      assembly_rapsearch2_blast_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_m8,
      assembly_rapsearch2_reassigned_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_reassigned_m8,
      assembly_rapsearch2_hitsummary2_tab = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_hitsummary2_tab,
      assembly_refined_rapsearch2_counts_with_dcr_json = BlastContigs_refined_rapsearch2_out.assembly_refined_rapsearch2_counts_with_dcr_json,
      assembly_rapsearch2_contig_summary_json = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_contig_summary_json,
      assembly_rapsearch2_blast_top_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_top_m8
  }

  call CombineJson {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      assembly_gsnap_blast_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_m8,
      assembly_gsnap_reassigned_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_reassigned_m8,
      assembly_gsnap_hitsummary2_tab = BlastContigs_refined_gsnap_out.assembly_gsnap_hitsummary2_tab,
      assembly_refined_gsnap_counts_with_dcr_json = BlastContigs_refined_gsnap_out.assembly_refined_gsnap_counts_with_dcr_json,
      assembly_gsnap_contig_summary_json = BlastContigs_refined_gsnap_out.assembly_gsnap_contig_summary_json,
      assembly_gsnap_blast_top_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_top_m8,
      assembly_rapsearch2_blast_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_m8,
      assembly_rapsearch2_reassigned_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_reassigned_m8,
      assembly_rapsearch2_hitsummary2_tab = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_hitsummary2_tab,
      assembly_refined_rapsearch2_counts_with_dcr_json = BlastContigs_refined_rapsearch2_out.assembly_refined_rapsearch2_counts_with_dcr_json,
      assembly_rapsearch2_contig_summary_json = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_contig_summary_json,
      assembly_rapsearch2_blast_top_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_top_m8
  }

  call GenerateAnnotatedFasta {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      assembly_gsnap_blast_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_m8,
      assembly_gsnap_reassigned_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_reassigned_m8,
      assembly_gsnap_hitsummary2_tab = BlastContigs_refined_gsnap_out.assembly_gsnap_hitsummary2_tab,
      assembly_refined_gsnap_counts_with_dcr_json = BlastContigs_refined_gsnap_out.assembly_refined_gsnap_counts_with_dcr_json,
      assembly_gsnap_contig_summary_json = BlastContigs_refined_gsnap_out.assembly_gsnap_contig_summary_json,
      assembly_gsnap_blast_top_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_top_m8,
      assembly_rapsearch2_blast_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_m8,
      assembly_rapsearch2_reassigned_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_reassigned_m8,
      assembly_rapsearch2_hitsummary2_tab = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_hitsummary2_tab,
      assembly_refined_rapsearch2_counts_with_dcr_json = BlastContigs_refined_rapsearch2_out.assembly_refined_rapsearch2_counts_with_dcr_json,
      assembly_rapsearch2_contig_summary_json = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_contig_summary_json,
      assembly_rapsearch2_blast_top_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_top_m8,
      cdhitdup_out_dedup1_fa_clstr = cdhitdup_out_dedup1_fa_clstr,
      cdhitdup_out_dedup1_fa = cdhitdup_out_dedup1_fa,
      cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv
  }

  call GenerateTaxidFasta {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      assembly_refined_annotated_merged_fa = GenerateAnnotatedFasta.assembly_refined_annotated_merged_fa,
      assembly_refined_unidentified_fa = GenerateAnnotatedFasta.assembly_refined_unidentified_fa,
      assembly_gsnap_blast_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_m8,
      assembly_gsnap_reassigned_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_reassigned_m8,
      assembly_gsnap_hitsummary2_tab = BlastContigs_refined_gsnap_out.assembly_gsnap_hitsummary2_tab,
      assembly_refined_gsnap_counts_with_dcr_json = BlastContigs_refined_gsnap_out.assembly_refined_gsnap_counts_with_dcr_json,
      assembly_gsnap_contig_summary_json = BlastContigs_refined_gsnap_out.assembly_gsnap_contig_summary_json,
      assembly_gsnap_blast_top_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_top_m8,
      assembly_rapsearch2_blast_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_m8,
      assembly_rapsearch2_reassigned_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_reassigned_m8,
      assembly_rapsearch2_hitsummary2_tab = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_hitsummary2_tab,
      assembly_refined_rapsearch2_counts_with_dcr_json = BlastContigs_refined_rapsearch2_out.assembly_refined_rapsearch2_counts_with_dcr_json,
      assembly_rapsearch2_contig_summary_json = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_contig_summary_json,
      assembly_rapsearch2_blast_top_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_top_m8,
      lineage_db = lineage_db
  }

  call GenerateTaxidLocator {
    input:
      docker_image_id = docker_image_id,
      dag_branch = dag_branch,
      s3_wd_uri = s3_wd_uri,
      assembly_refined_taxid_annot_fasta = GenerateTaxidFasta.assembly_refined_taxid_annot_fasta
  }

  output {
    File assembly_out_assembly_contigs_fasta = RunAssembly.assembly_contigs_fasta
    File assembly_out_assembly_scaffolds_fasta = RunAssembly.assembly_scaffolds_fasta
    File assembly_out_assembly_read_contig_sam = RunAssembly.assembly_read_contig_sam
    File assembly_out_assembly_contig_stats_json = RunAssembly.assembly_contig_stats_json
    File? assembly_out_count = RunAssembly.output_read_count
    File coverage_out_assembly_contig_coverage_json = GenerateCoverageStats.assembly_contig_coverage_json
    File coverage_out_assembly_contig_coverage_summary_csv = GenerateCoverageStats.assembly_contig_coverage_summary_csv
    File? coverage_out_count = GenerateCoverageStats.output_read_count
    File gsnap_accessions_out_assembly_nt_refseq_fasta = DownloadAccessions_gsnap_accessions_out.assembly_nt_refseq_fasta
    File? gsnap_accessions_out_count = DownloadAccessions_gsnap_accessions_out.output_read_count
    File rapsearch2_accessions_out_assembly_nr_refseq_fasta = DownloadAccessions_rapsearch2_accessions_out.assembly_nr_refseq_fasta
    File? rapsearch2_accessions_out_count = DownloadAccessions_rapsearch2_accessions_out.output_read_count
    File refined_gsnap_out_assembly_gsnap_blast_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_m8
    File refined_gsnap_out_assembly_gsnap_reassigned_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_reassigned_m8
    File refined_gsnap_out_assembly_gsnap_hitsummary2_tab = BlastContigs_refined_gsnap_out.assembly_gsnap_hitsummary2_tab
    File refined_gsnap_out_assembly_refined_gsnap_counts_with_dcr_json = BlastContigs_refined_gsnap_out.assembly_refined_gsnap_counts_with_dcr_json
    File refined_gsnap_out_assembly_gsnap_contig_summary_json = BlastContigs_refined_gsnap_out.assembly_gsnap_contig_summary_json
    File refined_gsnap_out_assembly_gsnap_blast_top_m8 = BlastContigs_refined_gsnap_out.assembly_gsnap_blast_top_m8
    File? refined_gsnap_out_count = BlastContigs_refined_gsnap_out.output_read_count
    File refined_rapsearch2_out_assembly_rapsearch2_blast_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_m8
    File refined_rapsearch2_out_assembly_rapsearch2_reassigned_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_reassigned_m8
    File refined_rapsearch2_out_assembly_rapsearch2_hitsummary2_tab = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_hitsummary2_tab
    File refined_rapsearch2_out_assembly_refined_rapsearch2_counts_with_dcr_json = BlastContigs_refined_rapsearch2_out.assembly_refined_rapsearch2_counts_with_dcr_json
    File refined_rapsearch2_out_assembly_rapsearch2_contig_summary_json = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_contig_summary_json
    File refined_rapsearch2_out_assembly_rapsearch2_blast_top_m8 = BlastContigs_refined_rapsearch2_out.assembly_rapsearch2_blast_top_m8
    File? refined_rapsearch2_out_count = BlastContigs_refined_rapsearch2_out.output_read_count
    File refined_taxon_count_out_assembly_refined_taxon_counts_with_dcr_json = CombineTaxonCounts.assembly_refined_taxon_counts_with_dcr_json
    File? refined_taxon_count_out_count = CombineTaxonCounts.output_read_count
    File contig_summary_out_assembly_combined_contig_summary_json = CombineJson.assembly_combined_contig_summary_json
    File? contig_summary_out_count = CombineJson.output_read_count
    File refined_annotated_out_assembly_refined_annotated_merged_fa = GenerateAnnotatedFasta.assembly_refined_annotated_merged_fa
    File refined_annotated_out_assembly_refined_unidentified_fa = GenerateAnnotatedFasta.assembly_refined_unidentified_fa
    File? refined_annotated_out_count = GenerateAnnotatedFasta.output_read_count
    File refined_taxid_fasta_out_assembly_refined_taxid_annot_fasta = GenerateTaxidFasta.assembly_refined_taxid_annot_fasta
    File? refined_taxid_fasta_out_count = GenerateTaxidFasta.output_read_count
    File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_nt_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_nt_fasta
    File refined_taxid_locator_out_assembly_refined_taxid_locations_nt_json = GenerateTaxidLocator.assembly_refined_taxid_locations_nt_json
    File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_nr_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_nr_fasta
    File refined_taxid_locator_out_assembly_refined_taxid_locations_nr_json = GenerateTaxidLocator.assembly_refined_taxid_locations_nr_json
    File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_genus_nt_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_genus_nt_fasta
    File refined_taxid_locator_out_assembly_refined_taxid_locations_genus_nt_json = GenerateTaxidLocator.assembly_refined_taxid_locations_genus_nt_json
    File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_genus_nr_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_genus_nr_fasta
    File refined_taxid_locator_out_assembly_refined_taxid_locations_genus_nr_json = GenerateTaxidLocator.assembly_refined_taxid_locations_genus_nr_json
    File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_family_nt_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_family_nt_fasta
    File refined_taxid_locator_out_assembly_refined_taxid_locations_family_nt_json = GenerateTaxidLocator.assembly_refined_taxid_locations_family_nt_json
    File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_family_nr_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_family_nr_fasta
    File refined_taxid_locator_out_assembly_refined_taxid_locations_family_nr_json = GenerateTaxidLocator.assembly_refined_taxid_locations_family_nr_json
    File refined_taxid_locator_out_assembly_refined_taxid_locations_combined_json = GenerateTaxidLocator.assembly_refined_taxid_locations_combined_json
    File? refined_taxid_locator_out_count = GenerateTaxidLocator.output_read_count
  }
}
