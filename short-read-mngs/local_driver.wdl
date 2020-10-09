version 1.0
# local_driver.wdl: this top-level workflow runs the four stages (host_filter, non_host_alignment,
# postprocess, experimental) in sequence. The IDseq back-end invokes those four WDLs separately
# for various reasons, which is effectively the same as running this locally.

import "host_filter.wdl" as stage1
import "non_host_alignment.wdl" as stage2
import "postprocess.wdl" as stage3
import "experimental.wdl" as stage4

workflow idseq_short_read_mngs {
    input {
        String docker_image_id
        File fastqs_0
        File? fastqs_1
        File non_host_rapsearch2_index
        File non_host_gsnap_index
        String non_host_gsnap_genome_name = "nt_k16"
        String s3_wd_uri = ""
    }
    call stage1.idseq_host_filter as host_filter {
        input:
        fastqs_0 = fastqs_0,
        fastqs_1 = fastqs_1,
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri
    }
    call stage2.idseq_non_host_alignment as non_host_alignment {
        input:
        host_filter_out_gsnap_filter_1_fa = host_filter.gsnap_filter_out_gsnap_filter_1_fa,
        host_filter_out_gsnap_filter_2_fa = host_filter.gsnap_filter_out_gsnap_filter_2_fa,
        host_filter_out_gsnap_filter_merged_fa = host_filter.gsnap_filter_out_gsnap_filter_merged_fa,
        duplicate_cluster_sizes_tsv = host_filter.idseq_dedup_out_duplicate_cluster_sizes_tsv,
        idseq_dedup_out_duplicate_clusters_csv = host_filter.idseq_dedup_out_duplicate_clusters_csv,
        local_gsnap_index = non_host_gsnap_index,
        local_gsnap_genome_name = non_host_gsnap_genome_name,
        local_rapsearch2_index = non_host_rapsearch2_index,
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri
    }
    call stage3.idseq_postprocess as postprocess {
        input:
        host_filter_out_gsnap_filter_1_fa = host_filter.gsnap_filter_out_gsnap_filter_1_fa,
        host_filter_out_gsnap_filter_2_fa = host_filter.gsnap_filter_out_gsnap_filter_2_fa,
        host_filter_out_gsnap_filter_merged_fa = host_filter.gsnap_filter_out_gsnap_filter_merged_fa,
        duplicate_cluster_sizes_tsv = host_filter.idseq_dedup_out_duplicate_cluster_sizes_tsv,
        idseq_dedup_out_duplicate_clusters_csv = host_filter.idseq_dedup_out_duplicate_clusters_csv,
        gsnap_out_gsnap_m8 = non_host_alignment.gsnap_out_gsnap_m8,
        gsnap_out_gsnap_deduped_m8 = non_host_alignment.gsnap_out_gsnap_deduped_m8,
        gsnap_out_gsnap_hitsummary_tab = non_host_alignment.gsnap_out_gsnap_hitsummary_tab,
        gsnap_out_gsnap_counts_with_dcr_json = non_host_alignment.gsnap_out_gsnap_counts_with_dcr_json,
        rapsearch2_out_rapsearch2_m8 = non_host_alignment.rapsearch2_out_rapsearch2_m8,
        rapsearch2_out_rapsearch2_deduped_m8 = non_host_alignment.rapsearch2_out_rapsearch2_deduped_m8,
        rapsearch2_out_rapsearch2_hitsummary_tab = non_host_alignment.rapsearch2_out_rapsearch2_hitsummary_tab,
        rapsearch2_out_rapsearch2_counts_with_dcr_json = non_host_alignment.rapsearch2_out_rapsearch2_counts_with_dcr_json,
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri
    }
    call stage4.idseq_experimental as experimental {
        input:
        taxid_fasta_in_annotated_merged_fa = non_host_alignment.annotated_out_annotated_merged_fa,
        taxid_fasta_in_gsnap_hitsummary_tab = non_host_alignment.gsnap_out_gsnap_hitsummary_tab,
        taxid_fasta_in_rapsearch2_hitsummary_tab = non_host_alignment.rapsearch2_out_rapsearch2_hitsummary_tab,
        gsnap_m8_gsnap_deduped_m8 = non_host_alignment.gsnap_out_gsnap_m8,
        refined_gsnap_in_gsnap_reassigned_m8 = postprocess.refined_gsnap_out_assembly_gsnap_reassigned_m8,
        refined_gsnap_in_gsnap_hitsummary2_tab = postprocess.refined_gsnap_out_assembly_gsnap_hitsummary2_tab,
        refined_gsnap_in_gsnap_blast_top_m8 = postprocess.refined_gsnap_out_assembly_gsnap_blast_top_m8,
        contig_in_contig_coverage_json = postprocess.coverage_out_assembly_contig_coverage_json,
        contig_in_contig_stats_json = postprocess.assembly_out_assembly_contig_stats_json,
        contig_in_contigs_fasta = postprocess.assembly_out_assembly_contigs_fasta,
        fastqs_0 = fastqs_0,
        fastqs_1 = fastqs_1,
        nonhost_fasta_refined_taxid_annot_fasta = postprocess.refined_taxid_fasta_out_assembly_refined_taxid_annot_fasta,
        duplicate_clusters_csv = host_filter.idseq_dedup_out_duplicate_clusters_csv,
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri
    }
}
