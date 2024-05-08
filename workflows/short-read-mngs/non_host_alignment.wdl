version 1.0

task CombineTaxonCounts {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] counts_json_files
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.combine_taxon_counts \
    --step-class PipelineStepCombineTaxonCounts \
    --step-name taxon_count_out \
    --input-files '["~{sep='", "' counts_json_files}"]' \
    --output-files '["taxon_counts_with_dcr.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    String step_description_md = read_string("taxon_count_out.description.md")
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
    File czid_dedup_out_duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name non_host_alignment \
    --step-module idseq_dag.steps.generate_annotated_fasta \
    --step-class PipelineStepGenerateAnnotatedFasta \
    --step-name annotated_out \
    --input-files '[["~{sep='","' host_filter_out_gsnap_filter_fa}"], ["~{gsnap_m8}", "~{gsnap_deduped_m8}", "~{gsnap_hitsummary_tab}", "~{gsnap_counts_with_dcr_json}"], ["~{rapsearch2_m8}", "~{rapsearch2_deduped_m8}", "~{rapsearch2_hitsummary_tab}", "~{rapsearch2_counts_with_dcr_json}"], ["~{czid_dedup_out_duplicate_clusters_csv}"], ["~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '["annotated_merged.fa", "unidentified.fa"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  >>>
  output {
    String step_description_md = read_string("annotated_out.description.md")
    File annotated_merged_fa = "annotated_merged.fa"
    File unidentified_fa = "unidentified.fa"
    File? output_read_count = "annotated_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunAlignment_minimap2_out {
    input {
        String s3_wd_uri
        Array[File]+ fastas
        String db_path
        String minimap2_args 
        String docker_image_id
        Boolean? run_locally = false
        File? local_minimap2_index 
        String prefix
        String minimap2_wdl_version
    }

    command <<<
        # --step-name minimap2_out
        set -euxo pipefail

        if [[ "~{run_locally}" == true ]]; then
          minimap2-scatter ~{minimap2_args} "~{local_minimap2_index}" "~{sep=' ' fastas}" > "~{prefix}.paf"
        else
          python3 <<CODE
        import os
        from idseq_utils.batch_run_helpers import run_alignment

        run_alignment(
            input_dir="~{s3_wd_uri}",
            db_path="~{db_path}",
            result_path="gsnap.paf",
            aligner="minimap2",
            aligner_args="~{minimap2_args}",
            aligner_wdl_version="~{minimap2_wdl_version}",
            queries=["~{sep='", "' fastas}"],
        )
        CODE
        fi
        python3 /usr/local/lib/python3.6/dist-packages/idseq_utils/paf2blast6.py gsnap.paf
        mv *frompaf.m8 "gsnap.m8" # TODO: rewrite paf2blast6.py to output in this format
        minimap2-scatter --version > minimap2_version.txt
    >>>
    output {
        File out_paf = "gsnap.paf"
        File out_m8 = "gsnap.m8"
        File? version = "minimap2_version.txt"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunAlignment_diamond_out {
    input {
        String docker_image_id
        String s3_wd_uri
        Array[File]+ fastas
        String db_path
        String diamond_args 
        Boolean? run_locally = false
        File? local_diamond_index
        String prefix
        String diamond_wdl_version
    }

    command <<<
        # --step-name diamond_out
        set -euxo pipefail  
        if [[ "~{run_locally}" == true ]]; then 
          diamond makedb --in "~{local_diamond_index}" -d reference
          diamond blastx -d reference -q "~{sep=' ' fastas}" -o "rapsearch2.m8" "--~{diamond_args}"
        else
          python3 <<CODE
        import os 
        from idseq_utils.batch_run_helpers import run_alignment

        run_alignment(
            input_dir="~{s3_wd_uri}",
            db_path="~{db_path}",
            result_path="rapsearch2.m8",
            aligner="diamond",
            aligner_args="~{diamond_args}",
            aligner_wdl_version="~{diamond_wdl_version}",
            queries=["~{sep='", "' fastas}"],
        )
        CODE
        diamond --version > diamond_version.txt
        fi
    >>>

    output {
        File out_m8 = "rapsearch2.m8"
        File? version = "diamond_version.txt"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunCallHitsMinimap2 {
    input {
        File m8_file
        File lineage_db
        File taxon_blacklist
        File deuterostome_db
        File accession2taxid
        File duplicate_cluster_size
        String prefix 
        Int min_read_length = 36
        String docker_image_id
        String count_type = "NT"
        String s3_wd_uri
    }

    command <<<
        # --step-name minimap2_call_hits_out
        set -euxo pipefail
        python3 <<CODE
        from idseq_dag.util.m8 import call_hits_m8, generate_taxon_count_json_from_m8
        call_hits_m8(
            input_m8="~{m8_file}",
            lineage_map_path="~{lineage_db}",
            accession2taxid_dict_path="~{accession2taxid}",
            output_m8="gsnap.deduped.m8",
            output_summary="gsnap.hitsummary.tab",
            min_alignment_length=~{min_read_length},
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
        )
        generate_taxon_count_json_from_m8(
            blastn_6_path="gsnap.deduped.m8",
            hit_level_path="gsnap.hitsummary.tab",
            count_type="~{count_type}",
            lineage_map_path="~{lineage_db}",
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
            duplicate_cluster_sizes_path="~{duplicate_cluster_size}",
            output_json_file="gsnap_counts_with_dcr.json",
        )
        CODE
        >>>

    output {
        File deduped_out_m8 = "gsnap.deduped.m8"
        File hitsummary = "gsnap.hitsummary.tab"
        File counts_json = "gsnap_counts_with_dcr.json"
        File? output_read_count = "minimap2_output.count"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunCallHitsDiamond {
    input {
        File m8_file
        File lineage_db
        File taxon_blacklist
        File deuterostome_db
        File accession2taxid
        File duplicate_cluster_size
        String prefix 
        Int min_read_length = 0
        String docker_image_id
        String count_type = "NR"
        String s3_wd_uri
    }
    command <<<
        # --step-name diamond_call_hits_out
        set -euxo pipefail
        python3 <<CODE
        from idseq_dag.util.m8 import call_hits_m8, generate_taxon_count_json_from_m8
        call_hits_m8(
            input_m8="~{m8_file}",
            lineage_map_path="~{lineage_db}",
            accession2taxid_dict_path="~{accession2taxid}",
            output_m8="rapsearch2.deduped.m8",
            output_summary="rapsearch2.hitsummary.tab",
            min_alignment_length=~{min_read_length},
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
        )
        generate_taxon_count_json_from_m8(
            blastn_6_path="rapsearch2.deduped.m8",
            hit_level_path="rapsearch2.hitsummary.tab",
            count_type="~{count_type}",
            lineage_map_path="~{lineage_db}",
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
            duplicate_cluster_sizes_path="~{duplicate_cluster_size}",
            output_json_file="rapsearch2_counts_with_dcr.json",
        )
        CODE
        >>>

    output {
        File deduped_out_m8 = "rapsearch2.deduped.m8"
        File hitsummary = "rapsearch2.hitsummary.tab"
        File counts_json = "rapsearch2_counts_with_dcr.json"
        File? output_read_count = "minimap2_output.count"
    }

    runtime {
        docker: docker_image_id
    }
}

workflow czid_non_host_alignment {
  input {
    String docker_image_id
    String s3_wd_uri
    File host_filter_out_gsnap_filter_1_fa
    File? host_filter_out_gsnap_filter_2_fa
    File? host_filter_out_gsnap_filter_merged_fa
    File duplicate_cluster_sizes_tsv
    File czid_dedup_out_duplicate_clusters_csv
    String index_version = "2021-01-22"
    File lineage_db = "s3://czid-public-references/taxonomy/2021-01-22/taxid-lineages.db"
    File accession2taxid_db = "s3://czid-public-references/alignment_data/2021-01-22/accession2taxid.db"
    File taxon_blacklist = "s3://czid-public-references/taxonomy/2021-01-22/taxon_blacklist.txt"
    String index_dir_suffix = index_version
    Int min_read_length = 36
    File deuterostome_db = "s3://czid-public-references/taxonomy/2021-01-22/deuterostome_taxids.txt"
    Boolean use_deuterostome_filter = true
    Boolean use_taxon_whitelist = false
    Boolean alignment_scalability = false
    File? local_gsnap_index
    File? minimap2_local_db_path
    File? diamond_local_db_path
    String? local_gsnap_genome_name
    File? local_rapsearch2_index
    String minimap2_db = "s3://czid-public-references/minimap2-test/2021-01-22/nt_k14_w8_20/"
    String diamond_db = "s3://czid-public-references/diamond-test/2021-01-22/"
    String minimap2_args = "-cx sr --secondary=yes"
    String diamond_args = "mid-sensitive"
    String minimap2_prefix = "gsnap"
    String diamond_prefix = "rapsearch2"
    String minimap2_wdl_version = "v1.0.0"
    String diamond_wdl_version = "v1.0.0"

  }
  call RunAlignment_minimap2_out { 
    input:         
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      fastas = [select_first([host_filter_out_gsnap_filter_merged_fa, host_filter_out_gsnap_filter_1_fa])], #select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa]),
      db_path = minimap2_db,
      minimap2_args = minimap2_args,
      run_locally = defined(minimap2_local_db_path),
      local_minimap2_index = minimap2_local_db_path,
      prefix= minimap2_prefix,
      minimap2_wdl_version=minimap2_wdl_version
  }
  call RunCallHitsMinimap2{ 
    input:
      m8_file = RunAlignment_minimap2_out.out_m8,
      lineage_db = lineage_db,
      duplicate_cluster_size = duplicate_cluster_sizes_tsv,
      taxon_blacklist = taxon_blacklist,
      deuterostome_db = deuterostome_db,
      accession2taxid = accession2taxid_db,
      prefix = minimap2_prefix,
      min_read_length = min_read_length,
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
  }
  call RunAlignment_diamond_out {
    input: 
      fastas = [select_first([host_filter_out_gsnap_filter_merged_fa, host_filter_out_gsnap_filter_1_fa])], #select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa]),
      s3_wd_uri = s3_wd_uri,
      db_path = diamond_db,
      diamond_args = diamond_args,
      prefix = diamond_prefix,
      run_locally = defined(diamond_local_db_path),
      local_diamond_index = diamond_local_db_path,
      diamond_wdl_version=diamond_wdl_version,
      docker_image_id = docker_image_id
  }
  call RunCallHitsDiamond { 
    input:
      m8_file = RunAlignment_diamond_out.out_m8,
      lineage_db = lineage_db,
      duplicate_cluster_size = duplicate_cluster_sizes_tsv,
      taxon_blacklist = taxon_blacklist,
      deuterostome_db = deuterostome_db,
      accession2taxid = accession2taxid_db,
      prefix = diamond_prefix,
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
  }

  call CombineTaxonCounts {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      counts_json_files = [
        RunCallHitsMinimap2.counts_json,
        RunCallHitsDiamond.counts_json
      ]
  }

  call GenerateAnnotatedFasta {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      gsnap_m8 = RunAlignment_minimap2_out.out_m8,
      gsnap_deduped_m8 = RunCallHitsMinimap2.deduped_out_m8,
      gsnap_hitsummary_tab = RunCallHitsMinimap2.hitsummary,
      gsnap_counts_with_dcr_json = RunCallHitsMinimap2.counts_json,
      rapsearch2_m8 = RunAlignment_diamond_out.out_m8,
      rapsearch2_deduped_m8 = RunCallHitsDiamond.deduped_out_m8,
      rapsearch2_hitsummary_tab = RunCallHitsDiamond.hitsummary,
      rapsearch2_counts_with_dcr_json = RunCallHitsDiamond.counts_json,
      czid_dedup_out_duplicate_clusters_csv = czid_dedup_out_duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = duplicate_cluster_sizes_tsv
  }

  output {
    File gsnap_out_gsnap_m8 = RunAlignment_minimap2_out.out_m8
    File gsnap_out_gsnap_deduped_m8 = RunCallHitsMinimap2.deduped_out_m8
    File gsnap_out_gsnap_hitsummary_tab = RunCallHitsMinimap2.hitsummary
    File gsnap_out_gsnap_counts_with_dcr_json = RunCallHitsMinimap2.counts_json
    File? minimap2_version = RunAlignment_minimap2_out.version
    File? gsnap_out_count = RunCallHitsMinimap2.output_read_count
    File rapsearch2_out_rapsearch2_m8 = RunAlignment_diamond_out.out_m8
    File rapsearch2_out_rapsearch2_deduped_m8 = RunCallHitsDiamond.deduped_out_m8
    File rapsearch2_out_rapsearch2_hitsummary_tab = RunCallHitsDiamond.hitsummary
    File rapsearch2_out_rapsearch2_counts_with_dcr_json = RunCallHitsDiamond.counts_json
    File? diamond_version = RunAlignment_diamond_out.version
    File? rapsearch2_out_count = RunCallHitsDiamond.output_read_count
    File taxon_count_out_taxon_counts_with_dcr_json = CombineTaxonCounts.taxon_counts_with_dcr_json
    File? taxon_count_out_count = CombineTaxonCounts.output_read_count
    File annotated_out_annotated_merged_fa = GenerateAnnotatedFasta.annotated_merged_fa
    File annotated_out_unidentified_fa = GenerateAnnotatedFasta.unidentified_fa
    File? annotated_out_count = GenerateAnnotatedFasta.output_read_count
  }
}
