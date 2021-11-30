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
    String step_description_md = read_string("gsnap_out.description.md")
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
    String step_description_md = read_string("rapsearch2_out.description.md")
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
        Array[File]+ fastqs
        String input_dir
        String chunk_dir
        String db_path
        String minimap2_args 
        String docker_image_id
        String prefix
    }

    command <<<
        set -euxo pipefail
        echo STARTING
        export DEPLOYMENT_ENVIRONMENT=dev
        export AWS_REGION="us-west-2"
        export AWS_DEFAULT_REGION="us-west-2"
        for fastq in ~{sep=' ' fastqs}; do 
          s3parcp $fastq "~{input_dir}"
        done 
        python3 <<CODE
        from idseq_utils.run_minimap2 import run_minimap2
        fastqs = ["~{sep='", "' fastqs}"]
        run_minimap2("~{input_dir}", "~{chunk_dir}", "~{db_path}", "~{prefix}.paf", "~{minimap2_args}", *fastqs)
        CODE
        python3 /usr/local/lib/python3.6/dist-packages/idseq_utils/paf2blast6.py "~{prefix}".paf
        mv *frompaf.m8 "~{prefix}.m8" # TODO: rewrite paf2blast6.py to output in this format
    >>>
    output {
        File out_paf = "~{prefix}.paf"
        File out_m8 = "~{prefix}.m8"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunAlignment_diamond_out {
    input {
        Array[File]+ fastqs
        String input_dir
        String chunk_dir
        String db_path
        String diamond_args 
        String prefix
        String docker_image_id
    }

    command <<<
        set -euxo pipefail  
        echo STARTING
        export DEPLOYMENT_ENVIRONMENT=dev
        export AWS_REGION="us-west-2"
        export AWS_DEFAULT_REGION="us-west-2"
        for fastq in ~{sep=' ' fastqs}; do 
          s3parcp $fastq "~{input_dir}"
        done 
        python3 <<CODE
        from idseq_utils.run_diamond import run_diamond
        fastqs = ["~{sep='", "' fastqs}"]
        run_diamond("~{input_dir}", "~{chunk_dir}", "~{db_path}", "~{prefix}.m8", *fastqs)
        CODE
    >>>

    output {
        File out_m8 = "~{prefix}.m8"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunCallHits {
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
    }

    command <<<
        set -euxo pipefail
        python3 <<CODE
        from idseq_dag.util.m8 import call_hits_m8, generate_taxon_count_json_from_m8
        call_hits_m8("~{m8_file}", "~{lineage_db}", "~{accession2taxid}", "~{prefix}.deduped.m8", "~{prefix}.hitsummary.tab", ~{min_read_length}, "~{deuterostome_db}", None, "~{taxon_blacklist}")
        generate_taxon_count_json_from_m8( "~{prefix}.deduped.m8", "~{prefix}.hitsummary.tab", "NT", "~{lineage_db}", "~{deuterostome_db}", None, "~{taxon_blacklist}", "~{duplicate_cluster_size}", "~{prefix}_counts_with_dcr.json")
        CODE
        >>>

    output {
        File deduped_out_m8 = "~{prefix}.deduped.m8"
        File hitsummary = "~{prefix}.hitsummary.tab"
        File counts_json = "~{prefix}_counts_with_dcr.json"
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
    String index_version = "2021-01-22"
    File lineage_db = "s3://idseq-public-references/taxonomy/2021-01-22/taxid-lineages.db"
    File accession2taxid_db = "s3://idseq-public-references/alignment_data/2021-01-22/accession2taxid.db"
    File taxon_blacklist = "s3://idseq-public-references/taxonomy/2021-01-22/taxon_blacklist.txt"
    String index_dir_suffix = index_version
    Int min_read_length = 36
    File deuterostome_db = "s3://idseq-public-references/taxonomy/2021-01-22/deuterostome_taxids.txt"
    Boolean use_deuterostome_filter = true
    Boolean use_taxon_whitelist = false
    Boolean alignment_scalability = true
    File? local_gsnap_index
    String? local_gsnap_genome_name
    File? local_rapsearch2_index
    String alignment_input_dir = "s3://idseq-samples-development/samples/alignment-scalability-test/combined-test/1/"
    String minimap2_chunk_dir = "s3://idseq-samples-development/samples/alignment-scalability-test/combined-test/1/minimap2-chunks/"
    String diamond_chunk_dir = "s3://idseq-samples-development/samples/alignment-scalability-test/combined-test/1/diamond-chunks/"
    String minimap2_db = "s3://idseq-public-references/minimap2-test/2021-01-22/nt_k12_w5_20/"
    String diamond_db = "s3://idseq-public-references/diamond-test/2021-01-22/"
    String minimap2_args = "-cx sr --secondary=yes"
    String diamond_args = ""
    String minimap2_prefix = "gsnap"
    String diamond_prefix = "rapsearch2"

  }
  if (!alignment_scalability){
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
  }
  if (alignment_scalability) { 
    call RunAlignment_minimap2_out { 
      input:         
        fastqs = [select_first([host_filter_out_gsnap_filter_merged_fa, host_filter_out_gsnap_filter_1_fa])], #select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa]),
        input_dir = alignment_input_dir,
        chunk_dir = minimap2_chunk_dir,
        db_path = minimap2_db,
        minimap2_args = minimap2_args,
        prefix= minimap2_prefix,
        docker_image_id = docker_image_id
    }
    call RunCallHits { 
        input:
        m8_file = RunAlignment_minimap2_out.out_m8,
        lineage_db = lineage_db,
        duplicate_cluster_size = duplicate_cluster_sizes_tsv,
        taxon_blacklist = taxon_blacklist,
        deuterostome_db = deuterostome_db,
        accession2taxid = accession2taxid_db,
        prefix = minimap2_prefix,
        min_read_length = min_read_length,
        docker_image_id = docker_image_id
    }
    call RunAlignment_diamond_out {
      input: 
      fastqs = [select_first([host_filter_out_gsnap_filter_merged_fa, host_filter_out_gsnap_filter_1_fa])], #select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa]),
      input_dir = alignment_input_dir,
      chunk_dir = diamond_chunk_dir,
      db_path = diamond_db,
      diamond_args = diamond_args,
      prefix = diamond_prefix,
      docker_image_id = docker_image_id
    }
    call RunCallHits as RunCallHitsDiamond { 
        input:
        m8_file = RunAlignment_diamond_out.out_m8,
        lineage_db = lineage_db,
        duplicate_cluster_size = duplicate_cluster_sizes_tsv,
        taxon_blacklist = taxon_blacklist,
        deuterostome_db = deuterostome_db,
        accession2taxid = accession2taxid_db,
        prefix = diamond_prefix,
        min_read_length = 0,
        docker_image_id = docker_image_id
    }
  }



  call CombineTaxonCounts {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      counts_json_files = [
        select_first([RunAlignment_gsnap_out.gsnap_counts_with_dcr_json, RunCallHits.counts_json]),
        select_first([RunAlignment_rapsearch2_out.rapsearch2_counts_with_dcr_json, RunCallHitsDiamond.counts_json])
      ]
  }

  call GenerateAnnotatedFasta {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      host_filter_out_gsnap_filter_fa = select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa, host_filter_out_gsnap_filter_merged_fa]),
      gsnap_m8 = select_first([RunAlignment_gsnap_out.gsnap_m8, RunAlignment_minimap2_out.out_m8]),
      gsnap_deduped_m8 = select_first([RunAlignment_gsnap_out.gsnap_deduped_m8, RunCallHits.deduped_out_m8]),
      gsnap_hitsummary_tab = select_first([RunAlignment_gsnap_out.gsnap_hitsummary_tab, RunCallHits.hitsummary]),
      gsnap_counts_with_dcr_json = select_first([RunAlignment_gsnap_out.gsnap_counts_with_dcr_json, RunCallHits.counts_json]),
      rapsearch2_m8 = select_first([RunAlignment_rapsearch2_out.rapsearch2_m8, RunAlignment_diamond_out.out_m8]),
      rapsearch2_deduped_m8 = select_first([RunAlignment_rapsearch2_out.rapsearch2_deduped_m8, RunCallHitsDiamond.deduped_out_m8]),
      rapsearch2_hitsummary_tab = select_first([RunAlignment_rapsearch2_out.rapsearch2_hitsummary_tab, RunCallHitsDiamond.hitsummary]),
      rapsearch2_counts_with_dcr_json = select_first([RunAlignment_rapsearch2_out.rapsearch2_counts_with_dcr_json, RunCallHitsDiamond.counts_json]),
      idseq_dedup_out_duplicate_clusters_csv = idseq_dedup_out_duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = duplicate_cluster_sizes_tsv
  }

  output {
    File gsnap_out_gsnap_m8 = select_first([RunAlignment_gsnap_out.gsnap_m8, RunAlignment_minimap2_out.out_m8])
    File gsnap_out_gsnap_deduped_m8 = select_first([RunAlignment_gsnap_out.gsnap_deduped_m8, RunCallHits.deduped_out_m8])
    File gsnap_out_gsnap_hitsummary_tab = select_first([RunAlignment_gsnap_out.gsnap_hitsummary_tab, RunCallHits.hitsummary])
    File gsnap_out_gsnap_counts_with_dcr_json = select_first([RunAlignment_gsnap_out.gsnap_counts_with_dcr_json, RunCallHits.counts_json])
    File? gsnap_out_count = RunAlignment_gsnap_out.output_read_count
    File rapsearch2_out_rapsearch2_m8 = select_first([RunAlignment_rapsearch2_out.rapsearch2_m8, RunAlignment_diamond_out.out_m8])
    File rapsearch2_out_rapsearch2_deduped_m8 = select_first([RunAlignment_rapsearch2_out.rapsearch2_deduped_m8, RunCallHitsDiamond.deduped_out_m8])
    File rapsearch2_out_rapsearch2_hitsummary_tab = select_first([RunAlignment_rapsearch2_out.rapsearch2_hitsummary_tab, RunCallHitsDiamond.hitsummary])
    File rapsearch2_out_rapsearch2_counts_with_dcr_json = select_first([RunAlignment_rapsearch2_out.rapsearch2_counts_with_dcr_json, RunCallHitsDiamond.counts_json])
    File? rapsearch2_out_count = RunAlignment_rapsearch2_out.output_read_count
    File taxon_count_out_taxon_counts_with_dcr_json = CombineTaxonCounts.taxon_counts_with_dcr_json
    File? taxon_count_out_count = CombineTaxonCounts.output_read_count
    File annotated_out_annotated_merged_fa = GenerateAnnotatedFasta.annotated_merged_fa
    File annotated_out_unidentified_fa = GenerateAnnotatedFasta.unidentified_fa
    File? annotated_out_count = GenerateAnnotatedFasta.output_read_count
  }
}
