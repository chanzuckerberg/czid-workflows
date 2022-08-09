version 1.0

task DownloadFastas {
    input {
        Array[File]+ fastas
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        for fasta in ~{sep=' ' fastas}
        do
          cp "$fasta" "$(basename $fasta)"
        done
    >>>
    output {
        Array[File]+ loaded_fastas = glob("*.fasta")
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
        Boolean run_locally = false
        File? local_minimap2_index 
        String prefix
    }

    command <<<
        # --step-name minimap2_out
        set -euxo pipefail

        if [[ "~{run_locally}" == true ]]; then
          minimap2 ~{minimap2_args} "~{local_minimap2_index}" "~{sep=' ' fastas}" > "~{prefix}.paf"
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
            queries=["~{sep='", "' fastas}"],
        )
        CODE
        fi
        python3 /usr/local/lib/python3.6/dist-packages/idseq_utils/paf2blast6.py gsnap.paf
        mv *frompaf.m8 "gsnap.m8" # TODO: rewrite paf2blast6.py to output in this format
        minimap2 --version > minimap2_version.txt
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

task RunCallHitsMinimap2 {
    input {
        File m8_file
        File lineage_db
        File taxon_blacklist
        File deuterostome_db
        File accession2taxid
        Int min_read_length = 36
        String docker_image_id
        String count_type = "NT"
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

workflow czid_non_host_alignment {
  input {
    String docker_image_id
    String s3_wd_uri
    File host_filter_out_gsnap_filter_1_fa
    File? host_filter_out_gsnap_filter_merged_fa
    File lineage_db = "s3://czid-public-references/taxonomy/2021-01-22/taxid-lineages.db"
    File accession2taxid_db = "s3://czid-public-references/alignment_data/2021-01-22/accession2taxid.db"
    File taxon_blacklist = "s3://czid-public-references/taxonomy/2021-01-22/taxon_blacklist.txt"
    Int min_read_length = 36
    File deuterostome_db = "s3://czid-public-references/taxonomy/2021-01-22/deuterostome_taxids.txt"
    File? minimap2_local_db_path
    String minimap2_db = "s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nt_20_long/"
    String minimap2_args = "-x asm20 --secondary=yes"
    String minimap2_prefix = "gsnap"
  }

  # HACK: RunAlignment_minimap2_out depends on the fastas being in s3_wd_uri with their basename
  #   this step outputs the input fastas so swipe will upload them to s3_wd_uri so RunAlignment_minimap2_out
  #   can run
  call DownloadFastas {
    input:
      docker_image_id = docker_image_id,
      fastas = [select_first([host_filter_out_gsnap_filter_merged_fa, host_filter_out_gsnap_filter_1_fa])], #select_all([host_filter_out_gsnap_filter_1_fa, host_filter_out_gsnap_filter_2_fa]),
  }

  call RunAlignment_minimap2_out { 
    input:         
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      fastas = DownloadFastas.loaded_fastas,
      db_path = minimap2_db,
      minimap2_args = minimap2_args,
      run_locally = defined(minimap2_local_db_path),
      local_minimap2_index = minimap2_local_db_path,
      prefix= minimap2_prefix
  }

  call RunCallHitsMinimap2{ 
    input:
      m8_file = RunAlignment_minimap2_out.out_m8,
      lineage_db = lineage_db,
      taxon_blacklist = taxon_blacklist,
      deuterostome_db = deuterostome_db,
      accession2taxid = accession2taxid_db,
      min_read_length = min_read_length,
      docker_image_id = docker_image_id,
  }

  output {
    File gsnap_out_gsnap_m8 = RunAlignment_minimap2_out.out_m8
    File gsnap_out_gsnap_deduped_m8 = RunCallHitsMinimap2.deduped_out_m8
    File gsnap_out_gsnap_hitsummary_tab = RunCallHitsMinimap2.hitsummary
    File gsnap_out_gsnap_counts_with_dcr_json = RunCallHitsMinimap2.counts_json
    File? minimap2_version = RunAlignment_minimap2_out.version
    File? gsnap_out_count = RunCallHitsMinimap2.output_read_count
  }
}
