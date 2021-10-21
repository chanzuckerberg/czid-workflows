version 1.1

workflow alignment_scalability {
    input {
        File fastqs_0
        File? fastqs_1
        String input_dir
        String chunk_dir
        String db_path
        String docker_image_id
        File duplicate_cluster_size
        File lineage_db="s3://idseq-public-references/taxonomy/2020-04-20/taxid-lineages.db"
        File taxon_blacklist="s3://idseq-public-references/taxonomy/2020-04-20/taxon_blacklist.txt"
        File deuterostome_db="s3://idseq-public-references/taxonomy/2020-04-20/deuterostome_taxids.txt"
        File accession2taxid="s3://idseq-public-references/alignment_data/2020-04-20/accession2taxid.db"
        String prefix = "output"
        Int min_read_length = 36
        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call RunMinimap {
        input:
        fastqs = select_all([fastqs_0, fastqs_1]),
        input_dir = input_dir,
        chunk_dir = chunk_dir,
        db_path = db_path,
        docker_image_id = docker_image_id
    }
    call RunCallHits { 
        input:
        m8_file = RunMinimap.out_m8,
        lineage_db = lineage_db,
        duplicate_cluster_size = duplicate_cluster_size,
        taxon_blacklist = taxon_blacklist,
        deuterostome_db = deuterostome_db,
        accession2taxid = accession2taxid,
        prefix = prefix,
        min_read_length = min_read_length,
        docker_image_id = docker_image_id
    }

    output {
        File out_paf = RunMinimap.out_paf
        File out_m8 = RunMinimap.out_m8 
        File deduped_out_m8 = RunCallHits.deduped_out_m8 
        File hitsummary = RunCallHits.hitsummary
    }
}

task RunMinimap {
    input {
        Array[File]+ fastqs
        String input_dir
        String chunk_dir
        String db_path
        String docker_image_id
        String prefix = "output"
    }

    command <<<
        set -e 
        echo STARTING
        export DEPLOYMENT_ENVIRONMENT=dev
        export AWS_REGION="us-west-2"
        export AWS_DEFAULT_REGION="us-west-2"
        python3 <<CODE
        from idseq_utils.run_minimap2 import run_minimap2
        fastqs = ["~{sep='", "' fastqs}"]
        run_minimap2("~{input_dir}", "~{chunk_dir}", "~{db_path}", "~{prefix}.paf", *fastqs)
        CODE
        python3 /usr/local/lib/python3.6/dist-packages/idseq_utils/paf2blast6.py "~{prefix}".paf
    >>>

    output {
        File out_paf = "~{prefix}.paf"
        File out_m8 = "~{prefix}_frompaf.m8"
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
        set -e 
        python3 <<CODE
        from idseq_dag.util.m8 import call_hits_m8, generate_taxon_count_json_from_m8
        call_hits_m8("~{m8_file}", "~{lineage_db}", "~{accession2taxid}", "~{prefix}_frompaf_deduped.m8", "~{prefix}.hitsummary.tab", ~{min_read_length}, "~{deuterostome_db}", None, "~{taxon_blacklist}")
        generate_taxon_count_json_from_m8( "~{prefix}_frompaf_deduped.m8", "~{prefix}.hitsummary.tab", "NT", "~{lineage_db}", "~{deuterostome_db}", None, "~{taxon_blacklist}", "~{duplicate_cluster_size}", "~{prefix}.json")
        CODE
        >>>

    output {
        File deduped_out_m8 = "~{prefix}_frompaf_deduped.m8"
        File hitsummary = "~{prefix}.hitsummary.tab"
    }

    runtime {
        docker: docker_image_id
    }
}
