version development

workflow index_generation {
    input {
        # String index_name
        String ncbi_server = "https://ftp.ncbi.nih.gov"
        # Boolean write_to_db = false
        # String? environ
        # TODO: (alignment_config) remove after alignment config table is removed
        # String? s3_dir
        # File? previous_lineages
        Array[String] taxids_to_drop = ["9606"] # human
        Int nt_compression_k = 31
        Int nt_compression_scaled = 100000
        Float nt_compression_similarity_threshold = 0.9
        String docker_image_id
    }

    call DownloadNR {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id
    }

    call DownloadNT {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id
    }

    call DownloadAccession2Taxid {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id
    }

    # call DownloadTaxdump {
    #     input:
    #     ncbi_server = ncbi_server,
    #     docker_image_id = docker_image_id
    # }

    call GenerateIndexAccessions {
        input:
        nr = DownloadNR.nr,
        nt = DownloadNT.nt,
        accession2taxid = DownloadAccession2Taxid.accession2taxid,
        docker_image_id = docker_image_id
    }

    call CompressNT {
        input:
        nt = DownloadNT.nt,
        accession2taxid = GenerateIndexAccessions.accession2taxid_db,
        taxids_to_drop = taxids_to_drop,
        k = nt_compression_k,
        scaled = nt_compression_scaled,
        similarity_threshold = nt_compression_similarity_threshold,
        docker_image_id = docker_image_id,
        cpu = 32
    }

    call GenerateNTDB {
        input:
        nt = CompressNT.nt_compressed,
        docker_image_id = docker_image_id
    }

    # call GenerateNRDB {
    #     input:
    #     nr = DownloadNR.nr,
    #     docker_image_id = docker_image_id
    # }

    # call GenerateIndexDiamond {
    #     input:
    #     nr = DownloadNR.nr,
    #     docker_image_id = docker_image_id
    # }

    # call GenerateIndexLineages {
    #     input:
    #     taxdump = DownloadTaxdump.taxdump,
    #     index_name = index_name,
    #     previous_lineages = previous_lineages,
    #     docker_image_id = docker_image_id
    # }
    
    call GenerateIndexMinimap2 {
        input:
        nt = CompressNT.nt_compressed,
        docker_image_id = docker_image_id
    }


    # if (write_to_db && defined(environ) && defined(s3_dir)) {
    #     call LoadTaxonLineages {
    #         input:
    #         environ = environ,
    #         index_name = index_name,
    #         s3_dir = s3_dir,
    #         versioned_taxid_lineages_csv = GenerateIndexLineages.versioned_taxid_lineages_csv,
    #         docker_image_id = docker_image_id
    #     } 
    # }


    output {
        File nr = DownloadNR.nr
        File nt = CompressNT.nt_compressed
        File accession2taxid_db = GenerateIndexAccessions.accession2taxid_db
        File nt_loc_db = GenerateNTDB.nt_loc_db
        File nt_info_db = GenerateNTDB.nt_info_db
        # File nr_loc_db = GenerateNRDB.nr_loc_db
        # Directory diamond_index = GenerateIndexDiamond.diamond_index
        # File taxid_lineages_db = GenerateIndexLineages.taxid_lineages_db
        # File versioned_taxid_lineages_csv = GenerateIndexLineages.versioned_taxid_lineages_csv
        # File deuterostome_taxids = GenerateIndexLineages.deuterostome_taxids
        # File taxon_ignore_list = GenerateIndexLineages.taxon_ignore_list
        Directory minimap2_index = GenerateIndexMinimap2.minimap2_index
    }
}

task DownloadNR {
    input {
        String ncbi_server
        String docker_image_id
    }

    command <<<
        ncbi_download "~{ncbi_server}" blast/db/FASTA/nr.gz
    >>>

    output {
        File nr = "blast/db/FASTA/nr"
    }

    runtime {
        docker: docker_image_id
    }
}

task DownloadNT {
    input {
        String ncbi_server
        String docker_image_id
    }

    command <<<
        ncbi_download "~{ncbi_server}" blast/db/FASTA/nt.gz
        
    >>>

    output {
        File nt = "blast/db/FASTA/nt"
    }

    runtime {
        docker: docker_image_id
    }
}

task DownloadAccession2Taxid {
    input {
        String ncbi_server
        String docker_image_id
    }

    command <<<
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
    >>>

    output {
        Directory accession2taxid = "pub/taxonomy/accession2taxid"
    }

    runtime {
        docker: docker_image_id
    }
}

task DownloadTaxdump {
    input {
        String ncbi_server
        String docker_image_id
    }

    command <<<
        ncbi_download "~{ncbi_server}" pub/taxonomy/taxdump.tar.gz
    >>>

    output {
        File taxdump = "pub/taxonomy/taxdump.tar"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexAccessions {
    input {
        File nr
        File nt
        Directory accession2taxid
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        # Build index
        python3 /usr/local/bin/generate_accession2taxid.py \
            ~{accession2taxid}/nucl_wgs.accession2taxid \
            ~{accession2taxid}/nucl_gb.accession2taxid \
            ~{accession2taxid}/pdb.accession2taxid \
            ~{accession2taxid}/prot.accession2taxid.FULL \
            --nt_file ~{nt} \
            --nr_file ~{nr} \
            --accession2taxid_db accession2taxid.marisa \
    >>>

    output {
        File accession2taxid_db = "accession2taxid.marisa"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateNTDB {
    input {
        File nt
        String docker_image_id
    }

    command <<<
        python3 /usr/local/bin/generate_loc_db.py ~{nt} nt_loc.marisa nt_info.marisa
    >>>

    output {
        File nt_loc_db = "nt_loc.marisa"
        File nt_info_db = "nt_info.marisa"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateNRDB {
    input {
        File nr
        String docker_image_id
    }

    command <<<
        python3 /usr/local/bin/generate_loc_db.py ~{nr} nr_loc.marisa nr_info.marisa
    >>>

    output {
        File nr_loc_db = "nr_loc.marisa"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexDiamond {
    input {
        File nr
        Int chunksize = 5500000000
        String docker_image_id
    }

    command <<<
        # Ignore warning is needed because sometimes NR has sequences of only DNA characters which causes this to fail
        diamond makedb --ignore-warnings --in ~{nr} -d diamond_index_chunksize_~{chunksize} --scatter-gather -b ~{chunksize}
    >>>

    output {
        Directory diamond_index = "diamond_index_chunksize_~{chunksize}"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexLineages {
    input {
        File taxdump
        File? previous_lineages
        String index_name
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        # Build Indexes
        mkdir -p taxdump/taxdump
        tar xf ~{taxdump} -C ./taxdump/taxdump

        python3 /usr/local/bin/ncbitax2lin.py \
            --nodes-file taxdump/taxdump/nodes.dmp \
            --names-file taxdump/taxdump/names.dmp \
            --names-output-prefix names \
            --taxid-lineages-output-prefix taxid-lineages \
            --name-lineages-output-prefix lineages

        # Add names to lineages

        python3 /usr/local/bin/generate_lineage_csvs.py \
            names.csv.gz \
            taxid-lineages.csv.gz \
            ~{index_name} \
            named-taxid-lineages.csv.gz \
            versioned-taxid-lineages.csv.gz \
            ~{previous_lineages}

        # Build deuterostome list
        # decompress first and only read what we need to prevent pipefail
        # lineages.csv.gz generated by ncbitax2lin
        gzip -d lineages.csv.gz
        TAXID_COL_NUM=$(head -n1 "lineages.csv" | tr "," "\n" | grep -n 'tax_id' | cut -f1 -d":")
        if [[  "$TAXID_COL_NUM" == ""  ]]; then
          >&2 echo "ERROR - Couldn't find tax_id header column in lineages.csv"
          exit 1
        fi

        grep 'Chordata\|Echinodermata\|Hemichordata' lineages.csv | cut -f"$TAXID_COL_NUM" -d"," > deuterostome_taxids.txt
        # TODO: refine which taxa we want to ignore
        grep 'vector\|plasmid\|Plasposon\|replicon\|synthetic\|construct\|Artificial\|Recombinant\|insert\|cassette' lineages.csv | cut -f"$TAXID_COL_NUM" -d"," > taxon_ignore_list.txt
    >>>

    output {
        File taxid_lineages_db = "taxid-lineages.marisa"
        File versioned_taxid_lineages_csv = "versioned-taxid-lineages.csv.gz"
        File deuterostome_taxids = "deuterostome_taxids.txt"
        File taxon_ignore_list = "taxon_ignore_list.txt"
    }

    runtime {
        docker: docker_image_id
    }
}

task LoadTaxonLineages {
    input {
        String? environ
        String index_name
        String? s3_dir
        File versioned_taxid_lineages_csv
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        get_param () {
            aws ssm get-parameter --name "$1" --with-decryption | jq -r '.Parameter.Value'
        }

        HOST=$(get_param "/idseq-~{environ}-web/RDS_ADDRESS")
        USER=$(get_param "/idseq-~{environ}-web/DB_USERNAME")
        DATABASE="idseq_~{environ}"

        {
            echo "[client]"
            echo "protocol=tcp"
            echo "host=$HOST"
            echo "user=$USER"
        } > my.cnf

        # Add the password without making it a part of the command so we can print commands via set -x without exposing the password
        # Add the password= for the config
        echo "password=" >> my.cnf
        # Remove the newline at end of file
        truncate -s -1 my.cnf
        # Append the password after the =
        get_param "/idseq-~{environ}-web/db_password" >> my.cnf

        gzip -dc ~{versioned_taxid_lineages_csv} > "taxon_lineages_new.csv"
        COLS=$(head -n 1 "taxon_lineages_new.csv")
        mysql --defaults-extra-file=my.cnf -D "$DATABASE" -e "CREATE TABLE taxon_lineages_new LIKE taxon_lineages"
        mysqlimport --defaults-extra-file=my.cnf --verbose --local --columns="$COLS" --fields-terminated-by=',' --fields-optionally-enclosed-by='"' --ignore-lines 1 "$DATABASE" "taxon_lineages_new.csv"
        mysql --defaults-extra-file=my.cnf -D "$DATABASE" -e "RENAME TABLE taxon_lineages TO taxon_lineages_old, taxon_lineages_new To taxon_lineages"
        # TODO: remove old table once we feel safe
        # mysql --defaults-extra-file=my.cnf -D "$DATABASE" -e "DROP TABLE taxon_lineages_old"
        # TODO: (alignment_config) remove after alignment config table is removed
        mysql --defaults-extra-file=my.cnf -D "$DATABASE" -e "
            INSERT INTO alignment_configs(
                name,
                s3_nt_db_path,
                s3_nt_loc_db_path,
                s3_nr_db_path,
                s3_nr_loc_db_path,
                s3_lineage_path,
                s3_accession2taxid_path,
                s3_deuterostome_db_path,
                s3_nt_info_db_path,
                s3_taxon_blacklist_path,
                lineage_version,
                created_at,
                updated_at
            ) VALUES(
                '~{index_name}',
                '~{s3_dir}/nt',
                '~{s3_dir}/nt_loc.marisa',
                '~{s3_dir}/nr',
                '~{s3_dir}/nr_loc.marisa',
                '~{s3_dir}/taxid-lineages.marisa',
                '~{s3_dir}/accession2taxid.marisa',
                '~{s3_dir}/deuterostome_taxids.txt',
                '~{s3_dir}/nt_info.marisa',
                '~{s3_dir}/taxon_ignore_list.txt',
                '~{index_name}',
                NOW(),
                NOW()
            ); 
        "
    >>>

    output {}

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexMinimap2 {
    input {
        File nt
        Int k = 14 # Minimizer k-mer length default is 21 for short reads option
        Int w = 8 # Minimizer window size default is 11 for short reads option
        String I = "9999G" # Load at most NUM target bases into RAM for indexing
        Int t = 20 # number of threads, doesn't really work for indexing I don't think
        Int n_chunks = 20
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        # Split nt into 20
        seqkit split2 ~{nt} -p ~{n_chunks} --out-dir nt.split

        # Make output directory
        OUTDIR="nt_k~{k}_w~{w}_~{n_chunks}"
        mkdir $OUTDIR

        # Run minimap2 on each chunk
        for i in nt.split/*
        do
                path="${i##*_}"
                minimap2 -cx sr -k ~{k} -w ~{w} -I ~{I} -t ~{t} -d "${OUTDIR}/nt.part_${path}.idx" "$i"
        done
    >>>

    output {
        Directory minimap2_index = "nt_k~{k}_w~{w}_~{n_chunks}"
    }

    runtime {
        docker: docker_image_id
    }
}

task CompressNT {
    input {
        File nt
        File accession2taxid
        Array[String] taxids_to_drop
        Int k
        Int scaled
        Float similarity_threshold
        String docker_image_id
        Int cpu
    }

    command <<< 
        set -euxo pipefail

        # Sort NT by length with the longer sequences first
        #   This is needed because the compression algorithm iterates through NT in order only emitting
        #   sequences if they are not contained by what it has already seen. If a shorter sequence is
        #   contained by a longer sequence, and the shorter sequence were to come first, it would be emitted
        #   even though it is redundant to the longer sequence.
        seqkit sort --reverse --by-length --two-pass --threads ~{cpu} ~{nt} -o nt_sorted

        # Rust can't read marisa tries so we need to convert to csv
        python3  /usr/local/bin/accession2taxid_csv.py \
            --accession2taxid-path ~{accession2taxid} \
            ~{ if length(taxids_to_drop) > 0 then "--taxids-to-drop ~{sep(" ", taxids_to_drop)}" else "" } \
            --csv-filepath accession2taxid.csv

        ncbi-compress \
            --input-fasta nt_sorted \
            --accession-to-taxid-csv accession2taxid.csv \
            --output-fasta nt_compressed.fa \
            --k ~{k} \
            --scaled ~{scaled} \
            --similarity-threshold ~{similarity_threshold} \
    >>>

    output {
        File nt_compressed = "nt_compressed.fa"
    }

    runtime {
        docker: docker_image_id
        cpu: cpu
    }
}
