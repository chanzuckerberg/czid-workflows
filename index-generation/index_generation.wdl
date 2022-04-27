version development

workflow index_generation {
    input {
        String index_name
        String ncbi_server = "https://ftp.ncbi.nih.gov"
        File? previous_lineages = ""
        String docker_image_id
    }

    call DownloadIndexSources {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id
    }

    call GenerateIndexAccessions {
        input:
        nr = DownloadIndexSources.nr,
        nt = DownloadIndexSources.nt,
        accession2taxid = DownloadIndexSources.accession2taxid,
        docker_image_id = docker_image_id
    }

    call GenerateIndexDiamond {
        input:
        nr = DownloadIndexSources.nr,
        docker_image_id = docker_image_id
    }

    call GenerateIndexLineages {
        input:
        taxdump = DownloadIndexSources.taxdump,
        index_name = index_name,
        docker_image_id = docker_image_id
    }

    call GenerateIndexMinimap2 {
        input:
        nt = DownloadIndexSources.nt,
        docker_image_id = docker_image_id
    }

    output {
        File nr = DownloadIndexSources.nr
        File nt = DownloadIndexSources.nt
        Directory accession2taxid = DownloadIndexSources.accession2taxid
        File taxdump = DownloadIndexSources.taxdump
        File accession2taxid_gz = GenerateIndexAccessions.accession2taxid_gz
        File accession2taxid_wgs = GenerateIndexAccessions.accession2taxid_wgs
        File accession2taxid_db = GenerateIndexAccessions.accession2taxid_db
        File taxid2wgs_accession_db = GenerateIndexAccessions.taxid2wgs_accession_db
        File nt_loc_db = GenerateIndexAccessions.nt_loc_db
        File nt_info_db = GenerateIndexAccessions.nt_info_db
        File nr_loc_db = GenerateIndexAccessions.nr_loc_db
        File nr_info_db = GenerateIndexAccessions.nr_info_db
        Directory diamond_index = GenerateIndexDiamond.diamond_index
        File taxid_lineages_db = GenerateIndexLineages.taxid_lineages_db
        File taxid_lineages_csv = GenerateIndexLineages.taxid_lineages_csv
        File names_csv = GenerateIndexLineages.names_csv
        File named_taxid_lineages_csv = GenerateIndexLineages.named_taxid_lineages_csv
        File versioned_taxid_lineages_csv = GenerateIndexLineages.versioned_taxid_lineages_csv
        File phage_list_csv = GenerateIndexLineages.phage_list_csv
        File deuterostome_taxids = GenerateIndexLineages.deuterostome_taxids
        File taxon_ignore_list = GenerateIndexLineages.taxon_ignore_list
        Directory minimap2_index = GenerateIndexMinimap2.minimap2_index
    }
}

task DownloadIndexSources {
    input {
        String ncbi_server
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        download () {
            mkdir -p $(dirname $1)
            wget -P $(dirname $1) -cnv ~{ncbi_server}/$1
            wget -P $(dirname $1) -cnv ~{ncbi_server}/$1.md5

            if [[ $(md5sum $1 | cut -f 1 -d' ') != $(cat $1.md5 | sed -e 's/^\s*//' | cut -f 1 -d' ') ]]
            then
                exit 1
            fi
        }
        
        download blast/db/FASTA/nt.gz
        download blast/db/FASTA/nr.gz
        download pub/taxonomy/taxdump.tar.gz
        download pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz
        download pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz
        download pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz
        download pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        download pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
        download pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
        download pub/taxonomy/accession2taxid/prot.accession2taxid.gz
        download pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

        pigz -dc blast/db/FASTA/nt.gz > nt
        pigz -dc blast/db/FASTA/nr.gz > nr
    >>>

    output {
        File nr = "nr"
        File nt = "nt"
        Directory accession2taxid = "pub/taxonomy/accession2taxid"
        File taxdump = "pub/taxonomy/taxdump.tar.gz"
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

        cp -r ~{accession2taxid} accession2taxid

        # Build index
        python3 /usr/local/bin/generate_accession2taxid.py \
            accession2taxid/nucl_wgs.accession2taxid.gz \
            accession2taxid/nucl_gb.accession2taxid.gz \
            accession2taxid/pdb.accession2taxid.gz \
            accession2taxid/prot.accession2taxid.FULL.gz \
            --nt_file ~{nt} \
            --nr_file ~{nr} \
            --output_gz accession2taxid.gz \
            --output_wgs_gz accession2taxid_wgs.gz \
            --accession2taxid_db accession2taxid.db \
            --taxid2wgs_accession_db taxid2wgs_accession.db
        python3 /usr/local/bin/generate_loc_db.py ~{nt} nt_loc.db nt_info.db
        python3 /usr/local/bin/generate_loc_db.py ~{nr} nr_loc.db nr_info.db
    >>>

    output {
        File accession2taxid_gz = "accession2taxid.gz"
        File accession2taxid_wgs = "accession2taxid_wgs.gz"
        File accession2taxid_db = "accession2taxid.db"
        File taxid2wgs_accession_db = "taxid2wgs_accession.db"
        File nt_loc_db = "nt_loc.db"
        File nt_info_db = "nt_info.db"
        File nr_loc_db = "nr_loc.db"
        File nr_info_db = "nr_info.db"
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
        diamond makedb --in ~{nr} -d diamond_index_chunksize_~{chunksize} --scatter-gather -b ~{chunksize}
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
        File? previous_lineages = ""
        String index_name
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        # Build Indexes
        git clone https://github.com/chanzuckerberg/ncbitax2lin.git
        cd ncbitax2lin
        mkdir -p taxdump/taxdump
        tar zxf ~{taxdump} -C ./taxdump/taxdump
        make 1>&2

        # Add names to lineages

        python3 /usr/local/bin/generate_lineage_csvs.py \
            names.csv.gz \
            taxid-lineages.csv.gz \
            ~{index_name} \
            named-taxid-lineages.csv.gz \
            versioned-taxid-lineages.csv.gz \
            ~{previous_lineages}

        python3 /usr/local/bin/generate_phage_list.py versioned-taxid-lineages.csv.gz versioned_phage_list.csv.gz

        # Build deuterostome list
        # decompress first and only read what we need to prevent pipefail
        # lineages.csv.gz generated by ncbitax2lin
        gzip -d lineages.csv.gz
        TAXID_COL_NUM=$(head -n1 "lineages.csv" | tr "," "\n" | grep -n 'tax_id' | cut -f1 -d":")
        if [[  "$TAXID_COL_NUM" == ""  ]]; then
          >&2 echo "ERROR - Couldn't find tax_id header column in lineages.csv"
          exit 1
        fi

        cat lineages.csv | grep 'Chordata\|Echinodermata\|Hemichordata' | cut -f"$TAXID_COL_NUM" -d"," > deuterostome_taxids.txt
        # TODO: refine which taxa we want to ignore
        cat lineages.csv | grep 'vector\|plasmid\|Plasposon\|replicon\|synthetic\|construct\|Artificial\|Recombinant\|insert\|cassette' | cut -f"$TAXID_COL_NUM" -d"," > taxon_ignore_list.txt
    >>>

    output {
        File taxid_lineages_db = "ncbitax2lin/taxid-lineages.db"
        File taxid_lineages_csv = "ncbitax2lin/taxid-lineages.csv.gz"
        File names_csv = "ncbitax2lin/names.csv.gz"
        File named_taxid_lineages_csv = "ncbitax2lin/named-taxid-lineages.csv.gz"
        File versioned_taxid_lineages_csv = "ncbitax2lin/versioned-taxid-lineages.csv.gz"
        File phage_list_csv = "ncbitax2lin/versioned_phage_list.csv.gz"
        File deuterostome_taxids = "ncbitax2lin/deuterostome_taxids.txt"
        File taxon_ignore_list = "ncbitax2lin/taxon_ignore_list.txt"
    }

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
                minimap2 -cx sr -k ~{k} -w ~{w} -I ~{I} -t ~{t} -d $OUTDIR/"nt.part_"$path".idx" $i
        done
    >>>

    output {
        Directory minimap2_index = "nt_k~{k}_w~{w}_~{n_chunks}"
    }

    runtime {
        docker: docker_image_id
    }
}
