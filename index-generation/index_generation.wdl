version 1.1

workflow index_generation {
    input {
        String index_name
        String docker_image_id
    }

    call DownloadIndexSources {
        input:
        index_name = index_name
    }

    call GenerateIndexAccessions {
        input:
        index_name = index_name,
        dummy = DownloadIndexSources.dummy
    }

    call GenerateIndexDiamond {
        input:
        index_name = index_name,
        dummy = DownloadIndexSources.dummy
    }

    call GenerateIndexLineages {
        input:
        index_name = index_name,
        dummy = DownloadIndexSources.dummy
    }

    call GenerateIndexMinimap2 {
        input:
        index_name = index_name,
        dummy = DownloadIndexSources.dummy
    }

    output {}
}

task DownloadIndexSources {
    input {
        String index_name
    }

    command <<<
        python3 /usr/local/bin/download_index_sources.py
    >>>

    output {
        String dummy = "dummy"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexAccessions {
    input {
        String index_name
        String dummy
    }

    command <<<
        #!/bin/bash -ex

        set -o pipefail

        # Get $INDEX_NAME
        s3parcp "s3://$BUCKET/ncbi-sources/index_name"
        INDEX_NAME=$(cat index_name)

        # Download nr
        s3parcp --checksum "s3://$BUCKET/ncbi-sources/$INDEX_NAME/nr"

        # Download nt
        s3parcp --checksum "s3://$BUCKET/ncbi-sources/$INDEX_NAME/nt"

        # Download accession2taxid
        s3parcp --checksum --recursive "s3://$BUCKET/ncbi-sources/$INDEX_NAME/accession2taxid/" accession2taxid/

        # Build index
        python3 /usr/local/bin/generate_accession2taxid.py \
            accession2taxid/nucl_wgs.accession2taxid.gz \
            accession2taxid/nucl_gb.accession2taxid.gz \
            accession2taxid/pdb.accession2taxid.gz \
            accession2taxid/prot.accession2taxid.FULL.gz \
            --nt_file nt \
            --nr_file nr \
            --output_gz accession2taxid.gz \
            --output_wgs_gz accession2taxid_wgs.gz \
            --accession2taxid_db accession2taxid.db \
            --taxid2wgs_accession_db taxid2wgs_accession.db
        python3 /usr/local/bin/generate_loc_db.py nt nt_loc.db nt_info.db
        python3 /usr/local/bin/generate_loc_db.py nr nr_loc.db nr_info.db

        # Output dir in S3
        S3_OUTPUT_DIR="s3://$BUCKET/alignment_data/$INDEX_NAME/"

        # Upload index
        s3parcp --checksum accession2taxid.gz "$S3_OUTPUT_DIR"
        s3parcp --checksum accession2taxid_wgs.gz "$S3_OUTPUT_DIR"
        s3parcp --checksum accession2taxid.db "$S3_OUTPUT_DIR"
        s3parcp --checksum taxid2wgs_accession.db "$S3_OUTPUT_DIR"
        s3parcp --checksum nt_loc.db "$S3_OUTPUT_DIR"
        s3parcp --checksum nt_info.db "$S3_OUTPUT_DIR"
        s3parcp --checksum nr_loc.db "$S3_OUTPUT_DIR"
        s3parcp --checksum nr_info.db "$S3_OUTPUT_DIR"

        # Generate and store stats report
        python3 /usr/local/bin/compare_index_files.py "$S3_OUTPUT_DIR" > alignment_data_files_stats.txt
        s3parcp --checksum alignment_data_files_stats.txt "$S3_OUTPUT_DIR"
    >>>

    output {}

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexDiamond {
    input {
        String index_name
        String dummy
    }

    command <<<
        #!/bin/bash

        set -eo pipefail
        CHUNK_SIZE=5500000000 # Not sure if there's a better way of doing this. 

        # Get $INDEX_NAME
        s3parcp "s3://$BUCKET/ncbi-sources/index_name"
        INDEX_NAME=$(cat index_name)


        # Download nt
        s3parcp --checksum "s3://$BUCKET/ncbi-sources/$INDEX_NAME/nr
        OUTDIR=ref
        #Run diamond
        diamond/diamond makedb --in nr -d $OUTDIR --scatter-gather -b $CHUNK_SIZE
        # Output directory in S3. NOTE: no trailing slash.
        S3_OUTPUT_DIR="s3://$BUCKET/alignment_indexes/$INDEX_NAME"
        # Upload index
        s3parcp --checksum --recursive $OUTDIR "$S3_OUTPUT_DIR/$OUTDIR"
    >>>

    output {}

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexLineages {
    input {
        String index_name
        String dummy
    }

    command <<<
        #!/bin/bash -ex

        set -o pipefail

        # Get $INDEX_NAME
        s3parcp "s3://$BUCKET/ncbi-sources/index_name"
        INDEX_NAME=$(cat index_name)

        # Download taxdump
        s3parcp --checksum "s3://$BUCKET/ncbi-sources/$INDEX_NAME/taxdump.tar.gz"

        # Build Indexes
        git clone https://github.com/chanzuckerberg/ncbitax2lin.git
        cd ncbitax2lin
        mkdir -p taxdump/taxdump
        tar zxf ../taxdump.tar.gz -C ./taxdump/taxdump
        make

        # Output directory in S3
        S3_OUTPUT_DIR="s3://$BUCKET/taxonomy/$INDEX_NAME/"

        # Upload indexes
        s3parcp --checksum taxid-lineages.db "$S3_OUTPUT_DIR"
        s3parcp --checksum taxid-lineages.csv.gz "$S3_OUTPUT_DIR"
        s3parcp --checksum names.csv.gz "$S3_OUTPUT_DIR"

        # Build deuterostome list
        # decompress first and only read what we need to prevent pipefail
        gzip -d lineages.csv.gz
        TAXID_COL_NUM=$(head -n1 "lineages.csv" | tr "," "\n" | grep -n 'tax_id' | cut -f1 -d":")
        if [[  "$TAXID_COL_NUM" == ""  ]]; then
          >&2 echo "ERROR - Couldn't find tax_id header column in lineages.csv"
          exit 1
        fi

        cat lineages.csv |  grep 'Chordata\|Echinodermata\|Hemichordata' | cut -f"$TAXID_COL_NUM" -d"," > deuterostome_taxids.txt

        # These keywords are based on reverse engineering the original taxon_blacklist at
        # s3://czid-public-references/taxonomy/2018-04-01-utc-1522569777-unixtime__2018-04-04-utc-1522862260-unixtime/taxon_blacklist.txt
        # KEYWORDS='\bvector\b|\bplasmid\b|\bPlasposon\b|\breplicon\b|\bsynthetic\b|\bconstruct\b|\bArtificial\b|\bRecombinant\b|\binsert\b|\bcassette\b'
        # egrep "$KEYWORDS" lineages.csv | cut -f"$TAXID_COL_NUM" -d"," >
        # taxon_blacklist.txt
        # TODO: (gdingle): Decide on what to include in updated taxon_blacklist. Until
        # then, we will simply continue to use the original. See
        # https://github.com/chanzuckerberg/idseq/pull/131/ and
        # https://jira.czi.team/browse/IDSEQ-2760 .
        BLACKLIST_S3_PATH=${BLACKLIST_S3_PATH:-"s3://czid-public-references/taxonomy/2018-04-01-utc-1522569777-unixtime__2018-04-04-utc-1522862260-unixtime/taxon_blacklist.txt"}
        s3parcp $BLACKLIST_S3_PATH

        # Upload deuterostome and taxon_blacklist.txt lists
        s3parcp --checksum deuterostome_taxids.txt "$S3_OUTPUT_DIR"
        s3parcp --checksum taxon_blacklist.txt "$S3_OUTPUT_DIR"

        # Generate and store stats report
        python3 /usr/local/bin/compare_index_files.py "$S3_OUTPUT_DIR" > taxonomy_files_stats.txt
        s3parcp --checksum taxonomy_files_stats.txt "$S3_OUTPUT_DIR"
    >>>

    output {}

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexMinimap2 {
    input {
        String index_name
        String dummy
    }

    command <<<
        #!/bin/bash

        set -eo pipefail

        k=12 # Minimizer k-mer length default is 21 for short reads option
        w=8 # Minimizer window size default is 11 for short reads option

        I=9999G # Load at most NUM target bases into RAM for indexing
        t=20 # number of threads, doesn't really work for indexing I don't think
        CHUNKS=20

        # Get $INDEX_NAME
        s3parcp "s3://$BUCKET/ncbi-sources/index_name"
        INDEX_NAME=$(cat index_name)

        # Download nt
        s3parcp --checksum "s3://$BUCKET/ncbi-sources/$INDEX_NAME/nt"

        # Split nt into 20
        seqkit split2 nt -p $CHUNKS

        # Make output directory
        OUTDIR="nt_k"$k"_w"$w"_"$CHUNKS
        mkdir $OUTDIR

        # Run minimap2 on each chunk
        for i in nt.split/nt*
        do
                path="${i##*_}"
                minimap2 -cx sr -k $k -w $w -I $I -t $t -d $OUTDIR/"genome_"$path".mmi" $i
        done


        # Output directory in S3. NOTE: no trailing slash.
        S3_OUTPUT_DIR="s3://$BUCKET/alignment_indexes/$INDEX_NAME"

        # Upload index
        s3parcp --checksum --recursive $OUTDIR "$S3_OUTPUT_DIR/$OUTDIR"
    >>>

    output {}

    runtime {
        docker: docker_image_id
    }
}
