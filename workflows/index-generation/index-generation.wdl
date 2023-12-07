version development

workflow index_generation {
    input {
        String index_name
        String ncbi_server = "https://ftp.ncbi.nih.gov"
        # Boolean write_to_db = false
        String? environ
        # TODO: (alignment_config) remove after alignment config table is removed
        String? s3_dir
        File? previous_lineages
        Int nt_compression_k = 31
        Int nt_compression_scaled = 1000
        Float nt_compression_similarity_threshold = 0.5
        Int nr_compression_k = 31
        Int nr_compression_scaled = 100
        Float nr_compression_similarity_threshold = 0.5
        String s3_accession_mapping_prefix = "" # old accession mapping files
        String docker_image_id
        String old_nt_s3_path = "" # old nt file
        String old_nr_s3_path = "" # old nr file
        Boolean skip_protein_compression = false
        Boolean skip_nuc_compression = false
        Boolean logging_enabled = false
    }

    call DownloadNR {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id,
        old_nr_s3_path = old_nr_s3_path
    }
    call DownloadNT {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id,
        old_nt_s3_path = old_nt_s3_path
    }
    call DownloadAccession2Taxid {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id,
        s3_accession_mapping_prefix = s3_accession_mapping_prefix

    }
    call DownloadTaxdump {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id
    }
    if (!skip_protein_compression) {
        call SplitFastaBySeqLengthAndSort as SplitFastaBySeqLengthAndSortNR {
            input:
            fasta = DownloadNR.nr,
            combined_sorted_path = "combined_sorted_nr.fa",
            cpu = 64,
            docker_image_id = docker_image_id
        }
        call BreakApartByTaxid as BreakApartByTaxidNR {
            input:
            sorted_fasta = SplitFastaBySeqLengthAndSortNR.sorted,
            accession2taxid = DownloadAccession2Taxid.accession2taxid,
            reads_by_taxid_path = "reads_by_taxid_nr",
            cpu = 64,
            docker_image_id = docker_image_id
        }

        call CompressNR {
            input:
            reads_by_taxid = BreakApartByTaxidNR.reads_by_taxid,
            k = nr_compression_k,
            scaled = nr_compression_scaled,
            similarity_threshold = nr_compression_similarity_threshold,
            logging_enabled = logging_enabled,
            docker_image_id = docker_image_id,
        }
    }

    if (!skip_nuc_compression) {
        call SplitFastaBySeqLengthAndSort as SplitFastaBySeqLengthAndSortNT {
            input:
            fasta = DownloadNT.nt,
            combined_sorted_path = "combined_sorted_nt.fa",
            cpu = 64,
            docker_image_id = docker_image_id
        }
        call BreakApartByTaxid as BreakApartByTaxidNT {
            input:
            sorted_fasta = SplitFastaBySeqLengthAndSortNT.sorted,
            accession2taxid = DownloadAccession2Taxid.accession2taxid,
            reads_by_taxid_path = "reads_by_taxid_nt",
            cpu = 64,
            docker_image_id = docker_image_id
        }
        call CompressNT {
            input:
            reads_by_taxid = BreakApartByTaxidNT.reads_by_taxid,
            k = nt_compression_k,
            scaled = nt_compression_scaled,
            similarity_threshold = nt_compression_similarity_threshold,
            logging_enabled = logging_enabled,
            docker_image_id = docker_image_id,
        }
    }

    if (skip_nuc_compression) {
        call GenerateNTDB as GenerateNTDBNoCompression {
            input:
            nt = DownloadNT.nt,
            docker_image_id = docker_image_id
        }
        call GenerateIndexMinimap2 as GenerateIndexMinimap2NoCompression {
            input:
            nt = DownloadNT.nt,
            docker_image_id = docker_image_id
        }
    }

    if (!skip_nuc_compression) {
        call GenerateNTDB as GenerateNTDBWCompression {
            input:
            nt = select_first([CompressNT.nt_compressed]),
            docker_image_id = docker_image_id
        }
        call GenerateIndexMinimap2 as GenerateIndexMinimap2WCompression {
            input:
            nt = select_first([CompressNT.nt_compressed]),
            docker_image_id = docker_image_id
        }
    }

    if (skip_protein_compression) {
        call GenerateNRDB as GenerateNRDBNoCompression {
            input:
            nr = DownloadNR.nr,
            docker_image_id = docker_image_id
        }
        call GenerateIndexDiamond as GenerateIndexDiamondNoCompression {
            input:
            nr = DownloadNR.nr,
            docker_image_id = docker_image_id
        }
    }
    if (!skip_protein_compression) {
        call GenerateNRDB as GenerateNRDBWCompression {
            input:
            nr = select_first([CompressNR.nr_compressed]),
            docker_image_id = docker_image_id
        }
        call GenerateIndexDiamond as GenerateIndexDiamondWCompression{
            input:
            nr = select_first([CompressNR.nr_compressed]),
            docker_image_id = docker_image_id
        }
    }

    if (skip_protein_compression && !skip_nuc_compression) {
        call GenerateIndexAccessions as GenerateIndexAccessionsNoCompressProtein {
            input:
            nr = DownloadNR.nr,
            nt = select_first([CompressNT.nt_compressed]),
            accession2taxid = DownloadAccession2Taxid.accession2taxid,
            docker_image_id = docker_image_id
        }
    }

    if (!skip_protein_compression && skip_nuc_compression) {
        call GenerateIndexAccessions as GenerateIndexAccessionsNoCompressNuc {
            input:
            nr = select_first([CompressNR.nr_compressed]),
            nt = DownloadNT.nt,
            accession2taxid = DownloadAccession2Taxid.accession2taxid,
            docker_image_id = docker_image_id
        }
    }
    if (!skip_nuc_compression && !skip_protein_compression) {
        call GenerateIndexAccessions as GenerateIndexAccessionsWCompression {
            input:
            nr = select_first([CompressNR.nr_compressed]),
            nt = select_first([CompressNT.nt_compressed]),
            accession2taxid = DownloadAccession2Taxid.accession2taxid,
            docker_image_id = docker_image_id
        }
    }

    if (skip_nuc_compression && skip_protein_compression) {
        call GenerateIndexAccessions as GenerateIndexAccessionsNoCompression {
            input:
            nr = DownloadNR.nr,
            nt = DownloadNT.nt,
            accession2taxid = DownloadAccession2Taxid.accession2taxid,
            docker_image_id = docker_image_id
        }
    }


    call GenerateIndexLineages {
        input:
        taxdump = DownloadTaxdump.taxdump,
        index_name = index_name,
        previous_lineages = previous_lineages,
        docker_image_id = docker_image_id
    }


    # if (write_to_db && defined(environ) && defined(s3_dir)) {
    call LoadTaxonLineages {
            input:
            environ = environ,
            index_name = index_name,
            s3_dir = s3_dir,
            minimap2_index = if (skip_nuc_compression) then GenerateIndexMinimap2NoCompression.minimap2_index else GenerateIndexMinimap2WCompression.minimap2_index,
            diamond_index = if (skip_protein_compression) then GenerateIndexDiamondNoCompression.diamond_index else GenerateIndexDiamondWCompression.diamond_index,
            nt = if (skip_nuc_compression) then DownloadNT.nt else CompressNT.nt_compressed,
            nr = if (skip_protein_compression) then DownloadNR.nr else CompressNR.nr_compressed,
            nt_loc_db = if (skip_nuc_compression) then GenerateNTDBNoCompression.nt_loc_db else GenerateNTDBWCompression.nt_loc_db,
            nt_info_db = if (skip_nuc_compression) then GenerateNTDBNoCompression.nt_info_db else GenerateNTDBWCompression.nt_info_db,
            nr_loc_db = if (skip_protein_compression) then GenerateNRDBNoCompression.nr_loc_db else GenerateNRDBWCompression.nr_loc_db,
            accession2taxid_db = if (skip_protein_compression && !skip_nuc_compression) then GenerateIndexAccessionsNoCompressProtein.accession2taxid_db
                        else if (!skip_protein_compression && skip_nuc_compression) then GenerateIndexAccessionsNoCompressNuc.accession2taxid_db
                        else if (!skip_nuc_compression && !skip_protein_compression)then GenerateIndexAccessionsWCompression.accession2taxid_db
                        else GenerateIndexAccessionsNoCompression.accession2taxid_db,
            docker_image_id = docker_image_id,

    }

    output {
        File? nr = LoadTaxonLineages.nr_
        File? nt = LoadTaxonLineages.nt_
        File? accession2taxid_db = LoadTaxonLineages.accession2taxid_db_
        File? nt_loc_db = LoadTaxonLineages.nt_loc_db_
        File? nt_info_db = LoadTaxonLineages.nt_info_db_
        File? nr_loc_db = LoadTaxonLineages.nr_loc_db_
        Directory? minimap2_index = LoadTaxonLineages.minimap2_index_
        Directory? diamond_index = LoadTaxonLineages.diamond_index_
        File taxid_lineages_db = GenerateIndexLineages.taxid_lineages_db
        File versioned_taxid_lineages_csv = GenerateIndexLineages.versioned_taxid_lineages_csv
        File deuterostome_taxids = GenerateIndexLineages.deuterostome_taxids
        File taxon_ignore_list = GenerateIndexLineages.taxon_ignore_list
        File changed_taxa_log = GenerateIndexLineages.changed_taxa_log
        File deleted_taxa_log = GenerateIndexLineages.deleted_taxa_log
        File new_taxa_log = GenerateIndexLineages.new_taxa_log
    }
}

task DownloadNR {
    input {
        String ncbi_server
        String docker_image_id
        String old_nr_s3_path
    }

    command <<<
        if ["~{old_nr_s3_path}" == ""]; then
            ncbi_download "~{ncbi_server}" blast/db/FASTA/nr.gz
        else
            aws s3 cp ~{old_nr_s3_path} blast/db/FASTA/nr
        fi
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
        String old_nt_s3_path
    }

    command <<<
        if ["~{old_nt_s3_path}" == ""]; then
            ncbi_download "~{ncbi_server}" blast/db/FASTA/nt.gz
        else
            aws s3 cp ~{old_nt_s3_path} blast/db/FASTA/nt
        fi
        
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
        String s3_accession_mapping_prefix
    }

    command <<<
        if ["~{s3_accession_mapping_prefix}" == ""]; then
            ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
            ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
            ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
            ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
        else
            aws s3 cp ~{s3_accession_mapping_prefix}/accession2taxid/nucl_gb.accession2taxid.gz pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
            aws s3 cp ~{s3_accession_mapping_prefix}/accession2taxid/nucl_wgs.accession2taxid.gz pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
            aws s3 cp ~{s3_accession_mapping_prefix}/accession2taxid/pdb.accession2taxid.gz pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
            aws s3 cp ~{s3_accession_mapping_prefix}/accession2taxid/prot.accession2taxid.gz pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

            gunzip pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
            gunzip pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
            gunzip pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
            gunzip pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
        fi

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
        # File nr_info_db = "nr_info.marisa"
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
        cd diamond_index_chunksize_~{chunksize}
        for dmnd_file in *.dmnd; do
            output=$(diamond dbinfo -d "$dmnd_file")

            # Extract the number of Letters and Sequences from the output
            letters=$(echo "$output" | grep "Letters" | awk '{print $NF}')
            sequences=$(echo "$output" | grep "Sequences" | awk '{print $NF}')

            # Extract the chunk identifier from the current dmnd file name
            chunk_number="${dmnd_file##*_}"
            chunk_number="${chunk_number%.dmnd}"

            # Rename the file
            mv "$dmnd_file" "${chunk_number}-${sequences}-${letters}.dmnd"
        done

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
        File changed_taxa_log = "changed_lineage_taxa.csv.gz"
        File deleted_taxa_log = "deleted_taxa.csv.gz"
        File new_taxa_log = "new_taxa.csv.gz"
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
        # File versioned_taxid_lineages_csv
        Directory? minimap2_index
        Directory? diamond_index
        File? nt
        File? nr
        File? nt_loc_db
        File? nt_info_db
        File? nr_loc_db
        File? accession2taxid_db
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

        # don't perform taxon lineage update
        # gzip -dc versioned_taxid_lineages_csv > "taxon_lineages_new.csv"
        # COLS=$(head -n 1 "taxon_lineages_new.csv")
        # mysql --defaults-extra-file=my.cnf -D "$DATABASE" -e "CREATE TABLE taxon_lineages_new LIKE taxon_lineages"
        # mysqlimport --defaults-extra-file=my.cnf --verbose --local --columns="$COLS" --fields-terminated-by=',' --fields-optionally-enclosed-by='"' --ignore-lines 1 "$DATABASE" "taxon_lineages_new.csv"
        # mysql --defaults-extra-file=my.cnf -D "$DATABASE" -e "RENAME TABLE taxon_lineages TO taxon_lineages_old, taxon_lineages_new To taxon_lineages"

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
                minimap2_short_db_path,
                diamond_db_path,
                lineage_version,
                created_at,
                updated_at
            ) VALUES(
                '~{index_name}',
                '~{s3_dir}/index-generation-2/nt_compressed.fa',
                '~{s3_dir}/index-generation-2/nt_loc.marisa',
                '~{s3_dir}/index-generation-2/nr_compressed.fa',
                '~{s3_dir}/index-generation-2/nr_loc.marisa',
                's3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/taxid-lineages.marisa',
                '~{s3_dir}/index-generation-2/accession2taxid.marisa', # patch for now
                's3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/deuterostome_taxids.txt',
                '~{s3_dir}/index-generation-2/nt_info.marisa',
                's3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/taxon_ignore_list.txt',
                '~{s3_dir}/index-generation-2/nt_k14_w8_20/',
                '~{s3_dir}/index-generation-2/diamond_index_chunksize_5500000000/',
                '2021-01-22',
                NOW(),
                NOW()
            ); 
        "
    >>>

    output {
        File? nt_=nt
        File? nr_=nr
        File? accession2taxid_db_=accession2taxid_db
        File? nt_loc_db_=nt_loc_db
        File? nt_info_db_=nt_info_db
        File? nr_loc_db_=nr_loc_db
        Directory? minimap2_index_=minimap2_index
        Directory? diamond_index_=diamond_index
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

task SplitFastaBySeqLengthAndSort {
    input {
        File fasta
        String combined_sorted_path
        Int cpu
        Int threads = if cpu * 0.6 < 1 then 1 else floor(cpu * 0.6)
        String docker_image_id
    }

    command <<<
        # Sort NT/NR by length with the longer sequences first
        #   This is needed because the downstream compression algorithm iterates through NT/NR in order only emitting
        #   sequences if they are not contained by what it has already seen. If a shorter sequence is
        #   contained by a longer sequence, and the shorter sequence were to come first, it would be emitted
        #   even though it is redundant to the longer sequence.
        set -euxo pipefail

        total_seqs=$(grep ">" ~{fasta} | wc -l)
        ncbi-compress  break-up-fasta-by-sequence-length  \
            --input-fasta ~{fasta}  \
            --output-dir outputs  \
            --total-sequence-count ${total_seqs}  \
            --bin-size 25 \

        # apt-get install -y parallel
        # parallel -j ~{threads} 'seqkit sort --reverse --by-length --two-pass {} -o sorted_{};' ::: outputs/*.fa

        # cd outputs
        # for file in *.fa; do
        #     seqkit sort --reverse --by-length --two-pass --threads ~{threads} $file -o sorted_$file
        #     # rm $file
        #     # rm $file.seqkit.fai
        #     echo "finished sorting $file"
        # done
        # cd ..

        # Combine the sorted files with longest sequences at the top
        ls -r outputs/*.fa | xargs cat > ~{combined_sorted_path}
    >>>
    output {
        File sorted = combined_sorted_path
    }
    runtime {
        docker: docker_image_id
        cpu: cpu
    }
}

task BreakApartByTaxid {
    input {
        File sorted_fasta
        Directory accession2taxid
        String reads_by_taxid_path
        Int cpu
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        ncbi-compress break-into-individual-taxids-only \
            --input-fasta ~{sorted_fasta} \
            --accession-mapping-files \
                ~{accession2taxid}nucl_wgs.accession2taxid \
                ~{accession2taxid}nucl_gb.accession2taxid \
                ~{accession2taxid}pdb.accession2taxid \
            --output-dir ~{reads_by_taxid_path}
    >>>

    output {
        Directory reads_by_taxid = reads_by_taxid_path
    }

    runtime {
        docker: docker_image_id
        cpu: cpu
    }
}

task CompressNT {
    input {
        Directory reads_by_taxid
        Int k
        Int scaled
        Float similarity_threshold
        Boolean logging_enabled
        String docker_image_id
    }

    command <<< 
        set -euxo pipefail

        if [ "~{logging_enabled}" ]; then
            ncbi-compress fasta-compress-from-taxid-dir \
                --input-fasta-dir ~{reads_by_taxid} \
                --output-fasta nt_compressed.fa \
                --k ~{k} \
                --scaled ~{scaled} \
                --similarity-threshold ~{similarity_threshold} \
                --enable-sequence-retention-logging \
                --logging-contained-in-tree-fn nt_contained_in_tree.tsv \
                --logging-contained-in-chunk-fn nt_contained_in_chunk.tsv
        else
            ncbi-compress fasta-compress-from-taxid-dir \
                --input-fasta-dir ~{reads_by_taxid} \
                --output-fasta nt_compressed.fa \
                --k ~{k} \
                --scaled ~{scaled} \
                --similarity-threshold ~{similarity_threshold}
        fi
    >>>

    output {
        File nt_compressed = "nt_compressed.fa"
        File nt_contained_in_tree = "nt_contained_in_tree.tsv"
        File nt_contained_in_chunk = "nt_contained_in_chunk.tsv"
    }

    runtime {
        docker: docker_image_id
        cpu: 64
        memory: "488G"
    }
}

task CompressNR {
    input {
        Directory reads_by_taxid
        Int k
        Int scaled
        Float similarity_threshold
        Boolean logging_enabled
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        if [ "~{logging_enabled}" ]; then
            ncbi-compress fasta-compress-from-taxid-dir \
                --input-fasta-dir ~{reads_by_taxid} \
                --output-fasta nr_compressed.fa \
                --k ~{k} \
                --scaled ~{scaled} \
                --similarity-threshold ~{similarity_threshold} \
                --is-protein-fasta \
                --enable-sequence-retention-logging \
                --logging-contained-in-tree-fn nr_contained_in_tree.tsv \
                --logging-contained-in-chunk-fn nr_contained_in_chunk.tsv
        else
            ncbi-compress fasta-compress-from-taxid-dir \
                --input-fasta-dir ~{reads_by_taxid} \
                --output-fasta nr_compressed.fa \
                --k ~{k} \
                --scaled ~{scaled} \
                --similarity-threshold ~{similarity_threshold} \
                --is-protein-fasta
        fi
    >>>

    output {
        File nr_compressed = "nr_compressed.fa"
        File nr_contained_in_tree = "nr_contained_in_tree.tsv"
        File nr_contained_in_chunk = "nr_contained_in_chunk.tsv"
    }

    runtime {
        docker: docker_image_id
        cpu: 64
        memory: "488G"
    }
}
