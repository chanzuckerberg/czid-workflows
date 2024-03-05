version development


workflow index_generation {
    input {
        String index_name
        String ncbi_server = "https://ftp.ncbi.nih.gov"
        File? previous_lineages

        # Compression Parameters
        Int nt_compression_k = 31
        Int nt_compression_scaled = 1000
        Float nt_compression_similarity_threshold = 0.5

        Int nr_compression_k = 31
        Int nr_compression_scaled = 100
        Float nr_compression_similarity_threshold = 0.5

        Boolean skip_protein_compression = false
        Boolean skip_generate_nr_assets = false
        Boolean skip_nuc_compression = false
        Boolean skip_generate_nt_assets = false
        Boolean logging_enabled = false

        String? provided_nt
        String? provided_nr
        String provided_accession2taxid_nucl_gb = "~{ncbi_server}/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
        String provided_accession2taxid_nucl_wgs = "~{ncbi_server}/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"
        String provided_accession2taxid_pdb = "~{ncbi_server}/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz"
        String provided_accession2taxid_prot = "~{ncbi_server}/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz"
        String provided_taxdump = "~{ncbi_server}/pub/taxonomy/taxdump.tar.gz"

        String docker_image_id
    }

    Array[String?] possibly_zipped_files = [
        provided_accession2taxid_nucl_gb,
        provided_accession2taxid_nucl_wgs,
        provided_accession2taxid_pdb,
        provided_accession2taxid_prot,
        provided_taxdump,
        provided_nt,
        provided_nr
    ]

    # Download files if they are not provided
    scatter (file in possibly_zipped_files) {
        if (defined(file) && select_first([file]) != "") {
            # if filename ends with gz
            String file_ = select_first([file]) # this is safe because we know it's defined and not empty
            if (sub(basename(file_), "\\.gz$", "") != basename(file_)) {
                call UnzipFile {
                    input:
                    zipped_file = file_,
                    cpu = 8,
                    docker_image_id = docker_image_id,
                }
            }
        }

        
        File unzipped_file = select_first([UnzipFile.file, file])
    }
    File unzipped_accession2taxid_nucl_gb = unzipped_file[0]
    File unzipped_accession2taxid_nucl_wgs = unzipped_file[1]
    File unzipped_accession2taxid_pdb = unzipped_file[2]
    File unzipped_accession2taxid_prot = unzipped_file[3]
    File unzipped_taxdump = unzipped_file[4]
    File provided_nt_unzipped = unzipped_file[5]
    File provided_nr_unzipped = unzipped_file[6]

    Boolean is_nt_provided = defined(provided_nt_unzipped)
    Boolean is_nr_provided = defined(provided_nr_unzipped)

    if (!is_nt_provided) {
        call DownloadDatabase as DownloadNT {
            input:
            database_type = "nt",
            docker_image_id = docker_image_id
        }
    }
    if (!is_nr_provided) {
        call DownloadDatabase as DownloadNR {
            input:
            database_type = "nr",
            docker_image_id = docker_image_id
        }
    }

    File nt_download = select_first([provided_nt_unzipped, DownloadNT.database])
    File nr_download = select_first([provided_nr_unzipped, DownloadNR.database])

    if (!skip_protein_compression) {
        call CompressDatabase as CompressNR {
            input:
            database_type = "nr",
            fasta = nr_download,
            accession2taxid_files = [unzipped_accession2taxid_pdb, unzipped_accession2taxid_prot],
            k = nr_compression_k,
            scaled = nr_compression_scaled,
            similarity_threshold = nr_compression_similarity_threshold,
            logging_enabled = logging_enabled,
            docker_image_id = docker_image_id,
        }
    }
    File nr_or_compressed = select_first([CompressNR.compressed, nr_download])

    if (!skip_nuc_compression) {
        call CompressDatabase as CompressNT {
            input:
            database_type = "nt",
            fasta = nt_download,
            accession2taxid_files = [unzipped_accession2taxid_nucl_wgs, unzipped_accession2taxid_nucl_gb],
            k = nt_compression_k,
            scaled = nt_compression_scaled,
            similarity_threshold = nt_compression_similarity_threshold,
            logging_enabled = logging_enabled,
            docker_image_id = docker_image_id,
        }
    }
    File nt_or_compressed = select_first([CompressNT.compressed, nt_download])

    if (!skip_generate_nt_assets) {
        call GenerateLocDB as GenerateNTLocDB {
            input:
            db_fasta = nt_or_compressed,
            database_type = "nt",
            docker_image_id = docker_image_id
        }

        call GenerateInfoDB as GenerateNTInfoDB {
            input:
            db_fasta = nt_or_compressed,
            database_type = "nt",
            docker_image_id = docker_image_id
        }

        call GenerateIndexMinimap2 as GenerateIndexMinimap2 {
            input:
            nt = nt_or_compressed,
            docker_image_id = docker_image_id
        }
    }

    if (!skip_generate_nr_assets) {
        call GenerateLocDB as GenerateNRLocDB {
            input:
            db_fasta = nr_or_compressed,
            database_type = "nr",
            docker_image_id = docker_image_id
        }

        call GenerateIndexDiamond as GenerateIndexDiamond {
            input:
            nr = nr_or_compressed,
            docker_image_id = docker_image_id
        }
    }

    call GenerateIndexAccessions {
        input:
        nr = nr_or_compressed,
        nt = nt_or_compressed,
        accession2taxid_files = [
            unzipped_accession2taxid_nucl_wgs,
            unzipped_accession2taxid_nucl_gb,
            unzipped_accession2taxid_pdb,
            unzipped_accession2taxid_prot,
        ],
        docker_image_id = docker_image_id
    }

    call GenerateIndexLineages {
        input:
        taxdump = unzipped_taxdump,
        index_name = index_name,
        previous_lineages = previous_lineages,
        docker_image_id = docker_image_id
    }

    output {
        File nr = nr_or_compressed
        File nt = nt_or_compressed
        File accession2taxid_nucl_wgs = unzipped_accession2taxid_nucl_wgs
        File accession2taxid_nucl_gb = unzipped_accession2taxid_nucl_gb
        File accession2taxid_pdb = unzipped_accession2taxid_pdb
        File accession2taxid_prot = unzipped_accession2taxid_prot
        File? nt_loc_db = GenerateNTLocDB.loc_db
        File? nt_info_db = GenerateNTInfoDB.info_db
        File? nr_loc_db = GenerateNRLocDB.loc_db
        Directory? minimap2_index = GenerateIndexMinimap2.minimap2_index
        Directory? diamond_index = GenerateIndexDiamond.diamond_index
        File taxid_lineages_db = GenerateIndexLineages.taxid_lineages_db
        File versioned_taxid_lineages_csv = GenerateIndexLineages.versioned_taxid_lineages_csv
        File deuterostome_taxids = GenerateIndexLineages.deuterostome_taxids
        File taxon_ignore_list = GenerateIndexLineages.taxon_ignore_list
        File changed_taxa_log = GenerateIndexLineages.changed_taxa_log
        File deleted_taxa_log = GenerateIndexLineages.deleted_taxa_log
        File new_taxa_log = GenerateIndexLineages.new_taxa_log
        File? nr_contained_in_tree = CompressNR.contained_in_tree
        File? nr_contained_in_chunk = CompressNR.contained_in_chunk
        File? nt_contained_in_tree = CompressNT.contained_in_tree
        File? nt_contained_in_chunk = CompressNT.contained_in_chunk
        File accession2taxid_db = GenerateIndexAccessions.accession2taxid_db
    #     Directory? nt_split_apart_taxid_dir = CompressNT.split_apart_taxid_dir
    #     Directory? nt_sorted_taxid_dir = CompressNT.sorted_taxid_dir
    #     Directory? nr_split_apart_taxid_dir = CompressNR.split_apart_taxid_dir
    #     Directory? nr_sorted_taxid_dir = CompressNR.sorted_taxid_dir
    }
}

task DownloadDatabase {
    input {
        String database_type # nt or nr
        String docker_image_id
        Int threads = 64
    }

    command <<<
        set -euxo pipefail

        update_blastdb.pl --decompress ~{database_type} --num_threads ~{threads}
        blastdbcmd -db ~{database_type} -entry all -out ~{database_type}.fsa

    >>>

    output {
        File database = "~{database_type}.fsa"
    }

    runtime {
        docker: docker_image_id
        cpu: 64
    }

}

task UnzipFile {
    input {
        String zipped_file
        Int cpu
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        output="~{basename(zipped_file)}"
        curl ~{zipped_file} -o $output
        pigz -p ~{cpu} -dc $output > ~{sub(basename(zipped_file), "\\.gz$", "")}

	#pigz -p ~{cpu} -dc ~{zipped_file} > ~{sub(basename(zipped_file), "\\.gz$", "")}
    >>>

    output {
        File file = sub(basename(zipped_file), "\\.gz$", "")
    }

    runtime {
        docker: docker_image_id
        cpu: cpu
    }
}

task GenerateIndexAccessions {
    input {
        File nr
        File nt
        Array[File] accession2taxid_files
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        # Build index
        python3 /usr/local/bin/generate_accession2taxid.py ~{sep=" " accession2taxid_files} \
            --nt_file ~{nt} \
            --nr_file ~{nr} \
            --accession2taxid_db accession2taxid.marisa
    >>>

    output {
        File accession2taxid_db = "accession2taxid.marisa"
    }

    runtime {
        docker: docker_image_id
        cpu: 8
    }
}

task GenerateLocDB {
    input {
        File db_fasta
        String database_type # nt or nr
        String docker_image_id
    }
    command <<<
        python3 /usr/local/bin/generate_ncbi_db_index.py loc ~{db_fasta} ~{database_type}_loc.marisa
    >>>
    output {
        File loc_db = "~{database_type}_loc.marisa"
    }
    runtime {
        docker: docker_image_id
    }
}

task GenerateInfoDB {
    input {
        File db_fasta
        String database_type # nt or nr
        String docker_image_id
    }
    command <<<
        python3 /usr/local/bin/generate_ncbi_db_index.py info ~{db_fasta} ~{database_type}_info.marisa
    >>>
    output {
        File info_db = "~{database_type}_info.marisa"
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
        cpu: 16
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
        cpu: 8
    }
}

task GenerateIndexMinimap2 {
    input {
        File nt
        Int k = 14 # Minimizer k-mer length default is 21 for short reads option
        Int w = 8 # Minimizer window size default is 11 for short reads option
        String I = "9999G" # Load at most NUM target bases into RAM for indexing
        Int t = 20 # number of threads, doesn't really work for indexing I don't think
        Int n_chunks = 50
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
        cpu: 16
    }
}

task CompressDatabase {
    input {
        String database_type = "nt" # nt or nr
        File fasta
        Array[File] accession2taxid_files
        Int k
        Int scaled
        Float similarity_threshold
        Boolean logging_enabled
        String docker_image_id
    }

    command <<< 
        set -euxo pipefail

        READS_BY_TAXID_PATH=reads_by_taxid
        SPLIT_APART_TAXID_DIR_NAME="split_apart_taxid_~{database_type}"
        SORTED_TAXID_DIR_NAME="sorted_taxid_~{database_type}"

        # It is critical that this split happens in the same step as the compression
        # If the directory is a step output it will be uploaded, which takes an enormous amount of time
        ncbi-compress break-into-individual-taxids-only \
            --input-fasta ~{fasta} \
            --accession-mapping-files ~{sep=" " accession2taxid_files} \
            --output-dir $READS_BY_TAXID_PATH

        ncbi-compress sort-taxid-dir-by-sequence-length \
            --input-taxid-dir $READS_BY_TAXID_PATH \
            --output-taxid-dir $SORTED_TAXID_DIR_NAME

        mkdir $SPLIT_APART_TAXID_DIR_NAME

        if [ "~{logging_enabled}" == "true" ]; then
            ncbi-compress fasta-compress-from-taxid-dir ~{if database_type == "nr" then "--is-protein-fasta" else ""} \
                --input-fasta-dir $SORTED_TAXID_DIR_NAME \
                --output-fasta ~{database_type}_compressed.fa \
                --k ~{k} \
                --scaled ~{scaled} \
                --similarity-threshold ~{similarity_threshold} \
                --split-apart-taxid-dir-name $SPLIT_APART_TAXID_DIR_NAME \
                --enable-sequence-retention-logging \
                --logging-contained-in-tree-fn ~{database_type}_contained_in_tree.tsv \
                --logging-contained-in-chunk-fn ~{database_type}_contained_in_chunk.tsv
        else
            ncbi-compress fasta-compress-from-taxid-dir  ~{if database_type == "nr" then "--is-protein-fasta" else ""} \
                --input-fasta-dir $SORTED_TAXID_DIR_NAME \
                --output-fasta ~{database_type}_compressed.fa \
                --k ~{k} \
                --scaled ~{scaled} \
                --similarity-threshold ~{similarity_threshold} \
                --split-apart-taxid-dir-name $SPLIT_APART_TAXID_DIR_NAME
        fi

        # shuffle compressed fasta to distribute the accessions evenly across the file
        # this is important for spreading SC2 accessions (and any other large taxid) over
        # the chunked minimap2 and diamond indexes which impacts alignment time.
        ncbi-compress shuffle-fasta \
            --input-fasta ~{database_type}_compressed.fa \
            --output-fasta ~{database_type}_compressed_shuffled.fa

        # Remove to save space, intermediate files are not cleaned up within a run
        rm -rf $READS_BY_TAXID_PATH
        rm -rf $SPLIT_APART_TAXID_DIR_NAME
    >>>

    output {
        File compressed = "~{database_type}_compressed_shuffled.fa"
        File? contained_in_tree = "~{database_type}_contained_in_tree.tsv"
        File? contained_in_chunk = "~{database_type}_contained_in_chunk.tsv"
        # Directory sorted_taxid_dir = "sorted_taxid_~{database_type}"
    }

    runtime {
        docker: docker_image_id
        cpu: 88
    }
}

