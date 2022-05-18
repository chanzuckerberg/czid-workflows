version development

workflow index_generation {
    input {
        String index_name
        String ncbi_server = "ftp://ftp.ncbi.nih.gov"
        File? previous_lineages = ""
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

    call DownloadTaxdump {
        input:
        ncbi_server = ncbi_server,
        docker_image_id = docker_image_id
    }

    call GenerateIndexAccessions {
        input:
        nr = DownloadNR.nr,
        nt = DownloadNT.nt,
        accession2taxid = DownloadAccession2Taxid.accession2taxid,
        docker_image_id = docker_image_id
    }

    call GenerateNTDB {
        input:
        nt = DownloadNT.nt,
        docker_image_id = docker_image_id
    }

    call GenerateNRDB {
        input:
        nr = DownloadNR.nr,
        docker_image_id = docker_image_id
    }

    call GenerateIndexLineages {
        input:
        taxdump = DownloadTaxdump.taxdump,
        index_name = index_name,
        docker_image_id = docker_image_id
    }

    call ChunkNT {
        input:
        nt = DownloadNT.nt,
        docker_image_id = docker_image_id
    }

    scatter (nt_chunk in ChunkNT.nt_chunks) {
        call GenerateIndexMinimap2Chunk {
            input:
            nt_chunk = nt_chunk,
            docker_image_id = docker_image_id
        }
    }

    call ChunkNR {
        input:
        nr = DownloadNR.nr,
        docker_image_id = docker_image_id
    }

    scatter (nr_chunk in ChunkNR.nr_chunks) {
        call GenerateIndexDiamondChunk {
            input:
            nr = DownloadNR.nr,
            docker_image_id = docker_image_id
        }
    }

    output {
        File nr = DownloadNR.nr
        File nt = DownloadNT.nt
        Directory accession2taxid = DownloadAccession2Taxid.accession2taxid
        File taxdump = DownloadTaxdump.taxdump
        File accession2taxid_gz = GenerateIndexAccessions.accession2taxid_gz
        File accession2taxid_wgs = GenerateIndexAccessions.accession2taxid_wgs
        File accession2taxid_db = GenerateIndexAccessions.accession2taxid_db
        File taxid2wgs_accession_db = GenerateIndexAccessions.taxid2wgs_accession_db
        File nt_loc_db = GenerateNTDB.nt_loc_db
        File nt_info_db = GenerateNTDB.nt_info_db
        File nr_loc_db = GenerateNRDB.nr_loc_db
        File nr_info_db = GenerateNRDB.nr_info_db
        File taxid_lineages_db = GenerateIndexLineages.taxid_lineages_db
        File taxid_lineages_csv = GenerateIndexLineages.taxid_lineages_csv
        File names_csv = GenerateIndexLineages.names_csv
        File named_taxid_lineages_csv = GenerateIndexLineages.named_taxid_lineages_csv
        File versioned_taxid_lineages_csv = GenerateIndexLineages.versioned_taxid_lineages_csv
        File deuterostome_taxids = GenerateIndexLineages.deuterostome_taxids
        File taxon_ignore_list = GenerateIndexLineages.taxon_ignore_list
        Array[File] minimap2_index = GenerateIndexMinimap2Chunk.minimap2_index
        Array[File] diamond_index = GenerateIndexDiamondChunk.diamond_index
    }
}

task DownloadNR {
    input {
        String ncbi_server
        String docker_image_id
    }

    command <<<
        ncbi_download "~{ncbi_server}" blast/db/FASTA/nr.gz
        pigz -dc blast/db/FASTA/nr.gz > nr
    >>>

    output {
        File nr = "nr"
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
        pigz -dc blast/db/FASTA/nt.gz > nt
        
    >>>

    output {
        File nt = "nt"
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
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/pdb.accession2taxid.gz
        ncbi_download "~{ncbi_server}" pub/taxonomy/accession2taxid/prot.accession2taxid.gz
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
        Int parallelism = 16
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
            --parallelism ~{parallelism} \
            --nt_file ~{nt} \
            --nr_file ~{nr} \
            --output_gz accession2taxid.gz \
            --output_wgs_gz accession2taxid_wgs.gz \
            --accession2taxid_db accession2taxid.db \
            --taxid2wgs_accession_db taxid2wgs_accession.db
    >>>

    output {
        File accession2taxid_gz = "accession2taxid.gz"
        File accession2taxid_wgs = "accession2taxid_wgs.gz"
        File accession2taxid_db = "accession2taxid.db"
        File taxid2wgs_accession_db = "taxid2wgs_accession.db"
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
        python3 /usr/local/bin/generate_loc_db.py ~{nt} nt_loc.db nt_info.db
    >>>

    output {
        File nt_loc_db = "nt_loc.db"
        File nt_info_db = "nt_info.db"
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
        python3 /usr/local/bin/generate_loc_db.py ~{nr} nr_loc.db nr_info.db
    >>>

    output {
        File nr_loc_db = "nr_loc.db"
        File nr_info_db = "nr_info.db"
    }

    runtime {
        docker: docker_image_id
    }
}

task ChunkNR {
    input {
        File nr
        Int n_chunks = 20
        String docker_image_id
    }

    command <<<
        seqkit split2 ~{nr} -p ~{n_chunks} --out-dir nr.split
    >>>

    output {
        Array[File] nr_chunks = glob("nr.split/*")
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexDiamondChunk {
    input {
        File nr_chunk
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        chunk_path="~{nr_chunk}"
        chunk_number="${chunk_path##*_}"
        # Ignore warning is needed because sometimes NR has sequences of only DNA characters which causes this to fail
        diamond makedb --ignore-warnings --in ~{nr_chunk} --db "diamond_index_part_${chunk_number}"
    >>>

    output {
        File diamond_index = glob("diamond_index_part_*")[0]
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
        File deuterostome_taxids = "ncbitax2lin/deuterostome_taxids.txt"
        File taxon_ignore_list = "ncbitax2lin/taxon_ignore_list.txt"
    }

    runtime {
        docker: docker_image_id
    }
}

task ChunkNT {
    input {
        File nt
        Int n_chunks = 20
        String docker_image_id
    }

    command <<<
        seqkit split2 ~{nt} -p ~{n_chunks} --out-dir nt.split
    >>>

    output {
        Array[File] nt_chunks = glob("nt.split/*")
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateIndexMinimap2Chunk {
    input {
        File nt_chunk
        Int k = 14 # Minimizer k-mer length default is 21 for short reads option
        Int w = 8 # Minimizer window size default is 11 for short reads option
        String I = "9999G" # Load at most NUM target bases into RAM for indexing
        Int t = 20 # number of threads, doesn't really work for indexing I don't think
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        chunk_path="~{nt_chunk}"
        chunk_number="${chunk_path##*_}"
        minimap2 -cx sr -k ~{k} -w ~{w} -I ~{I} -t ~{t} -d "nt.part_"$chunk_number".idx" ~{nt_chunk}
    >>>

    output {
        File minimap2_index = glob("nt.part_*.idx")[0]
    }

    runtime {
        docker: docker_image_id
    }
}
