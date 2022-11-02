version development
# This WDL drives the host_filter_indexing.wdl workflow based on a JSON file specifying the raw
# data sources for all CZ ID host species (genome FASTA, transcripts GTF, transcripts FASTA).
import "host_filter_indexing.wdl" as main

struct HostRawData {
    File genome_fasta_gz
    File? transcripts_gtf_gz
    Array[File] transcripts_fasta_gz
}

struct HostIndex {
    File bowtie2_index_tar
    File hisat2_index_tar
    File kallisto_idx
    File minimap2_mmi

    File original_genome_fasta_gz
    File? original_transcripts_gtf_gz
    Array[File] original_transcripts_fasta_gz
    File original_ERCC_fasta_gz
    Array[File] original_other_fasta_gz
}

workflow host_filter_indexing_driver {
    input {
        File raw_data
        File ERCC_fasta_gz
        String docker

        # Populate in order to index only the given set of hosts. Otherwise (default) generate all.
        Array[String] host = []
    }

    Map[String,HostRawData] raw_data_map = read_json(raw_data)

    Array[String] hosts = if length(host) > 0 then host else keys(raw_data_map)   
    scatter (host1 in hosts) {
        HostRawData host1data = raw_data_map[host1]
        call main.host_filter_indexing as run {
            input:
            genome_name = host1,
            genome_fasta_gz = host1data.genome_fasta_gz,
            transcripts_gtf_gz = host1data.transcripts_gtf_gz,
            transcripts_fasta_gz = host1data.transcripts_fasta_gz,
            ERCC_fasta_gz, docker
        }
        Pair[String,HostIndex] entry = (host1, HostIndex {
            bowtie2_index_tar: run.bowtie2_index_tar,
            hisat2_index_tar: run.hisat2_index_tar,
            kallisto_idx: run.kallisto_idx,
            minimap2_mmi: run.minimap2_mmi,
            original_genome_fasta_gz: run.original_genome_fasta_gz,
            original_transcripts_gtf_gz: run.original_transcripts_gtf_gz,
            original_transcripts_fasta_gz: run.original_transcripts_fasta_gz,
            original_ERCC_fasta_gz: run.original_ERCC_fasta_gz,
            original_other_fasta_gz: run.original_other_fasta_gz
        })
    }

    output {
        Map[String,HostIndex] indexes = as_map(entry)
    }
}
