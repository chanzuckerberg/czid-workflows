version 1.1

workflow AMR {
    input {
        Array[File] non_host_reads
	File contigs
        String docker_image_id
        File card_json = "s3://czid-public-references/test/AMRv2/card.json"

        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call RunRgiBwtKma {
        input:
        non_host_reads = non_host_reads,
        card_json = card_json, 
        docker_image_id = docker_image_id
    }
    call RunRgiMain {
        input:
        contigs = contigs,
        card_json = card_json, 
        docker_image_id = docker_image_id
    }
    call RunRgiKmerBwt { 
        input:
        output_sorted_length_100 = RunRgiBwtKma.output_sorted_length_100,
        card_json = card_json,
        docker_image_id = docker_image_id
    }
    call RunRgiKmerMain { 
        input:
        main_output_json = RunRgiMain.output_json,
        card_json = card_json, 
        docker_image_id = docker_image_id
    }

}
task RunRgiKmerMain {
    input {
        File main_output_json
        File card_json
        String docker_image_id
    }
    command <<< 
	source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        rgi card_annotation -i "~{card_json}" > card_annotation.log
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta
        rgi kmer_query --rgi -k 61 -i "~{main_output_json}" --output output.rgi.main.kmerspecies --local
    >>>
    output { 
        Array[File] output_kmer_main = glob("output.rgi.main.kmerspecies*")
        File species_calling = "output.rgi.main.kmerspecies_61mer_analysis_rgi_summary.txt"
    }
    runtime {
        docker: docker_image_id
    }

}

task RunRgiKmerBwt { 
    input {
        File output_sorted_length_100
        File card_json
        String docker_image_id
    }
    command <<<
	source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        rgi card_annotation -i "~{card_json}" > card_annotation.log
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta
	rgi kmer_query --bwt -k 61 -i "~{output_sorted_length_100}" --output output.rgi.kma.kmerspecies --local
	>>>
    output {
        Array[File] output_kmer_bwt = glob("output.rgi.kma.kmerspecies*")
        File kma_species_calling = "output.rgi.kma.kmerspecies_61mer_analysis.gene.txt"
    }
    runtime {
        docker: docker_image_id
    }

}
task RunRgiMain { 
    input { 
        File contigs
        File card_json
        String docker_image_id
    }
    command <<<
        source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        rgi card_annotation -i "~{card_json}" > card_annotation.log
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta
        rgi main -i "~{contigs}" -o output.rgi.main -t contig -a BLAST --clean

    >>>
    output {
        Array[File] output_main = glob("output.rgi.main*")
        File output_json = "output.rgi.main.json"
        File main_amr_results = "output.rgi.main.txt"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunRgiBwtKma {
    input {
        Array[File] non_host_reads
        File card_json
        String docker_image_id
    }

    command <<<
	source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        rgi card_annotation -i "~{card_json}" > card_annotation.log 
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta  
    	rgi bwt -1 "~{non_host_reads[0]}" -2 "~{non_host_reads[0]}" -a kma -o output_kma.rgi.kma --clean
    >>>

    output {
        Array[File] output_kma = glob("output_kma.rgi.kma*")
        File output_sorted_length_100 = "output_kma.rgi.kma.sorted.length_100.bam"
        File kma_amr_results = "output_kma.rgi.kma.allele_mapping_data.txt"
    }

    runtime {
        docker: docker_image_id
    }
}
