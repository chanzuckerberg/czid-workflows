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
        File output_kma = "output_kma.rgi.kma"
    }

    runtime {
        docker: docker_image_id
    }
}
