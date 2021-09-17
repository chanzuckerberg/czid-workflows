version 1.1

workflow alignment_scalability {
    input {
        File fastqs_0
        File? fastqs_1
        String input_dir
        String chunk_dir
        String db_path
        String docker_image_id

        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call RunDiamond {
        input:
        fastqs = select_all([fastqs_0, fastqs_1]),
        input_dir = input_dir,
        chunk_dir = chunk_dir,
        db_path = db_path,
        docker_image_id = docker_image_id
    }


    output {
        File out_tsv = RunDiamond.out_tsv
    }
}

task RunDiamond {
    input {
        Array[File]+ fastqs
        String input_dir
        String chunk_dir
        String db_path
        String docker_image_id
    }

    command <<<
        echo STARTING
        python3 <<CODE
        from idseq_utils.run_diamond import run_diamond
        fastqs = ["~{sep='", "' fastqs}"]
        run_diamond("~{input_dir}", "~{chunk_dir}", "~{db_path}", "out.tsv", *fastqs)
        CODE
    >>>

    output {
        File out_tsv = "out.tsv"
    }

    runtime {
        docker: docker_image_id
    }
}
