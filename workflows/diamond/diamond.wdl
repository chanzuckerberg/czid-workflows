version 1.1

workflow diamond {
    input {
        File query_0
        File? query_1
        String extra_args = "--mid-sensitive"
        File db_chunk
        String docker_image_id
    }

    call RunDiamond {
        input:
        query_0 = query_0,
        query_1 = query_1,
        extra_args = extra_args,
        db_chunk = db_chunk,
        docker_image_id = docker_image_id
    }
    output {} # explicitly mark empty output for intermediate_output tagging
}

task RunDiamond {
    input {
        File query_0
        File? query_1
        String extra_args
        File db_chunk
        String docker_image_id
    }

    command <<<
        python3 /usr/local/bin/diamond_scatter.py blastx-chunk --db ~{db_chunk} --query ~{query_0} ~{if defined(query_1) then '--query ~{query_1}' else ''} --out-dir chunks --diamond-args="~{extra_args}"
    >>>

    output {
        Array[File] ref_blocks = glob("chunks/ref_block*")
        Array[File] ref_dicts = glob("chunks/ref_dict*")
    }

    runtime {
        docker: docker_image_id
    }
}
