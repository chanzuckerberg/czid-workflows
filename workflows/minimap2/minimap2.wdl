version 1.1

workflow minimap {
    input {
        File query_0
        File? query_1
        String extra_args = "-cx sr"
        File db_chunk
        String docker_image_id
    }

    call RunMinimap2 {
        input:
        query_0 = query_0,
        query_1 = query_1,
        extra_args = extra_args,
        db_chunk = db_chunk,
        docker_image_id = docker_image_id
    }

    output {
        Array[File] chunks = RunMinimap2.chunks
    }
}

task RunMinimap2 {
    input {
        File query_0
        File? query_1
        String extra_args
        File db_chunk
        String docker_image_id
    }

    command <<<
        mkdir chunks
        # example: nt.part_001 -< 001
        CHUNK=$(echo ~{db_chunk} | sed s/.*nt\.part_//)
        CPUS=$(nproc)

        minimap2 ~{extra_args} -t "$CPUS" --split-map "chunks/intermediate${CHUNK}" ~{db_chunk} ~{query_0} ~{if defined(query_1) then '~{query_1}' else ''}
    >>>

    output {
        Array[File] chunks = glob("chunks/*.idx")
    }

    runtime {
        docker: docker_image_id
    }
}
