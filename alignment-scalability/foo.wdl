version 1.1

workflow hello_world {
    input {
        File my_file
    }

    call AddStr as add_hello {
        input:
        in_file = my_file,
        str = "hello",
        docker_image_id = docker_image_id
    }

    call AddStr as add_world {
        input:
        input_file = add_hello.out_file,
        str = "world",
        docker_image_id = docker_image_id
    }


    output {
        File hello_world_file = add_world.out_file
    }
}

task AddStr {
    input {
        File in_file
        String str
        String docker_image_id
    }

    command <<<
        cat ~{input_file} > out_file
        cat ~{str} >> out_file
    >>>

    output {
        File out_file = "out_file"
    }

    runtime {
        docker: docker_image_id
    }
}
