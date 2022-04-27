version 1.1

workflow sql {
    input {
        String docker_image_id
    }

    call DummyQuery {
        input:
        docker_image_id = docker_image_id
    }

    output {
    }
}

task DummyQuery {
    input {
        String docker_image_id
    }

    command <<<
        python3 /usr/local/bin/dummy.py
    >>>

    output {
    }

    runtime {
        docker: docker_image_id
    }
}