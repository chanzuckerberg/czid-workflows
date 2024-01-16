version 1.1

workflow bulk_download {
    input {
        String action
        Array[File] files
        String docker_image_id
    }

    if (action == "concatenate") {
        call concatenate { 
            input:
                files = files
                docker_image_id = host_filtering_docker_image_id,
        }
    }

    if (action == "group") {
        call group { 
            input:
                files = files
                docker_image_id = host_filtering_docker_image_id,
        }
    }

    output {
        File? file = select_first([ concatenate.file, group.file ])
    }
}

task concatenate {
    input {
        String docker_image_id
        Array[File] files
    }
    command <<<
        set -euxo pipefail
        cat ~{sep=" " files} > concatenated.txt
    >>>
    output {
        File file = "concatenated.txt"
    }
    runtime {
        docker: docker_image_id
    }
}

task group {
    input {
        String docker_image_id
        Array[File] files
    }
    command <<<
        set -euxo pipefail
        zip ~{sep=" " files} > group.zip
    >>>
    output {
        File file = "group.zip"
    }
    runtime {
        docker: docker_image_id
    }
}
