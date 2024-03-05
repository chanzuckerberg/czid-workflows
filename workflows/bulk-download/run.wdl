version 1.1

workflow bulk_download {
    input {
        String action
        Array[File] files
        String docker_image_id = "czid-bulk-download"
    }

    if (action == "concatenate") {
        call concatenate { 
            input:
                files = files,
                docker_image_id = docker_image_id
        }
    }

    if (action == "zip") {
        call zip { 
            input:
                files = files,
                docker_image_id = docker_image_id
        }
    }

    output {
        File? file = select_first([ concatenate.file, zip.file ])
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

task zip {
    input {
        String docker_image_id
        Array[File] files
    }
    command <<<
        set -euxo pipefail

        # Don't store full path of original files in the .zip file
        zip --junk-paths result.zip ~{sep=" " files}
    >>>
    output {
        File file = "result.zip"
    }
    runtime {
        docker: docker_image_id
    }
}
