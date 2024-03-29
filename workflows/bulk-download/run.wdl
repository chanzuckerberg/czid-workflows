version 1.1

struct BulkDownloads { 
    File file_path
    String output_name
}

workflow bulk_download {
    input {
        String action
        Array[BulkDownloads] files
        String docker_image_id = "czid-bulk-download"
    }

    scatter (file in files) {
        call rename {
            input: 
                input_file = file.file_path,
                output_name = file.output_name,
                docker_image_id = docker_image_id
        }
    }

    if (action == "concatenate") {
        call concatenate { 
            input:
                files = rename.file,
                docker_image_id = docker_image_id
        }
    }

    if (action == "zip") {
        call zip { 
            input:
                files = rename.file,
                docker_image_id = docker_image_id
        }
    }

    output {
        File? file = select_first([ concatenate.file, zip.file ])
    }
}

task rename { 
    input {
        String docker_image_id
        File input_file
        String output_name
    }
    command <<<
        set -euxo pipefail
        cp "~{input_file}" "~{output_name}"
    >>>
    output {
        File file = "~{output_name}"
    }
    runtime { 
        docker: docker_image_id
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
