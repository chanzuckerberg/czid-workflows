version 1.1

struct BulkDownloads { 
    File file_path
    String output_name
}

workflow bulk_download {
    input {
        String action
        Array[BulkDownloads] files
        String concatenated_output_name = "concatenated.txt"
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
                concatenated_output_name = concatenated_output_name,
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
        String concatenated_output_name = "concatenated.txt"
        Array[File] files
    }
    command <<<
        set -euxo pipefail
        cat ~{sep=" " files} > "~{concatenated_output_name}"
    >>>
    output {
        File file = "~{concatenated_output_name}"
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
        if [[ "~{select_first(files)}" == *.zip ]]; then
            mkdir zip_folders
            for f in ~{sep=" " files}; do unzip "$f" -d zip_folders/$(basename "${f%.zip}"); done
            cd zip_folders
            zip ../result.zip *
        else
            zip --junk-paths result.zip ~{sep=" " files}
        fi 
    >>>
    output {
        File file = "result.zip"
    }
    runtime {
        docker: docker_image_id
    }
}
