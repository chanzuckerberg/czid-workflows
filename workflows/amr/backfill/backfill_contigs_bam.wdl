version 1.1

workflow amr {
    input {
        File execution_array  # json file with workflow information
        String docker_image_id
        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call backfillContigsBam { 
        input: 
        execution_array = execution_array
        docker_image_id = docker_image_id
    }
}

task backfillContigsBam {
    input {
        File execution_array
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        
    >>>

    output {
        File execution_log = "backfill_log.tsv"
    }

    runtime {
        docker: docker_image_id
    }
}
