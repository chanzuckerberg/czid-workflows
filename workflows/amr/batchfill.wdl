version 1.1

workflow amr {
    input {
        File backfill_data
        String docker_image_id
        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call backfillContigsBam {
        input: 
        backfill_data = backfill_data,
        docker_image_id = docker_image_id
    }
}

task backfillContigsBam {
    input {
        File backfill_data
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        python3 /scripts/backfill_contigs_bam.py "~{backfill_data}"
    >>>

    output {
        File output_log = "backfill_log.tsv"
    }

    runtime {
        docker: docker_image_id
    }
}
