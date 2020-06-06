version 1.0
# positive test cases for ../main/host_filter.wdl

import "../../main/host_filter.wdl"
import "../util/assert.wdl"

workflow test_host_filter {
    input {
        # fixtures supplied by run_tests.t
        String docker_image_id
        File test_fastq_gz
        File run_validate_input_multiline_fastq

        # control-flow flags to exercise different code paths
        Boolean invalid_fastq = false
    }

    # run the workflow; here we're just calling it once but we could call it in various
    # configurations, perhaps inside a scatter
    File test_fastq = if !invalid_fastq then test_fastq_gz else run_validate_input_multiline_fastq
    call host_filter.RunValidateInput {
        input:
        docker_image_id     = docker_image_id,
        aws_region          = "local",
        deployment_env      = "local",
        dag_branch          = "",
        s3_wd_uri           = "s3://idseq-workflows-local/test_host_filter/",
        max_input_fragments = 75000000,
        file_ext            = "fastq",
        fastqs              = [test_fastq]
    }

    # check output validity using an "assert" WDL task (more such tasks could be added in
    # util/assert.wdl)
    Int validate_count = read_json(RunValidateInput.validate_input_summary_json)["50-500"]
    call assert.IsTrue {
        input:
        cond = validate_count == 250,
        msg = 'RunValidateInput.validate_input_summary_json["50-500"] != 250'
    }

    # error messages are checked in run_tests.t
}
