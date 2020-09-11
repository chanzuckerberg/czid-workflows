import os
import json


def test_RunValidateInput_invalid(
    repo_dir, short_read_mngs_bench3_viral_outputs, miniwdl_inputs_outputs, miniwdl_run, RunFailed_stderr_msg
):
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = miniwdl_inputs_outputs(
        os.path.join(short_read_mngs_bench3_viral_outputs["dir"], "call-host_filter/call-RunValidateInput")
    )
    # override fastqs to invalid test article
    inputs["fastqs"] = [os.path.join(os.path.dirname(__file__), "test_RunValidateInput_invalid.fastq")]

    # run the task with the manipulated inputs, expecting an error exit status
    outp = miniwdl_run(
        os.path.join(repo_dir, "short-read-mngs/host_filter.wdl"),
        "--task",
        "RunValidateInput",
        "-i",
        json.dumps(inputs),
        returncode=1,
    )

    # verify the error is as expected
    err = RunFailed_stderr_msg(outp)
    assert err["wdl_error_message"]
    assert err["error"] == "InvalidFileFormatError"
