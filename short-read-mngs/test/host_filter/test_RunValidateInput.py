import os
import re
import gzip
import json
import tempfile


def test_RunValidateInput_invalid(util, short_read_mngs_bench3_viral_outputs):
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"], "call-host_filter/call-RunValidateInput"
        )
    )
    # override fastqs to invalid test article
    inputs["fastqs"] = [
        os.path.join(os.path.dirname(__file__), "test_RunValidateInput_invalid.fastq")
    ]

    # run the task with the manipulated inputs, expecting an error exit status
    outp = util.miniwdl_run(
        util.repo_dir() / "short-read-mngs/host_filter.wdl",
        "--task",
        "RunValidateInput",
        "-i",
        json.dumps(inputs),
        returncode=1,
    )

    # verify the error is as expected
    err = util.wdl_error_message(outp)
    assert err["wdl_error_message"]
    assert err["error"] == "InvalidFileFormatError"


def test_RunValidateInput_strip_bad_csv_characters(util, short_read_mngs_bench3_viral_outputs):
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"], "call-host_filter/call-RunValidateInput"
        )
    )

    with tempfile.NamedTemporaryFile('w') as input_fastq, gzip.open(inputs["fastqs"][0], 'r') as good_fastq:
        for i, line in enumerate(good_fastq):
            if i == 0:
                clean_line = line
                dirty_line = clean_line.strip() + "=+-@|\n"
                input_fastq.write(dirty_line)
            else:
                input_fastq.write(line)

        # override fastqs to test
        inputs["fastqs"] = [input_fastq.name]

        # run the task with the manipulated inputs, expecting an error exit status
        outp = util.miniwdl_run(
            util.repo_dir() / "short-read-mngs/host_filter.wdl",
            "--task",
            "RunValidateInput",
            "-i",
            json.dumps(inputs),
        )

        bad = re.compile('[=+-@|]')
        with open(outp["outputs"]["RunValidateInput.valid_input1_fastq"]) as o:
            for line in o:
                assert not bad.match(line[1:]), "found bad csv character in line"
