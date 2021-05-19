import os
import csv
import json


def test_RunIDSeqDedup_safe_csv(util, short_read_mngs_bench3_viral_outputs):
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"], "call-host_filter/call-RunIDSeqDedup"
        )
    )

    outp = util.miniwdl_run(
        util.repo_dir() / "short-read-mngs/host_filter.wdl",
        "--task",
        "RunIDSeqDedup",
        "-i",
        json.dumps(inputs),
    )

    dups = outp["outputs"]["RunIDSeqDedup.duplicate_clusters_csv"]

    # check we have an initial space to prevent CSV injection
    with open(dups) as f:
        for row in csv.reader(f):
            for elem in row:
                assert elem[0] == " ", f"cell does not have initial space '{elem}'"

    # check we can parse our CSV with skipinitialspace
    with open(dups) as f:
        for row in csv.reader(f, skipinitialspace=True):
            for elem in row:
                assert elem[0] != " ", f"initial space was not stripped from cell '{elem}'"
