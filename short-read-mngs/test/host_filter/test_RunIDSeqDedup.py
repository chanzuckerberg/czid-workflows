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

    with open(dups) as f:
        for row in csv.reader(f):
            for elem in row:
                # check for quotes that prevent csv injection
                assert elem[0] == "'" and elem[-1] == "'", f"line not surrounded with ': {elem}"
