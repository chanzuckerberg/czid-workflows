import os
import csv
import json
from tempfile import NamedTemporaryFile


def test_RunIDSeqDedup_safe_csv(util, short_read_mngs_bench3_viral_outputs):
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"], "call-host_filter/call-RunIDSeqDedup"
        )
    )

    with NamedTemporaryFile(prefix=os.path.dirname(__file__), mode="w") as input_file:
        quote_count = 10
        special_char_rows = 0
        for line in open(inputs["priceseq_fa"][0]):
            if line[0] == ">" or line[0] == "@":
                if special_char_rows < quote_count:
                    input_file.write(f"{line[0]}={line[1:]}")
                    special_char_rows += 1
                else:
                    input_file.write(line)
            else:
                input_file.write(line)

        input_file.seek(0)
        assert special_char_rows == quote_count

        inputs["priceseq_fa"] = [input_file.name]

        outp = util.miniwdl_run(
            util.repo_dir() / "short-read-mngs/host_filter.wdl",
            "--task",
            "RunIDSeqDedup",
            "-i",
            json.dumps(inputs),
        )

        dups = outp["outputs"]["RunIDSeqDedup.duplicate_clusters_csv"]

        found_quotes = 0
        # check we have an quotes before special characters space to prevent CSV injection
        with open(dups) as f:
            for row in csv.reader(f):
                for i, elem in enumerate(row):
                    if elem[0] == "'":
                        if i == 1:
                            found_quotes += 1
                        continue
                    assert elem[0].isalnum(), f"cell starts with a special character '{elem}'"
        assert found_quotes == quote_count
