import os
import json
import csv
import tempfile


def test_RunAlignmentBlacklist(util, short_read_mngs_bench3_viral_outputs):
    task_name = "RunAlignment_gsnap_out"
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"],
            "call-non_host_alignment",
            f"call-{task_name}",
        )
    )

    outp = util.miniwdl_run(
        util.repo_dir() / "short-read-mngs/non_host_alignment.wdl",
        "--task",
        task_name,
        "-i",
        json.dumps(inputs),
    )

    with open(os.path.join(outp["dir"], outp["outputs"][f"{task_name}.gsnap_hitsummary_tab"])) as f:
        taxids = set(row[2] for row in csv.reader(f, delimiter="\t"))

    assert "37124" in taxids, "taxid should be in hitsummary unless filtered out"
    assert "1273712" in taxids, "taxid should be in hitsummary unless filtered out"

    with tempfile.NamedTemporaryFile(prefix=os.path.dirname(__file__), mode="w") as blacklist_file:
        blacklist_file.writelines(["37124\n", "1273712\n"])
        blacklist_file.seek(0)
        blacklist_file.writelines
        inputs["taxon_blacklist"] = blacklist_file.name

        outp = util.miniwdl_run(
            util.repo_dir() / "short-read-mngs/non_host_alignment.wdl",
            "--task",
            task_name,
            "-i",
            json.dumps(inputs),
        )

        hitsummary = os.path.join(outp["dir"], outp["outputs"][f"{task_name}.gsnap_hitsummary_tab"])
        deduped = os.path.join(outp["dir"], outp["outputs"][f"{task_name}.gsnap_deduped_m8"])

        with open(hitsummary) as f:
            taxids = set(row[2] for row in csv.reader(f, delimiter="\t"))

        assert "37124" not in taxids, "taxid should be filtered out"
        assert "1273712" not in taxids, "taxid should be filtered out"

        with open(hitsummary) as hf, open(deduped) as df:
            rows = zip(csv.reader(hf, delimiter="\t"), csv.reader(df, delimiter="\t"))
            assert all(
                hrow[0] == drow[0] for hrow, drow in rows
            ), "hitsummary and deduped output should be aligned"


def test_RunAlignmentDeuterostomeFilter(util, short_read_mngs_bench3_viral_outputs):
    task_name = "RunAlignment_gsnap_out"
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"],
            "call-non_host_alignment",
            f"call-{task_name}",
        )
    )

    outp = util.miniwdl_run(
        util.repo_dir() / "short-read-mngs/non_host_alignment.wdl",
        "--task",
        task_name,
        "-i",
        json.dumps(inputs),
    )

    with open(os.path.join(outp["dir"], outp["outputs"][f"{task_name}.gsnap_hitsummary_tab"])) as f:
        taxids = set(row[2] for row in csv.reader(f, delimiter="\t"))

    assert "37124" in taxids, "taxid should be in hitsummary unless filtered out"
    assert "1273712" in taxids, "taxid should be in hitsummary unless filtered out"

    with tempfile.NamedTemporaryFile(
        prefix=os.path.dirname(__file__), mode="w"
    ) as deuterostome_file:
        deuterostome_file.writelines(["37124\n", "1273712\n"])
        deuterostome_file.seek(0)
        inputs["deuterostome_db"] = deuterostome_file.name
        inputs["use_deuterostome_filter"] = True

        outp = util.miniwdl_run(
            util.repo_dir() / "short-read-mngs/non_host_alignment.wdl",
            "--task",
            task_name,
            "-i",
            json.dumps(inputs),
        )

        hitsummary = os.path.join(outp["dir"], outp["outputs"][f"{task_name}.gsnap_hitsummary_tab"])
        deduped = os.path.join(outp["dir"], outp["outputs"][f"{task_name}.gsnap_deduped_m8"])

        with open(hitsummary) as f:
            taxids = set(row[2] for row in csv.reader(f, delimiter="\t"))

        assert "37124" not in taxids, "taxid should be filtered out"
        assert "1273712" not in taxids, "taxid should be filtered out"

        with open(hitsummary) as hf, open(deduped) as df:
            rows = zip(csv.reader(hf, delimiter="\t"), csv.reader(df, delimiter="\t"))
            assert all(
                hrow[0] == drow[0] for hrow, drow in rows
            ), "hitsummary and deduped output should be aligned"
