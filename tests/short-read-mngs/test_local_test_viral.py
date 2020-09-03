import os
import pytest
import subprocess


def test_bench3_viral(repo_dir, tmpdir):
    """
    Test the full short-read-mngs workflow on the bench3 synthetic sample, using smaller viral
    reference databases (~6GB total download)
    """
    if "DOCKER_IMAGE_ID" not in os.environ:
        pytest.skip("set env var DOCKER_IMAGE_ID for test_bench3_viral")
    subprocess.run(
        [
            "miniwdl",
            "run",
            os.path.join(repo_dir, "short-read-mngs/local_driver.wdl"),
            "docker_image_id=" + os.environ.get("DOCKER_IMAGE_ID"),
            "fastqs_0="
            + os.path.join(
                repo_dir,
                "tests/short-read-mngs/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R1.fastq.gz",
            ),
            "fastqs_1="
            + os.path.join(
                repo_dir,
                "tests/short-read-mngs/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R2.fastq.gz",
            ),
            "-i",
            os.path.join(repo_dir, "tests/short-read-mngs/local_test_viral.yml"),
            "--verbose",
            "--dir",
            tmpdir,
        ],
        check=True,
    )

    # TODO: load outputs and sanity check
