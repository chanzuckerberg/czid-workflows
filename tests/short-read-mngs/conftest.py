import os
import pytest
import subprocess
import json
import atexit


@pytest.fixture(scope="session")
def miniwdl_run_short_read_mngs(repo_dir, miniwdl_run):
    """
    miniwdl_run_short_read_mngs(*inputs, returncode) runs short-read-mngs/local_driver.wdl with the given input
    arguments
    """

    def kappa(*args, **kwargs):
        return miniwdl_run(
            os.path.join(repo_dir, "short-read-mngs/local_driver.wdl"),
            "docker_image_id=" + os.environ.get("DOCKER_IMAGE_ID", "idseq-short-read-mngs"),
            *args,
            **kwargs,
        )

    return kappa


@pytest.fixture(scope="session")
def short_read_mngs_bench3_viral_outputs(repo_dir, miniwdl_run_short_read_mngs):
    """
    Run short-read-mngs/local_driver.wdl on the synthetic bench3 outputs, using the viral reference databases (~6GB
    total download)
    """
    return miniwdl_run_short_read_mngs(
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
    )
