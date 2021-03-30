import os
import pytest


@pytest.fixture(scope="session")
def short_read_mngs_bench3_viral_outputs(util):
    """
    Run short-read-mngs/local_driver.wdl on the synthetic bench3 outputs, using the viral reference databases (~6GB
    total download)
    """
    return util.miniwdl_run(
        util.repo_dir() / "short-read-mngs/local_driver.wdl",
        "docker_image_id=" + os.environ.get("DOCKER_IMAGE_ID", "idseq-short-read-mngs"),
        "fastqs_0=",
        util.repo_dir()
        / "short-read-mngs/test/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R1.fastq.gz",
        "fastqs_1=",
        util.repo_dir()
        / "short-read-mngs/test/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R2.fastq.gz",
        "-i",
        util.repo_dir() / "short-read-mngs/test/local_test_viral.yml",
    )
