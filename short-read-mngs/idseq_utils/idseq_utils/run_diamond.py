import os
import re
from os.path import basename, join
from subprocess import run
import logging
from urllib.parse import urlparse
from multiprocessing import Pool

from idseq_utils.diamond_scatter import blastx_join
from idseq_utils.batch_run_helpers import _run_batch_job, _db_chunks

log = logging.getLogger(__name__)

MAX_CHUNKS_IN_FLIGHT = 10


def _run_chunk(
    chunk_id: int, input_dir: str, chunk_dir: str, db_chunk: str, *queries: str
):
    deployment_environment = os.environ["DEPLOYMENT_ENVIRONMENT"]
    priority_name = os.environ.get("PRIORITY_NAME", "normal")
    alignment_algorithm = "diamond"
    provisioning_model = "EC2"
    pattern = r"s3://.+/samples/([0-9]+)/([0-9]+)/"
    m = re.match(pattern, input_dir)
    if m:
        project_id, sample_id = m.group(1), m.group(2)
    else:
        project_id, sample_id = "0", "0"

    job_name = (f"idseq-{deployment_environment}-{alignment_algorithm}-"
                f"project-{project_id}-sample-{sample_id}-part-{chunk_id}")
    job_queue = f"idseq-{deployment_environment}-{alignment_algorithm}-{provisioning_model}-{priority_name}"
    job_definition = f"idseq-{deployment_environment}-{alignment_algorithm}"

    environment = [
        {
            "name": "DB_CHUNK",
            "value": db_chunk,
        },
        {
            "name": "OUTPUT_DIR",
            "value": chunk_dir,
        },
    ]

    for i, query in enumerate(queries):
        environment.append(
            {
                "name": f"QUERY_{i}",
                "value": join(input_dir, basename(query)),
            }
        )

    _run_batch_job(
        job_name=job_name,
        job_queue=job_queue,
        job_definition=job_definition,
        environment=environment,
        chunk_id=chunk_id,
        alignment_algorithm=alignment_algorithm,
        retries=2,
    )


def run_diamond(
    input_dir: str, chunk_dir: str, db_path: str, result_path: str, *queries: str
):
    parsed_url = urlparse(db_path, allow_fragments=False)
    bucket = parsed_url.netloc
    prefix = parsed_url.path.lstrip("/")
    chunks = (
        [chunk_id, input_dir, chunk_dir, f"s3://{bucket}/{db_chunk}", *queries]
        for chunk_id, db_chunk in enumerate(_db_chunks(bucket, prefix))
    )
    with Pool(MAX_CHUNKS_IN_FLIGHT) as p:
        p.starmap(_run_chunk, chunks)
    run(["s3parcp", "--recursive", chunk_dir, "chunks"], check=True)
    blastx_join("chunks", result_path, *queries)
