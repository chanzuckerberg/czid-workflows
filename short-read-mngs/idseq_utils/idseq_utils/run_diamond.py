import os
import re
from os.path import basename, join
from subprocess import run
import time
import random
import boto3
import json
import requests
import logging
from urllib.parse import urlparse
from multiprocessing import Pool
from botocore.exceptions import ClientError

from .diamond_scatter import blastx_join

log = logging.getLogger(__name__)

ALIGNMENT_ALGORITHM = "diamond"
MAX_CHUNKS_IN_FLIGHT = 10


def get_batch_job_desc_bucket():
    try:
        account_id = boto3.client("sts").get_caller_identity()["Account"]
    except ClientError:
        account_id = requests.get(
            "http://169.254.169.254/latest/dynamic/instance-identity/document"
        ).json()["accountId"]
    return f"aegea-batch-jobs-{account_id}"


class BatchJobFailed(Exception):
    pass


def _log_alignment_batch_job_status(
    job_id, job_queue, job_definition, chunk_id, status
):
    log.info(
        "alignment_batch_job_status",
        extra={
            "job_id": job_id,
            "chunk_id": chunk_id,
            "job_queue": job_queue,
            "job_definition": job_definition,
            "status": status,
            "alignment_algorithm": "diamond",
        },
    )


def _get_job_status(job_id):
    batch_job_desc_bucket = boto3.resource("s3").Bucket(get_batch_job_desc_bucket())
    key = f"job_descriptions/{job_id}"
    try:
        job_desc_object = batch_job_desc_bucket.Object(key)
        return json.loads(job_desc_object.get()["Body"].read())["status"]
    except ClientError as e:
        if e.response["Error"]["Code"] == "NoSuchKey":
            # Warn that the object is missing so any issue with the s3 mechanism can be identified
            log.debug("missing_job_description_ojbect", extra={key: key})
            # Return submitted because a missing job status probably means it hasn't been added yet
            return "SUBMITTED"
        else:
            raise e


def _run_batch_job(job_name, job_queue, job_definition, environment, chunk_id, retries):
    client = boto3.client("batch")
    response = client.submit_job(
        jobName=job_name,
        jobQueue=job_queue,
        jobDefinition=job_definition,
        containerOverrides={
            "environment": environment,
        },
        retryStrategy={"attempts": retries},
    )
    job_id = response["jobId"]
    _log_alignment_batch_job_status(
        job_id, job_queue, job_definition, chunk_id, "SUBMITTED"
    )

    delay = 60 + random.randint(
        -60 // 2, 60 // 2
    )  # Add some noise to de-synchronize chunks
    status = "SUBMITTED"
    # the job this is monitoring has an timeout and the job this runs in has a timeout
    while True:
        try:
            status = _get_job_status(job_id)
        except ClientError as e:
            # If we get throttled, randomly wait to de-synchronize the requests
            if e.response["Error"]["Code"] == "TooManyRequestsException":
                log.warn("describe_jobs_rate_limit_error", extra={"job_id": job_id})
                # Possibly implement a backoff here if throttling becomes an issue
            else:
                log.error(
                    "unexpected_client_error_while_polling_job_status",
                    extra={"job_id": job_id},
                )
                raise e

        if status == "SUCCEEDED":
            _log_alignment_batch_job_status(
                job_id, job_queue, job_definition, chunk_id, status
            )
            return job_id
        if status == "FAILED":
            log.error(
                "alignment_batch_job_failed",
                extra={
                    "job_id": job_id,
                    "chunk_id": chunk_id,
                    "alignment_algorithm": ALIGNMENT_ALGORITHM,
                },
            )
            _log_alignment_batch_job_status(
                job_id, job_queue, job_definition, chunk_id, status
            )
            raise BatchJobFailed("chunk alignment failed")
        time.sleep(delay)


def _db_chunks(bucket: str, prefix):
    s3_client = boto3.client("s3")
    paginator = s3_client.get_paginator("list_objects_v2")
    log.debug("db chunks")

    result = paginator.paginate(
        Bucket=bucket,
        Prefix=prefix,
    )

    for page in result:
        for obj in page["Contents"]:
            yield obj["Key"]


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
