import json
import logging
import os
import random
import re
import requests
import time
from multiprocessing import Pool
from subprocess import run
from typing import Dict, List, Literal
from urllib.parse import urlparse

from idseq_utils.diamond_scatter import blastx_join
from idseq_utils.minimap2_scatter import minimap2_merge

import boto3
from botocore.exceptions import ClientError

log = logging.getLogger(__name__)

MAX_CHUNKS_IN_FLIGHT = 10


_batch_client = boto3.client("batch")
_s3_client = boto3.client("boto3")


class BatchJobFailed(Exception):
    pass


def _bucket_and_key(s3_path: str):
    parsed_url = urlparse(s3_path, allow_fragments=False)
    return parsed_url.netloc, parsed_url.path.lstrip("/")


def _get_batch_job_desc_bucket():
    try:
        account_id = boto3.client("sts").get_caller_identity()["Account"]
    except ClientError:
        account_id = requests.get(
            "http://169.254.169.254/latest/dynamic/instance-identity/document"
        ).json()["accountId"]
    return f"aegea-batch-jobs-{account_id}"


def _get_job_status(job_id):
    batch_job_desc_bucket = boto3.resource("s3").Bucket(_get_batch_job_desc_bucket())
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


def _run_batch_job(
    job_name: str,
    job_queue: str,
    job_definition: str,
    environment: Dict[str, str],
    retries: int,
):
    response = _batch_client.submit_job(
        jobName=job_name,
        jobQueue=job_queue,
        jobDefinition=job_definition,
        containerOverrides={
            "environment": [{
                "name": k,
                "value": v,
            } for k, v in environment.items()],
        },
        retryStrategy={"attempts": retries},
    )
    job_id = response["jobId"]

    def _log_status(status: str):
        level = logging.INFO if status != "FAILED" else logging.ERROR
        log.log(
            level,
            "batch_job_status",
            extra={
                "job_id": job_id,
                "job_name": job_name,
                "job_queue": job_queue,
                "job_definition": job_definition,
                "status": status,
            },
        )

    _log_status("SUBMITTED")

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
            _log_status(status)
            return job_id
        if status == "FAILED":
            _log_status(status)
            raise BatchJobFailed("chunk alignment failed")
        time.sleep(delay)


def _run_chunk(
    input_dir: str,
    chunk_dir: str,
    result_path: str,
    alignment_algorithm: Literal["diamond", "minimap2"],
    aligner_args: str,
    queries: List[str],
    chunk_id: int,
    db_chunk: str,
):
    deployment_environment = os.environ["DEPLOYMENT_ENVIRONMENT"]
    pattern = r"s3://.+/samples/([0-9]+)/([0-9]+)/"
    m = re.match(pattern, input_dir)
    if m:
        project_id, sample_id = m.group(1), m.group(2)
    else:
        project_id, sample_id = "0", "0"

    def _job_queue(provisioning_model: Literal["SPOT", "EC2"]):
        return f"idseq-{deployment_environment}-{alignment_algorithm}-{provisioning_model}-{priority_name}"

    priority_name = os.environ.get("PRIORITY_NAME", "normal")
    job_name = (f"idseq-{deployment_environment}-{alignment_algorithm}-"
                f"project-{project_id}-sample-{sample_id}-part-{chunk_id}")
    job_definition = f"idseq-swipe-{deployment_environment}-main"

    query_uris = [os.path.join(input_dir, os.path.basename(q)) for q in queries]
    inputs = {
        "query_0": query_uris[0],
        "extra_args": aligner_args,
        "db_chunk": db_chunk,
        "docker_image_id": "",  # TODO hardcode this based on alignment_algorithm
    }

    if len(query_uris) > 1:
        inputs["query_1"] = query_uris[1]

    input_bucket, _ = _bucket_and_key(chunk_dir)
    wdl_input_key = os.path.join(chunk_dir, f"{chunk_id}-input.json")
    wdl_output_key = os.path.join(chunk_dir, f"{chunk_id}-output.json")

    _s3_client.put_object(
        Bucket=input_bucket,
        Key=wdl_input_key,
        Body=json.dumps(inputs).encode(),
        ContentType="application/json",
    )

    wdl_workflow_uri = ""  # TODO harcode this based on alignment_algorithm
    wdl_input_uri = f"s3://{input_bucket}/{wdl_input_key}"
    wdl_output_uri = f"s3://{input_bucket}/{wdl_output_key}"

    environment = {
        "WDL_WORKFLOW_URI": wdl_workflow_uri,
        "WDL_INPUT_URI": wdl_input_uri,
        "WDL_OUTPUT_URI": wdl_output_uri,
        "SFN_EXECUTION_ID": os.getenv("SFN_EXECUTION_ID", "ALIGNMENT"),
        "SFN_CURRENT_STATE": os.getenv("SFN_CURRENT_STATE", "ALIGNMENT"),
    }

    try:
        _run_batch_job(
            job_name=job_name,
            job_queue=_job_queue("SPOT"),
            job_definition=job_definition,
            environment=environment,
            retries=2,
        )
    except BatchJobFailed:
        _run_batch_job(
            job_name=job_name,
            job_queue=_job_queue("EC2"),
            job_definition=job_definition,
            environment=environment,
            retries=1,
        )


def _db_chunks(bucket: str, prefix):
    s3_client = boto3.client("s3")
    paginator = s3_client.get_paginator("list_objects_v2")
    log.debug("db chunks")

    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page["Contents"]:
            yield obj["Key"]


def run_alignment(
    input_dir: str,
    chunk_dir: str,
    db_path: str,
    result_path: str,
    alignment_algorithm: Literal["diamond", "minimap2"],
    aligner_args: str,
    queries: List[str],
):
    bucket, prefix = _bucket_and_key(db_path)
    chunks = (
        [chunk_id, input_dir, chunk_dir, f"s3://{bucket}/{db_chunk}", aligner_args, *queries]
        for chunk_id, db_chunk in enumerate(_db_chunks(bucket, prefix))
    )
    with Pool(MAX_CHUNKS_IN_FLIGHT) as p:
        p.starmap(_run_chunk, chunks)
    run(["s3parcp", "--recursive", chunk_dir, "chunks"], check=True)
    if alignment_algorithm == "diamond":
        blastx_join("chunks", result_path, aligner_args, *queries)
    else:
        minimap2_merge("chunks", result_path, aligner_args, *queries)
