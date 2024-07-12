import hashlib
import json
import logging
import os
import random
import re
import requests
import shutil
import time
from os import listdir
from multiprocessing import Pool
from subprocess import run
from typing import Dict, List, Optional
from urllib.parse import urlparse

from idseq_utils.diamond_scatter import blastx_join
from idseq_utils.minimap2_scatter import minimap2_merge

import boto3
from botocore.exceptions import ClientError
from botocore.config import Config

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
log = logging.getLogger(__name__)

MAX_CHUNKS_IN_FLIGHT = 30  # TODO: remove this constant, currently does nothing since we have at most 30 index chunks

# mitigation for TooManyRequestExceptions
config = Config(
    retries={
        "max_attempts": 20,
        "mode": "adaptive",
    }
)


_batch_client = boto3.client("batch", config=config)
# the retries are less necessary for S3 because rate limiting is more generous but we have hit the limit
#   and it is costly for this job to be re-run
_s3_client = boto3.client("s3", config=config)

try:
    account_id = boto3.client("sts").get_caller_identity()["Account"]
except ClientError:
    token = requests.put(
        "http://169.254.169.254/latest/api/token",
        headers={"X-aws-ec2-metadata-token-ttl-seconds": "21600"},
    ).text
    account_id = requests.get(
        "http://169.254.169.254/latest/dynamic/instance-identity/document",
        headers={"X-aws-ec2-metadata-token": token},
    ).json()["accountId"]


class BatchJobFailed(Exception):
    pass


def _bucket_and_key(s3_path: str):
    parsed_url = urlparse(s3_path, allow_fragments=False)
    return parsed_url.netloc, parsed_url.path.lstrip("/")


def _get_job_status(job_id, use_batch_api=False):
    if use_batch_api:
        jobs = _batch_client.describe_jobs(jobs=[job_id])["jobs"]
        if not jobs:
            log.debug(f"missing_job_description_from_api: {job_id}")
            return "SUBMITTED"
        return jobs[0]["status"]
    batch_job_desc_bucket = boto3.resource("s3").Bucket(
        f"aegea-batch-jobs-{account_id}"
    )
    key = f"job_descriptions/{job_id}"
    try:
        job_desc_object = batch_job_desc_bucket.Object(key)
        return json.loads(job_desc_object.get()["Body"].read())["status"]
    except ClientError as e:
        if e.response["Error"]["Code"] == "NoSuchKey":
            # Warn that the object is missing so any issue with the s3 mechanism can be identified
            log.debug(f"missing_job_description_object key: {key}")
            # Return submitted because a missing job status probably means it hasn't been added yet
            return "SUBMITTED"
        else:
            raise e


class BatchJobCache:
    """
    BatchJobCache saves job IDs so the coordinator can re-attach to running batch jobs when the coordinator fails

    The output should always be the same if the inputs are the same, however we also incorporate the batch_args
    into the cache because a retry on spot vs on demand will result in a different batch queue.
    """
    def __init__(self, bucket: str, prefix: str, inputs: Dict[str, str]):
        self.bucket = bucket
        self.prefix = prefix
        self.inputs = inputs

    def _key(self, batch_args: Dict) -> str:
        hash = hashlib.sha256()
        cache_dict = {"inputs": self.inputs, "batch_args": batch_args}
        hash.update(json.dumps(cache_dict, sort_keys=True).encode())
        return os.path.join(self.prefix, hash.hexdigest())

    def get(self, batch_args: Dict) -> Optional[str]:
        try:
            resp = _s3_client.get_object(Bucket=self.bucket, Key=self._key(batch_args))
            return resp["Body"].read().decode()
        except ClientError as e:
            if e.response["Error"]["Code"] == "NoSuchKey":
                return None
            else:
                raise e

    def put(self, batch_args: Dict, job_id: str):
        _s3_client.put_object(
            Bucket=self.bucket,
            Key=self._key(batch_args),
            Body=job_id.encode(),
            Tagging="AlignmentCoordination=True",
        )


def _run_batch_job(
    job_name: str,
    job_queue: str,
    job_definition: str,
    environment: Dict[str, str],
    retries: int,
    cache: BatchJobCache,
):
    submit_args = {
        "jobName": job_name,
        "jobQueue": job_queue,
        "jobDefinition": job_definition,
        "containerOverrides": {
            "environment": [{"name": k, "value": v} for k, v in environment.items()],
            "memory": 130816,
            "vcpus": 24,
        },
        "retryStrategy": {"attempts": retries},
    }

    def _log_status(status: str):
        level = logging.INFO if status != "FAILED" else logging.ERROR
        log.log(
            level,
            "batch_job_status " + json.dumps(
                {
                    "job_id": job_id,
                    "job_name": job_name,
                    "job_queue": job_queue,
                    "job_definition": job_definition,
                    "status": status,
                    "environment": environment,
                }
            ),
        )

    job_id = cache.get(submit_args)
    if job_id:
        log.info(f"reattach to batch job: {job_id}")
    else:
        response = _batch_client.submit_job(**submit_args)
        job_id = response["jobId"]
        cache.put(submit_args, job_id)
        _log_status("SUBMITTED")

    delay = 60 + random.randint(
        -60 // 2, 60 // 2
    )  # Add some noise to de-synchronize chunks
    status = "SUBMITTED"
    # the job this is monitoring has an timeout and the job this runs in has a timeout
    i = 0
    while True:
        try:
            status = _get_job_status(job_id, use_batch_api=(i > 0 and i % 30 == 0))
        except ClientError as e:
            # If we get throttled, randomly wait to de-synchronize the requests
            if e.response["Error"]["Code"] == "TooManyRequestsException":
                log.warn(f"describe_jobs_rate_limit_error for job_id: {job_id}")
                # Possibly implement a backoff here if throttling becomes an issue
            else:
                log.error(
                    f"unexpected_client_error_while_polling_job_status for job_id: {job_id}",
                )
                raise e

        if status == "SUCCEEDED":
            _log_status(status)
            return job_id
        if status == "FAILED":
            _log_status(status)
            raise BatchJobFailed("chunk alignment failed")
        time.sleep(delay)
        i += 1


def _run_chunk(
    input_dir: str,
    chunk_dir: str,
    aligner: str,
    aligner_args: str,
    aligner_wdl_version: str,
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

    def _job_queue(provisioning_model: str):
        return f"idseq-{deployment_environment}-{aligner}-{provisioning_model}-{priority_name}"

    priority_name = os.environ.get("PRIORITY_NAME", "normal")
    job_name = (
        f"idseq-{deployment_environment}-{aligner}-"
        f"project-{project_id}-sample-{sample_id}-part-{chunk_id}"
    )
    job_definition = f"idseq-swipe-{deployment_environment}-main"

    query_uris = [os.path.join(input_dir, os.path.basename(q)) for q in queries]
    inputs = {
        "query_0": query_uris[0],
        "extra_args": aligner_args,
        "db_chunk": db_chunk,
        "docker_image_id": f"{account_id}.dkr.ecr.us-west-2.amazonaws.com/{aligner}:{aligner_wdl_version}",
    }

    if len(query_uris) > 1:
        inputs["query_1"] = query_uris[1]

    wdl_input_uri = os.path.join(chunk_dir, f"{chunk_id}-input.json")
    input_bucket, input_key = _bucket_and_key(wdl_input_uri)

    wdl_output_uri = os.path.join(chunk_dir, f"{chunk_id}-output.json")

    wdl_workflow_uri = f"s3://idseq-workflows/{aligner}-{aligner_wdl_version}/{aligner}.wdl"

    cache_prefix_uri = os.path.join(chunk_dir, "batch_job_cache/")
    cache_bucket, cache_prefix = _bucket_and_key(cache_prefix_uri)
    cache = BatchJobCache(cache_bucket, cache_prefix, inputs)

    _s3_client.put_object(
        Bucket=input_bucket,
        Key=input_key,
        Body=json.dumps(inputs).encode(),
        ContentType="application/json",
        Tagging="AlignmentCoordination=True",
    )

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
            cache=cache,
        )
    except BatchJobFailed:
        _run_batch_job(
            job_name=job_name,
            job_queue=_job_queue("EC2"),
            job_definition=job_definition,
            environment=environment,
            retries=1,
            cache=cache,
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
    db_path: str,
    result_path: str,
    aligner: str,
    aligner_args: str,
    aligner_wdl_version: str,
    queries: List[str],
):
    db_bucket, db_prefix = _bucket_and_key(db_path)
    chunk_dir = os.path.join(input_dir, f"{aligner}-chunks")
    chunk_bucket, chunk_prefix = _bucket_and_key(chunk_dir)
    chunks = (
        [
            input_dir,
            chunk_dir,
            aligner,
            aligner_args,
            aligner_wdl_version,
            queries,
            chunk_id,
            f"s3://{db_bucket}/{db_chunk}",
        ]
        for chunk_id, db_chunk in enumerate(_db_chunks(db_bucket, db_prefix))
    )
    with Pool(MAX_CHUNKS_IN_FLIGHT) as p:
        p.starmap(_run_chunk, chunks)
    run(["s3parcp", "--recursive", chunk_dir, "chunks"], check=True)
    if os.path.exists(os.path.join("chunks", "cache")):
        shutil.rmtree(os.path.join("chunks", "cache"))
    if os.path.exists(os.path.join("chunks", "batch_job_cache")):
        shutil.rmtree(os.path.join("chunks", "batch_job_cache"))
    for fn in listdir("chunks"):
        if fn.endswith("json"):
            os.remove(os.path.join("chunks", fn))
            try:
                _s3_client.put_object_tagging(
                    Bucket=chunk_bucket,
                    Key=os.path.join(chunk_prefix, fn),
                    Tagging={"TagSet": [{"Key": "AlignmentCoordination", "Value": "True"}]},
                )
            except ClientError as e:
                log.error(f"failed to tag 's3://{chunk_bucket}/{os.path.join(chunk_prefix, fn)}'")
                raise e
    if aligner == "diamond":
        blastx_join("chunks", result_path, aligner_args, *queries)
    else:
        minimap2_merge("chunks", result_path, aligner_args, *queries)
