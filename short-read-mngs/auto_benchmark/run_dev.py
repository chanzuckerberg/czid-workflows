#!/usr/bin/env python3
"""
Run benchmark.yml samples by submitting them to the idseq-dev SFN-WDL backend. The invoking session
must have (i) AWS_PROFILE set to an appropriate role to access idseq-dev, and (ii) an SSH key able
to clone chanzuckerberg/idseq.

The SFN-WDL runs will use the full databases and the latest -released- version of the WDL code, not
necessarily the current git checkout.

After the runs complete, one can harvest the results e.g.
    harvest.py idseq_bench_3=s3://idseq-samples-development/auto_benchmark/YYYYMMDD_HHmmSS_default_latest/results/short-read-mngs-V
"""  # noqa
import sys
import os
import argparse
import subprocess
import concurrent.futures
import threading
import time
import json
import requests
from datetime import datetime
from _util import load_benchmarks_yml

BENCHMARKS = load_benchmarks_yml()
BUCKET = "idseq-samples-development"
KEY_PREFIX = "auto_benchmark"


def main():
    parser = argparse.ArgumentParser(
        sys.argv[0], description="run benchmark samples on idseq-dev SFN-WDL"
    )
    parser.add_argument(
        "samples",
        metavar="SAMPLE",
        type=str,
        nargs="+",
        help="any of: " + ", ".join(BENCHMARKS["samples"].keys()),
    )
    parser.add_argument("--idseq", required=True, help="path to idseq monorepo")
    parser.add_argument(
        "--workflow-version",
        metavar="X.Y.Z",
        type=str,
        default=None,
        help="short-read-mngs version tag",
    )
    parser.add_argument(
        "--settings",
        metavar="ID",
        type=str,
        default="default",
        help="settings presets; any of: " + ",".join(BENCHMARKS["settings"].keys()),
    )

    args = parser.parse_args(sys.argv[1:])

    if args.workflow_version is None:
        github_api = "https://api.github.com"
        github_repo_api = f"{github_api}/repos/chanzuckerberg/idseq-workflows"
        github_refs_api = f"{github_repo_api}/git/matching-refs/tags/short-read-mngs"
        github_refs = sorted(
            os.path.basename(ref["url"]).split("-v", 1)[1]
            for ref in requests.get(github_refs_api).json()
        )
        args.workflow_version = github_refs[-1]

    args.workflow_version = args.workflow_version.lstrip("v")
    run_samples(**vars(args))


def run_samples(idseq, samples, workflow_version, settings):
    samples = set(samples)
    for sample_i in samples:
        assert sample_i in BENCHMARKS["samples"], f"unknown sample {sample_i}"

    failures = []
    results = []

    # formulate S3 directory for these benchmarking runs
    key_prefix = (
        f"{KEY_PREFIX}/{datetime.today().strftime('%Y%m%d_%H%M%S')}_{settings}_{workflow_version}"
    )

    # parallelize run_sample on a thread pool
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
        futures = {
            executor.submit(
                run_sample,
                idseq,
                workflow_version,
                settings,
                key_prefix,
                sample_i,
            ): sample_i
            for sample_i in samples
        }
        for future in concurrent.futures.as_completed(futures):
            exn = future.exception()
            if exn:
                failures.append((futures[future], exn))
            else:
                results.append(future.result())

    print("\n".join(f"{sample}={s3path}" for sample, s3path in results))

    if failures:
        for (failed_sample, exn) in failures:
            print(
                f"\nsample {failed_sample} failed; see logs under s3://{BUCKET}/{key_prefix}/{failed_sample}/\n",
                file=sys.stderr,
            )
            print(exn, file=sys.stderr)
        raise failures[0][1]


_timestamp_lock = threading.Lock()


def run_sample(idseq_repo, workflow_version, settings, key_prefix, sample):
    local_input = {
        **BENCHMARKS["settings"][settings],
        **BENCHMARKS["databases"]["full"],
        **BENCHMARKS["samples"][sample]["inputs"],
    }

    # copy inputs into the S3 directory where run_sfn.py will look for them
    sample_inputs = BENCHMARKS["samples"][sample]["inputs"]
    alt_inputs = next((k for k in sample_inputs if k.startswith("s3_")), None) is not None
    for k, v in sample_inputs.items():
        if (not alt_inputs and k.startswith("fastqs_")) or (
            alt_inputs and k.startswith("s3_fastqs_")
        ):
            cmd = ["aws", "s3", "cp", v, f"s3://{BUCKET}/{key_prefix}/{sample}/fastqs/"]
            print(" ".join(cmd), file=sys.stderr)
            subprocess.run(cmd, check=True, stdout=sys.stderr.buffer)

    # run_sfn.py
    cmd = [
        "/bin/bash",
        "-c",
        "source workflows/environment.dev"
        f" && scripts/run_sfn.py"
        f" --sample-dir s3://{BUCKET}/{key_prefix}/{sample}"
        f" --host-genome {local_input['host_filter.host_genome']}"
        f" --max-input-fragments {local_input['host_filter.max_input_fragments']}"
        f" --max-subsample-fragments {local_input['host_filter.max_subsample_fragments']}"
        f" --adapter-fasta {local_input['host_filter.adapter_fasta']}"
        f" --workflow-version {workflow_version}"
        f" --sfn-input '{json.dumps(local_input)}'"
    ]
    print(cmd, file=sys.stderr)
    with _timestamp_lock:
        time.sleep(1.1)
        subprocess.run(cmd, cwd=idseq_repo, check=True, stdout=sys.stderr.buffer)

    workflow_major_version = workflow_version.split(".")[0]
    return (
        sample,
        f"s3://{BUCKET}/{key_prefix}/{sample}/results/short-read-mngs-{workflow_major_version}/",
    )


if __name__ == "__main__":
    main()
