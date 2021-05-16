import multiprocessing
import os
import random
import shutil
import threading
import time
import traceback
import json
import re
import tempfile
from subprocess import run, PIPE
from urllib.parse import urlparse
from botocore.exceptions import ClientError

import boto3
import requests

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.exceptions import InsufficientReadsError
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.log as log
import idseq_dag.util.m8 as m8

from idseq_dag.util.s3 import fetch_reference
from idseq_dag.util.trace_lock import TraceLock

from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3, DEFAULT_WHITELIST_S3
from idseq_dag.util.m8 import NT_MIN_ALIGNMENT_LEN

MAX_CHUNKS_IN_FLIGHT = 16
CHUNK_MAX_ATTEMPTS = 3
GSNAP_CHUNK_SIZE = 240000
RAPSEARCH_CHUNK_SIZE = 320000


def get_batch_job_desc_bucket():
    try:
        account_id = boto3.client("sts").get_caller_identity()["Account"]
    except ClientError:
        account_id = requests.get("http://169.254.169.254/latest/dynamic/instance-identity/document").json()["accountId"]
    return f"aegea-batch-jobs-{account_id}"

def download_from_s3(session, src, dest):
    """
    We are transitioning away from s3 downloads of files within steps.
    This is only to be used for getting the state of batch jobs from the
    s3 cache. These files are small and dynamically generated based on jobs
    so this fetching is very different from downloading input files.
    """
    try:
        url = urlparse(src)
        bucket, key = url.netloc, url.path
        session.client("s3").download_file(bucket, key[1:], dest)
        return True
    except ClientError as e:
        if e.response["Error"]["Code"] == "404":
            return False
        raise e

class BatchJobFailed(Exception):
    pass

class PipelineStepRunAlignment(PipelineStep):
    """ Runs gsnap/rapsearch2 remotely.

    For GSNAP:
    ```
    gsnapl
    -A m8
    --batch=0
    --use-shared-memory=0
    --gmap-mode=none
    --npaths=100
    --ordered
    -t 48
    --max-mismatches=40
    -D {remote_index_dir}
    -d nt_k16
    {remote_input_files} > {multihit_remote_outfile}
    ```

    GSNAP documentation is available [here](http://research-pub.gene.com/gmap/).
    -t (threads): For example, r5d.24xlarge machines have 96 vCPUs. These are the machines used by the alignment batch compute environment.
    Each batch job is allotted 48 vcpus. Use 48 threads and each instance will be able to concurrently process 2 chunks

    For Rapsearch2:
    ```
    rapsearch
    -d {remote_index_dir}/nr_rapsearch
    -e -6
    -l 10
    -a T
    -b 0
    -v 50
    -z 24
    -q {remote_input_files}
    -o {multihit_remote_outfile}
    ```

    Rapsearch2 documentation is available [here](http://omics.informatics.indiana.edu/mg/RAPSearch2/).
    """

    @staticmethod
    def _alignment_algorithm_inputs(host_filter_outputs):
        # Destructuring-bind for host filter outputs.
        # The host filter outputs either 1 file for unpaired reads, let's call it R1;
        # or, 3 files for paired-end reads, which can be succinctly described
        # as [R1, R2, R1/1 + R2/2].  All are usually named gsnap_filter_*.
        try:
            # Unpaired?
            gsnap_filter_1, = host_filter_outputs
            return {
                "gsnap": [gsnap_filter_1],
                "rapsearch2": [gsnap_filter_1]
            }
        except:
            # Paired!
            gsnap_filter_1, gsnap_filter_2, gsnap_filter_merged = host_filter_outputs
            return {
                "gsnap": [gsnap_filter_1, gsnap_filter_2],
                "rapsearch2": [gsnap_filter_merged]
            }

    def validate_input_files(self):
        # first two files are gsnap_filter_1.fa and gsnap_filter_2.fa
        if not count.files_have_min_reads(self.input_files_local[0][:-1], 1):
            raise InsufficientReadsError("Insufficient reads")

    def __init__(self, *args, **kwrds):
        PipelineStep.__init__(self, *args, **kwrds)
        self.alignment_algorithm = self.additional_attributes.get("alignment_algorithm")
        assert self.alignment_algorithm in ("gsnap", "rapsearch2")
        self.chunks_in_flight_semaphore = threading.Semaphore(MAX_CHUNKS_IN_FLIGHT)
        self.chunks_result_dir_local = os.path.join(self.output_dir_local, "chunks")
        self.chunks_result_dir_s3 = os.path.join(self.output_dir_s3, "chunks")
        self.is_local_run = bool(self.additional_attributes.get("run_locally"))
        self.genome_name = self.additional_attributes.get("genome_name", "nt_k16")
        self.index = self.additional_files.get("index")
        if self.is_local_run:
            assert self.index, "local runs require an index to be passed in"
        else:
            assert not self.index, "passing in an index is not supported for remote runs"
            command.make_dirs(self.chunks_result_dir_local)
            self.batch_job_desc_bucket = get_batch_job_desc_bucket()

    def count_reads(self):
        pass

    def run(self):
        ''' Run alignmment remotely '''

        alignment_algorithm_inputs = PipelineStepRunAlignment._alignment_algorithm_inputs(self.input_files_local[0])
        duplicate_cluster_sizes_path, = self.input_files_local[1]
        output_m8, deduped_output_m8, output_hitsummary, output_counts_with_dcr_json = self.output_files_local()
        assert output_counts_with_dcr_json.endswith("_with_dcr.json"), self.output_files_local()

        if self.is_local_run:
            self.run_locally(alignment_algorithm_inputs[self.alignment_algorithm], output_m8)
        else:
            self.run_remotely(alignment_algorithm_inputs[self.alignment_algorithm], output_m8)

        # get database
        lineage_db = fetch_reference(self.additional_files["lineage_db"], self.ref_dir_local)
        accession2taxid_db = fetch_reference(self.additional_files["accession2taxid_db"], self.ref_dir_local, allow_s3mi=True)

        deuterostome_db = None
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = fetch_reference(self.additional_files["deuterostome_db"],
                                              self.ref_dir_local, allow_s3mi=True)

        blacklist_s3_file = self.additional_files.get('taxon_blacklist', DEFAULT_BLACKLIST_S3)
        taxon_blacklist = fetch_reference(blacklist_s3_file, self.ref_dir_local)

        taxon_whitelist = None
        if self.additional_attributes.get("use_taxon_whitelist"):
            taxon_whitelist = fetch_reference(self.additional_files.get("taxon_whitelist", DEFAULT_WHITELIST_S3),
                                              self.ref_dir_local)

        min_alignment_length = NT_MIN_ALIGNMENT_LEN if self.alignment_algorithm == 'gsnap' else 0
        m8.call_hits_m8(output_m8, lineage_db, accession2taxid_db,
                        deduped_output_m8, output_hitsummary, min_alignment_length,
                        deuterostome_db, taxon_whitelist, taxon_blacklist)

        db_type = 'NT' if self.alignment_algorithm == 'gsnap' else 'NR'

        m8.generate_taxon_count_json_from_m8(
            deduped_output_m8, output_hitsummary, db_type,
            lineage_db, deuterostome_db, taxon_whitelist, taxon_blacklist, duplicate_cluster_sizes_path,
            output_counts_with_dcr_json)

    def run_locally(self, input_fas, output_m8):
        with tempfile.TemporaryDirectory(prefix=self.alignment_algorithm) as tmpdir:
            index_path = os.path.join(tmpdir, "reference")
            os.mkdir(index_path)
            run(["tar", "-xzvf", self.index, "-C", index_path], check=True)

            if self.alignment_algorithm == "gsnap":
                # Hack to determine gsnap vs gsnapl
                error_message = run(
                    ['gsnapl.avx2-2018-10-26', '-D', index_path, '-d', self.genome_name],
                    input='>'.encode('utf-8'),
                    stderr=PIPE,
                    stdout=PIPE
                ).stderr
                # note, alignment uses a pinned version of gmap/gsnap
                gsnap_command = "gsnap.avx2-2018-10-26" if 'please run gsnap instead' in error_message.decode('utf-8') else "gsnapl.avx2-2018-10-26"
            else:
                gsnap_command = None

            # rapsearch2 expects the filename of the primary file, but still depends on the .info file being a sibbling
            if self.alignment_algorithm == "rapsearch2":
                for filename in os.listdir(index_path):
                    if not filename.endswith(".info"):
                        index_path = os.path.join(index_path, filename)

            cmd = self._get_command(
                index_path,
                input_fas,
                output_m8,
                threads=multiprocessing.cpu_count(),
                gsnap_command=gsnap_command
            )
            log.write(f"running command {cmd}")
            run(cmd, check=True)

    def run_remotely(self, input_fas, output_m8):
        # Split files into chunks for performance
        chunk_size = GSNAP_CHUNK_SIZE if self.alignment_algorithm == "gsnap" else RAPSEARCH_CHUNK_SIZE
        part_suffix, input_chunks = self.chunk_input(input_fas, chunk_size)
        self.chunk_count = len(input_chunks)

        # Process chunks
        chunk_output_files = [None] * self.chunk_count
        chunk_threads = []
        mutex = TraceLock("run_remotely", threading.RLock())
        # Randomize execution order for performance
        randomized = list(enumerate(input_chunks))
        random.shuffle(randomized)

        try:
            for n, chunk_input_files in randomized:
                self.chunks_in_flight_semaphore.acquire()
                self.check_for_errors(mutex, chunk_output_files, input_chunks, self.alignment_algorithm)
                t = threading.Thread(
                    target=PipelineStepRunAlignment.run_chunk_wrapper,
                    kwargs={
                        'chunks_in_flight_semaphore': self.chunks_in_flight_semaphore,
                        'chunk_output_files': chunk_output_files,
                        'n': n,
                        'mutex': mutex,
                        'target': self.run_chunk,
                        'kwargs': {
                            'part_suffix': part_suffix,
                            'input_files': chunk_input_files,
                            'lazy_run': False,
                        },
                    })
                t.start()
                chunk_threads.append(t)

        finally:
            # Check chunk completion
            for ct in chunk_threads:
                ct.join()

        self.check_for_errors(mutex, chunk_output_files, input_chunks, self.alignment_algorithm)

        assert None not in chunk_output_files
        # Concatenate the pieces and upload results
        self.concatenate_files(chunk_output_files, output_m8)

    def chunk_input(self, input_files, chunksize):
        """Chunk input files into pieces for performance and parallelism."""
        part_lists = []  # Lists of partial files
        known_nlines = None
        part_suffix = ""
        chunk_nlines = chunksize * 2

        for input_file in input_files:
            # Count number of lines in the file
            cmd_output = command.execute_with_output(
                command_patterns.SingleCommand(
                    cmd="wc",
                    args=[
                        "-l",
                        input_file
                    ]
                )
            )
            nlines = int(cmd_output.strip().split()[0])
            # Number of lines should be the same in paired files
            if known_nlines is not None:
                msg = "Mismatched line counts in supposedly paired files: {}".format(
                    input_files)
                assert nlines == known_nlines, msg
            known_nlines = nlines

            # Set number of pieces and names
            numparts = (nlines + chunk_nlines - 1) // chunk_nlines
            ndigits = len(str(numparts - 1))
            part_suffix = "-chunksize-%d-numparts-%d-part-" % (chunksize, numparts)
            out_prefix_base = os.path.basename(input_file) + part_suffix
            out_prefix = os.path.join(self.chunks_result_dir_local, out_prefix_base)

            # Split large file into smaller named pieces
            command.execute(
                command_patterns.SingleCommand(
                    cmd="split",
                    args=[
                        "-a",
                        ndigits,
                        "--numeric-suffixes",
                        "-l",
                        chunk_nlines,
                        input_file,
                        out_prefix
                    ]
                )
            )
            command.execute_with_retries(
                command_patterns.SingleCommand(
                    cmd="aws",
                    args=[
                        "s3",
                        "sync",
                        "--only-show-errors",
                        os.path.join(self.chunks_result_dir_local, ""),
                        os.path.join(self.chunks_result_dir_s3, ""),
                        "--exclude",
                        "*",
                        "--include",
                        out_prefix_base + "*"
                    ]
                )
            )

            # Get the partial file names
            partial_files = []
            paths = command.glob(
                glob_pattern=out_prefix + "*",
                strip_folder_names=True
            )
            partial_files.extend(paths)

            # Check that the partial files match our expected chunking pattern
            pattern = "{:0%dd}" % ndigits
            expected_partial_files = [(out_prefix_base + pattern.format(i))
                                      for i in range(numparts)]
            msg = "something went wrong with chunking: {} != {}".format(
                partial_files, expected_partial_files)
            assert expected_partial_files == partial_files, msg
            part_lists.append(partial_files)

        # Ex: [["input_R1.fasta-part-1", "input_R2.fasta-part-1"],
        # ["input_R1.fasta-part-2", "input_R2.fasta-part-2"],
        # ["input_R1.fasta-part-3", "input_R2.fasta-part-3"],...]
        input_chunks = [list(part) for part in zip(*part_lists)]
        return part_suffix, input_chunks

    @staticmethod
    def check_for_errors(mutex, chunk_output_files, input_chunks, alignment_algorithm):
        with mutex:
            if "error" in chunk_output_files:
                # We already have per-chunk retries to address transient (non-deterministic) system issues.
                # If a chunk fails on all retries, it must be due to a deterministic & reproducible problem
                # with the chunk data or command options, so we should not even attempt the other chunks.
                err_i = chunk_output_files.index("error")
                raise RuntimeError("All retries failed for {} chunk {}.".format(
                    alignment_algorithm, input_chunks[err_i]))

    @staticmethod
    def concatenate_files(chunk_output_files, output_m8):
        with log.log_context("run_alignment_remote.concatenate_files", {"chunk_output_files": chunk_output_files}):
            with open(output_m8, 'wb') as outf:
                for f in chunk_output_files:
                    with log.log_context("run_alignment_remote.concatenate_files#chunk", {"f": f}):
                        with open(f, 'rb') as fd:
                            shutil.copyfileobj(fd, outf)

    @staticmethod
    def run_chunk_wrapper(chunks_in_flight_semaphore, chunk_output_files, n, mutex, target, kwargs):
        result = "error"
        try:
            result = target(**kwargs)
        except:
            with mutex:
                log.write(traceback.format_exc())
        finally:
            with mutex:
                chunk_output_files[n] = result
            chunks_in_flight_semaphore.release()

    def _get_job_status(self, session, job_id):
        batch_job_desc_bucket = session.resource("s3").Bucket(self.batch_job_desc_bucket)
        key = f"job_descriptions/{job_id}"
        try:
            job_desc_object = batch_job_desc_bucket.Object(key)
            return json.loads(job_desc_object.get()["Body"].read())["status"]
        except ClientError as e:
            if e.response['Error']['Code'] == 'NoSuchKey':
                # Warn that the object is missing so any issue with the s3 mechanism can be identified
                log.log_event("missing_job_description_ojbect", values={key: key}, debug=True)
                # Return submitted because a missing job status probably means it hasn't been added yet
                return "SUBMITTED"
            else:
                raise e

    def _log_alignment_batch_job_status(self, job_id, job_queue, job_definition, chunk_id, status):
        log.log_event('alignment_batch_job_status', values={
            'job_id': job_id,
            'chunk_id': chunk_id,
            'job_queue': job_queue,
            'job_definition': job_definition,
            'status': status,
            'alignment_algorithm': self.alignment_algorithm,
        })

    def _run_batch_job(self, session, job_name, job_queue, job_definition, command, environment, chunk_id, retries):
        client = session.client("batch")
        response = client.submit_job(
            jobName=job_name,
            jobQueue=job_queue,
            jobDefinition=job_definition,
            containerOverrides={
                "command": command,
                "environment": environment,
            },
            retryStrategy={"attempts": retries}
        )
        job_id = response["jobId"]
        self._log_alignment_batch_job_status(job_id, job_queue, job_definition, chunk_id, 'SUBMITTED')

        chunks_in_flight = min(self.chunk_count, MAX_CHUNKS_IN_FLIGHT)  # use min in case we have fewer chunks than the maximum
        mean_delay = chunks_in_flight  # ~1 chunk per second to avoid throttling,
        delay = mean_delay + random.randint(-mean_delay // 2, mean_delay // 2)  # Add some noise to de-synchronize chunks
        status = "SUBMITTED"
        # the job this is monitoring has an timeout and the job this runs in has a timeout
        while True:
            try:
                status = self._get_job_status(session, job_id)
            except ClientError as e:
                # If we get throttled, randomly wait to de-synchronize the requests
                if e.response['Error']['Code'] == "TooManyRequestsException":
                    log.log_event("describe_jobs_rate_limit_error", values={"job_id": job_id}, warning=True)
                    # Possibly implement a backoff here if throttling becomes an issue
                else:
                    log.log_event("unexpected_client_error_while_polling_job_status", values={"job_id": job_id})
                    raise e

            if status == "SUCCEEDED":
                self._log_alignment_batch_job_status(job_id, job_queue, job_definition, chunk_id, status)
                return job_id
            if status == "FAILED":
                log.log_event("alignment_batch_job_failed", values={'job_id': job_id, 'chunk_id': chunk_id, 'alignment_algorithm': self.alignment_algorithm})
                self._log_alignment_batch_job_status(job_id, job_queue, job_definition, chunk_id, status)
                raise BatchJobFailed("chunk alignment failed")
            time.sleep(delay)

    def _validate_chunk_output(self, chunk_output_filename):
        cmd = "awk '{print NF}' | sort -nu | head -n 1"
        if self.alignment_algorithm == "rapsearch2":
            cmd = "grep -v '^#' | " + cmd
        with open(chunk_output_filename, "rb") as f:
            line_count_str = run(cmd, stdin=f, stdout=PIPE, check=True, shell=True).stdout
        if line_count_str:
            line_count = int(line_count_str)
            log.write(f"smallest number of columns observed in any line in output file {chunk_output_filename} was {line_count}")
            assert line_count == 12, "output file {chunk_output_filename} was not a valid m8 file"
        else:
            log.write(f"no hits in output file {chunk_output_filename}")

    def _get_command(self, index_path, input_paths, output_path, threads=None, gsnap_command="gsnapl"):
        if not threads:
            threads = 48 if self.alignment_algorithm == "gsnap" else 24
        if self.alignment_algorithm == "gsnap":
            return [gsnap_command,
                    "-A", "m8",
                    "--batch=0",
                    "--use-shared-memory=0",
                    "--gmap-mode=none",
                    "--npaths=100",
                    "--ordered",
                    "-t", str(threads),
                    "--max-mismatches=40",
                    "-D", index_path,
                    "-d", self.genome_name,
                    "-o", output_path,
                    ] + input_paths
        else:
            output_path = re.sub(r'\.m8$', '', output_path)
            return ["rapsearch",
                    "-d", index_path,
                    "-e", "-6",
                    "-l", "10",
                    "-a", "T",
                    "-b", "0",
                    "-v", "50",
                    "-z", str(threads),
                    "-o", output_path,
                    "-q"] + input_paths

    def run_chunk(self, part_suffix, input_files, lazy_run):
        """
        Dispatch a chunk to worker machines for distributed GSNAP or RAPSearch
        group machines and handle their execution.
        """

        chunk_id = int(input_files[0].split(part_suffix)[-1])
        multihit_basename = f"multihit-{self.alignment_algorithm}-out{part_suffix}{chunk_id}.m8"
        multihit_local_outfile = os.path.join(self.chunks_result_dir_local, multihit_basename)
        multihit_s3_outfile = os.path.join(self.chunks_result_dir_s3, multihit_basename)

        session = boto3.session.Session()
        if lazy_run and download_from_s3(session, multihit_s3_outfile, multihit_local_outfile):
            log.write(f"finished alignment for chunk {chunk_id} with {self.alignment_algorithm} by lazily fetching last result")
            return multihit_local_outfile

        deployment_environment = os.environ["DEPLOYMENT_ENVIRONMENT"]
        priority_name = os.environ.get("PRIORITY_NAME", "normal")

        index_dir_suffix = self.additional_attributes["index_dir_suffix"]

        pattern = r's3://.+/samples/([0-9]+)/([0-9]+)/'
        m = re.match(pattern, self.chunks_result_dir_s3)
        if m:
            project_id, sample_id = m.group(1), m.group(2)
        else:
            project_id, sample_id = '0', '0'

        provisioning_model = 'SPOT'

        # TODO: parameterize these https://jira.czi.team/browse/IDSEQ-2673
        job_name = f"idseq-{deployment_environment}-{self.alignment_algorithm}-project-{project_id}-sample-{sample_id}-part-{chunk_id}"
        job_queue = f"idseq-{deployment_environment}-{self.alignment_algorithm}-{provisioning_model}-{index_dir_suffix}-{priority_name}"
        job_definition = f"idseq-{deployment_environment}-{self.alignment_algorithm}"

        environment = [{
            'name': f"INPUT_PATH_{i}",
            'value': os.path.join(self.chunks_result_dir_s3, input_file),
        } for i, input_file in enumerate(input_files)] + [{
            'name': "OUTPUT_PATH",
            'value': multihit_s3_outfile,
        }]

        input_paths = [f"local_input_{i}" for i, _ in enumerate(input_files)]
        command = self._get_command("/references/reference", input_paths, "local_output")

        try:
            job_id = self._run_batch_job(
                session=session,
                job_name=job_name,
                job_queue=job_queue,
                job_definition=job_definition,
                command=command,
                environment=environment,
                chunk_id=chunk_id,
                retries=2
            )
        except BatchJobFailed:
            provisioning_model = 'EC2'
            # TODO: parameterize this https://jira.czi.team/browse/IDSEQ-2673
            job_queue = f"idseq-{deployment_environment}-{self.alignment_algorithm}-{provisioning_model}-{index_dir_suffix}-{priority_name}"
            job_id = self._run_batch_job(
                session=session,
                job_name=job_name,
                job_queue=job_queue,
                job_definition=job_definition,
                command=command,
                environment=environment,
                chunk_id=chunk_id,
                retries=1
            )

        for _ in range(12):
            if download_from_s3(session, multihit_s3_outfile, multihit_local_outfile):
                break
            time.sleep(10)
        else:
            log.log_event("chunk_result_missing_in_s3", values={
                'job_id': job_id,
                'chunk_id': chunk_id,
                'job_queue': job_queue,
                'job_definition': job_definition,
                'alignment_algorithm': self.alignment_algorithm,
            })
            raise Exception("Chunk result is missing from s3")

        log.log_event("alignment_batch_chunk_result_downloaded", values={
            'job_id': job_id,
            'chunk_id': chunk_id,
            'job_queue': job_queue,
            'job_definition': job_definition,
            'alignment_algorithm': self.alignment_algorithm,
        })

        self._validate_chunk_output(multihit_local_outfile)
        return multihit_local_outfile

    def step_description(self, require_docstrings=False):
        if (self.alignment_algorithm == "gsnap"):
            return """
                Runs gsnap remotely.

                ```
                gsnapl
                -A m8
                --batch=0
                --use-shared-memory=0
                --gmap-mode=none
                --npaths=100
                --ordered
                -t 48
                --max-mismatches=40
                -D {remote_index_dir}
                -d nt_k16
                {remote_input_files} > {multihit_remote_outfile}
                ```

                GSNAP documentation is available [here](http://research-pub.gene.com/gmap/).
            """
        elif self.alignment_algorithm == "rapsearch2":
            return """
                Runs rapsearch2 remotely.

                ```
                rapsearch
                -d {remote_index_dir}/nr_rapsearch
                -e -6
                -l 10
                -a T
                -b 0
                -s f
                -v 50
                -z 24
                -q {remote_input_files}
                -o {multihit_remote_outfile}
                ```

                Rapsearch2 documentation is available [here](http://omics.informatics.indiana.edu/mg/RAPSearch2/).
            """
        # If neither, then return default step_description method.
        return super(PipelineStepRunAlignment, self).step_description()
