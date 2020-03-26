import multiprocessing
import os
import random
import shlex
import shutil
import threading
import time
import traceback

from idseq_dag.engine.pipeline_step import InputFileErrors, PipelineStep

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.log as log
import idseq_dag.util.m8 as m8

from idseq_dag.util.s3 import fetch_from_s3, fetch_reference
from idseq_dag.util.server import ASGInstance, ChunkStatus, chunk_status_tracker
from idseq_dag.util.trace_lock import TraceLock

from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3, DEFAULT_WHITELIST_S3
from idseq_dag.util.m8 import NT_MIN_ALIGNMENT_LEN

MAX_CONCURRENT_CHUNK_UPLOADS = 4
CORRECT_NUMBER_OF_OUTPUT_COLUMNS = 12
CHUNK_MAX_TRIES = 3

# Please override this with gsnap_chunk_timeout or rapsearch_chunk_timeout in DAG json.
# Default is several sigmas beyond the pale and indicates the data has to be QC-ed better.
# Note(2020-01-10): Raised to 3 hrs to mitigate Rapsearch chunk timeouts after recent index update.
DEFAULT_CHUNK_TIMEOUT = 60 * 60 * 3

class PipelineStepRunAlignmentRemotely(PipelineStep):
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
    -t (threads): r5d.metal machines have 96 vCPUs. Use 48 threads and each process will be able to
    concurrently process 2 chunks (see attribute 'max_concurrent').

    For Rapsearch:
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
    def _service_inputs(host_filter_outputs):
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
            self.input_file_error = InputFileErrors.INSUFFICIENT_READS

    def __init__(self, *args, **kwrds):
        PipelineStep.__init__(self, *args, **kwrds)
        self.chunks_in_flight = threading.Semaphore(self.additional_attributes['chunks_in_flight'])
        self.chunks_result_dir_local = os.path.join(self.output_dir_local, "chunks")
        self.chunks_result_dir_s3 = os.path.join(self.output_dir_s3, "chunks")
        self.iostream_upload = multiprocessing.Semaphore(MAX_CONCURRENT_CHUNK_UPLOADS)
        command.make_dirs(self.chunks_result_dir_local)

    def count_reads(self):
        pass

    def run(self):
        ''' Run alignmment remotely '''

        service_inputs = PipelineStepRunAlignmentRemotely._service_inputs(self.input_files_local[0])
        cdhit_cluster_sizes_path, = self.input_files_local[1]
        output_m8, deduped_output_m8, output_hitsummary, output_counts_with_dcr_json = self.output_files_local()
        assert output_counts_with_dcr_json.endswith("_with_dcr.json"), self.output_files_local()

        service = self.additional_attributes["service"]
        self.run_remotely(service_inputs[service], output_m8, service)

        # get database
        lineage_db = fetch_reference(self.additional_files["lineage_db"], self.ref_dir_local)
        accession2taxid_db = fetch_reference(self.additional_files["accession2taxid_db"], self.ref_dir_local, allow_s3mi=True)

        min_alignment_length = NT_MIN_ALIGNMENT_LEN if service == 'gsnap' else 0
        m8.call_hits_m8(output_m8, lineage_db, accession2taxid_db,
                        deduped_output_m8, output_hitsummary, min_alignment_length)

        db_type = 'NT' if service == 'gsnap' else 'NR'
        evalue_type = 'log10' if service == 'rapsearch2' else 'raw'

        deuterostome_db = None
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = fetch_reference(self.additional_files["deuterostome_db"],
                                              self.ref_dir_local, allow_s3mi=True)

        blacklist_s3_file = self.additional_attributes.get('taxon_blacklist', DEFAULT_BLACKLIST_S3)
        taxon_blacklist = fetch_reference(blacklist_s3_file, self.ref_dir_local)

        taxon_whitelist = None
        if self.additional_attributes.get("use_taxon_whitelist"):
            taxon_whitelist = fetch_reference(self.additional_files.get("taxon_whitelist", DEFAULT_WHITELIST_S3),
                                              self.ref_dir_local)

        m8.generate_taxon_count_json_from_m8(
            deduped_output_m8, output_hitsummary, evalue_type, db_type,
            lineage_db, deuterostome_db, taxon_whitelist, taxon_blacklist, cdhit_cluster_sizes_path,
            output_counts_with_dcr_json)

    def run_remotely(self, input_fas, output_m8, service):
        key_path = self.fetch_key(os.environ['KEY_PATH_S3'])
        sample_name = self.output_dir_s3.rstrip('/').replace('s3://', '').replace('/', '-')
        chunk_size = int(self.additional_attributes["chunk_size"])
        index_dir_suffix = self.additional_attributes.get("index_dir_suffix")
        remote_username = "ec2-user"
        remote_home_dir = os.path.join("/home", remote_username)
        remote_index_dir = os.path.join(remote_home_dir, "references")

        if index_dir_suffix:
            remote_index_dir = os.path.join(remote_index_dir, index_dir_suffix)

        sample_remote_work_dir = os.path.join(remote_home_dir, "batch-pipeline-workdir", sample_name)

        # Split files into chunks for performance
        part_suffix, input_chunks = self.chunk_input(input_fas, chunk_size)

        # Process chunks
        chunk_output_files = [None] * len(input_chunks)
        chunk_threads = []
        mutex = TraceLock("run_remotely", threading.RLock())
        # Randomize execution order for performance
        randomized = list(enumerate(input_chunks))
        random.shuffle(randomized)

        try:
            for n, chunk_input_files in randomized:
                self.chunks_in_flight.acquire()
                self.check_for_errors(mutex, chunk_output_files, input_chunks, service)
                chunk_remote_work_dir = f"{sample_remote_work_dir}-chunk-{n}"
                t = threading.Thread(
                    target=PipelineStepRunAlignmentRemotely.run_chunk_wrapper,
                    args=[
                        self.chunks_in_flight, chunk_output_files, n, mutex, self.run_chunk,
                        [
                            part_suffix, remote_home_dir, remote_index_dir,
                            chunk_remote_work_dir, remote_username, chunk_input_files,
                            key_path, service, True
                        ]
                    ])
                t.start()
                chunk_threads.append(t)

        finally:
            # Check chunk completion
            for ct in chunk_threads:
                ct.join()
            try:
                chunk_status_tracker(service).log_stats(len(input_chunks))
            except:
                log.write(f"Problem dumping status report for {service}")
                log.write(traceback.format_exc())

        self.check_for_errors(mutex, chunk_output_files, input_chunks, service)

        assert None not in chunk_output_files
        # Concatenate the pieces and upload results
        self.concatenate_files(chunk_output_files, output_m8)

    def fetch_key(self, key_path_s3):
        key_path = fetch_from_s3(key_path_s3, self.output_dir_local)
        command.chmod(key_path, 0o400)
        return key_path

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
    def check_for_errors(mutex, chunk_output_files, input_chunks, service):
        with mutex:
            if "error" in chunk_output_files:
                # We already have per-chunk retries to address transient (non-deterministic) system issues.
                # If a chunk fails on all retries, it must be due to a deterministic & reproducible problem
                # with the chunk data or command options, so we should not even attempt the other chunks.
                err_i = chunk_output_files.index("error")
                raise RuntimeError("All retries failed for {} chunk {}.".format(
                    service, input_chunks[err_i]))

    @staticmethod
    def concatenate_files(chunk_output_files, output_m8):
        with log.log_context("run_alignment_remote.concatenate_files", {"chunk_output_files": chunk_output_files}):
            with open(output_m8, 'wb') as outf:
                for f in chunk_output_files:
                    with log.log_context("run_alignment_remote.concatenate_files#chunk", {"f": f}):
                        with open(f, 'rb') as fd:
                            shutil.copyfileobj(fd, outf)

    @staticmethod
    def run_chunk_wrapper(chunks_in_flight, chunk_output_files, n, mutex, target, args):
        result = "error"
        try:
            result = target(*args)
        except:
            with mutex:
                log.write(traceback.format_exc())
        finally:
            with mutex:
                chunk_output_files[n] = result
            chunks_in_flight.release()

    @staticmethod
    def __interpret_min_column_number_string(min_column_number_string,
                                             correct_number_of_output_columns,
                                             try_number):
        if min_column_number_string:
            min_column_number = float(min_column_number_string)
            log.write(
                "Try no. %d: Smallest number of columns observed in any line was %d"
                % (try_number, min_column_number))
        else:
            log.write("Try no. %d: No hits" % try_number)
            min_column_number = correct_number_of_output_columns
        return min_column_number

    @command.retry
    def __check_if_output_is_corrupt(self, service, key_path, remote_username, instance_ip,  # self unused
                                     multihit_remote_outfile, chunk_id, try_number):
        # Check if every row has correct number of columns (12) in the output
        # file on the remote machine
        if service == "gsnap":
            verification_command = "cat %s" % shlex.quote(multihit_remote_outfile)
        else:
            # For rapsearch, first remove header lines starting with '#'
            verification_command = "grep -v '^#' %s" % shlex.quote(multihit_remote_outfile)
        verification_command += " | awk '{print NF}' | sort -nu | head -n 1"
        min_column_number_string = command.execute_with_output(
            command.remote(verification_command, key_path, remote_username, instance_ip))
        min_column_number = PipelineStepRunAlignmentRemotely.__interpret_min_column_number_string(
            min_column_number_string, CORRECT_NUMBER_OF_OUTPUT_COLUMNS,
            try_number)
        error = None
        if min_column_number != CORRECT_NUMBER_OF_OUTPUT_COLUMNS:
            msg = "Chunk %s output corrupt; not copying to S3. min_column_number = %d -> expected = %d."
            msg += " Re-start pipeline to try again."
            error = msg % (chunk_id, min_column_number, CORRECT_NUMBER_OF_OUTPUT_COLUMNS)
        return error

    @command.retry
    def __copy_multihit_remote_outfile(self, key_path, remote_username, instance_ip,  # self unused
                                       multihit_remote_outfile, multihit_local_outfile):
        # Copy output from remote machine to local machine.  This needs to happen while we are
        # holding the machine reservation, i.e., inside the "with ASGInstnace" context.
        command.execute(
            command.scp(
                key_path, remote_username, instance_ip,
                multihit_remote_outfile,
                multihit_local_outfile
            )
        )

    @command.retry
    def __delete_remote_dir(self, remote_dir, key_path, remote_username, instance_ip):
        """
        Delete a directory on a remote machine
        This needs to happen while we are holding the machine reservation,
        i.e., inside the "with ASGInstnace" context.
        """
        rm_command = f"rm -rf {remote_dir}"
        command.execute(command.remote(rm_command, key_path, remote_username, instance_ip))

    def run_chunk(self, part_suffix, remote_home_dir, remote_index_dir, remote_work_dir,
                  remote_username, input_files, key_path, service, lazy_run):
        """Dispatch a chunk to worker machines for distributed GSNAP or RAPSearch
        group machines and handle their execution.
        """
        assert service in ("gsnap", "rapsearch2")

        chunk_id = int(input_files[0].split(part_suffix)[-1])
        multihit_basename = f"multihit-{service}-out{part_suffix}{chunk_id}.m8"
        multihit_local_outfile = os.path.join(self.chunks_result_dir_local, multihit_basename)
        multihit_remote_outfile = os.path.join(remote_work_dir, multihit_basename)
        multihit_s3_outfile = os.path.join(self.chunks_result_dir_s3, multihit_basename)

        def aws_cp_operation(input_fa):
            return "aws s3 cp --only-show-errors {src} {dest}".format(
                src=shlex.quote(os.path.join(self.chunks_result_dir_s3, input_fa)),
                dest=shlex.quote(os.path.join(remote_work_dir, input_fa))
            )

        download_input_from_s3 = " ; ".join(map(aws_cp_operation, input_files))

        # Clean up remote work directory before running
        #   This ensures that files from a failed previous run that may still be on the instance
        #   are removed so they don't corrupt the current run
        base_str = "rm -rf {remote_work_dir} ; mkdir -p {remote_work_dir} ; {download_input_from_s3} ; "
        environment = self.additional_attributes["environment"]

        # See step class docstrings for more parameter details.
        if service == "gsnap":
            commands = base_str + "{remote_home_dir}/bin/gsnapl -A m8 --batch=0 --use-shared-memory=0 --gmap-mode=none --npaths=100 --ordered -t 48 --max-mismatches=40 -D {remote_index_dir} -d nt_k16 {remote_input_files} > {multihit_remote_outfile}"
        else:
            commands = base_str + "/usr/local/bin/rapsearch -d {remote_index_dir}/nr_rapsearch -e -6 -l 10 -a T -b 0 -v 50 -z 24 -q {remote_input_files} -o {multihit_remote_outfile}"

        commands = commands.format(
            remote_work_dir=shlex.quote(remote_work_dir),
            download_input_from_s3=download_input_from_s3,
            remote_home_dir=shlex.quote(remote_home_dir),
            remote_index_dir=shlex.quote(remote_index_dir),
            remote_input_files=" ".join(shlex.quote(remote_work_dir + "/" + input_fa) for input_fa in input_files),
            multihit_remote_outfile=shlex.quote(multihit_remote_outfile) if service == "gsnap" else shlex.quote(multihit_remote_outfile[:-3])
            # Strip the .m8 for RAPSearch as it adds that
        )

        if lazy_run and fetch_from_s3(multihit_s3_outfile, multihit_local_outfile, okay_if_missing=True, allow_s3mi=False):
            log.write(f"finished alignment for chunk {chunk_id} with {service} by lazily fetching last result")
        else:
            chunk_timeout = int(self.additional_attributes.get(f"{service.lower()}_chunk_timeout",
                                                               DEFAULT_CHUNK_TIMEOUT))
            for try_number in range(1, CHUNK_MAX_TRIES + 1):
                log.write(f"waiting for {service} server for chunk {chunk_id}. Try #{try_number}")
                with ASGInstance(service, key_path, remote_username, environment, chunk_id, try_number,
                                 self.additional_attributes) as instance_ip:
                    # Try/Except block needs to be inside the ASGInstance context.
                    # A failure to acquire an ASGInstnace is and should be unrecoverable.
                    chunk_status = None
                    elapsed = 0.0
                    try:
                        t_start = time.time()
                        try:
                            command.execute(command.remote(commands, key_path, remote_username, instance_ip),
                                            timeout=chunk_timeout)
                        except:
                            chunk_status = ChunkStatus.CRASH
                            raise
                        finally:
                            elapsed = time.time() - t_start
                            if chunk_status == ChunkStatus.CRASH and elapsed >= chunk_timeout:
                                chunk_status = ChunkStatus.TIMEOUT

                        output_corrupt = self.__check_if_output_is_corrupt(
                            service, key_path, remote_username, instance_ip, multihit_remote_outfile,
                            chunk_id, try_number)

                        if output_corrupt:
                            chunk_status = ChunkStatus.CORRUPT_OUTPUT
                            assert not output_corrupt, output_corrupt

                        # Yay, chunk succeeded.  Copy from server and break out of retry loop.
                        try:
                            self.__copy_multihit_remote_outfile(
                                key_path, remote_username, instance_ip,
                                multihit_remote_outfile, multihit_local_outfile)
                            chunk_status = ChunkStatus.SUCCESS
                            break
                        except:
                            # If we failed to copy from the server, it's as bad as a crash in alignment.
                            chunk_status = ChunkStatus.CRASH
                            raise

                    except Exception as e:

                        # 1. No backoff needed here before retrying.  We rate limit chunk dispatch (the ASGInstance
                        # acquisition above is blocking).  ASGInstance acquisition also tries to ensure that every
                        # chunk flight gets its first try before any retry is dispatched.

                        # 2. If the reason we failed is timeout on the server, we don't retry.  The operator must decide
                        # whether to QC the data more, or use smaller chunk size.  In fact, we only retry for CRASH and
                        # CORRUPT_OUTPUT.

                        # 3. If this is the last attempt, we gotta re-raise the exception.

                        # 4. Elapsed time is only the time spent in alignment.  It excludes the time spent waiting to
                        # acquire ASGinstance.

                        log.log_event('alignment_remote_error', values={"chunk": chunk_id,
                                                                        "try_number": try_number,
                                                                        "CHUNK_MAX_TRIES": CHUNK_MAX_TRIES,
                                                                        "chunk_status": chunk_status,
                                                                        "elapsed": elapsed,
                                                                        "chunk_timeout": chunk_timeout,
                                                                        "exception": log.parse_exception(e)})
                        retrying_might_help = chunk_status in (ChunkStatus.CORRUPT_OUTPUT, ChunkStatus.CRASH)
                        if try_number < CHUNK_MAX_TRIES and retrying_might_help:
                            # Retry!
                            continue
                        else:
                            # End of the road.
                            raise
                    finally:
                        # None chunk_status indicates code bug above.  An exception has been raised already
                        # for it, and it says nothing about whether the alignment succeeded or not.
                        if chunk_status != None:
                            chunk_status_tracker(service).note_outcome(instance_ip, chunk_id, elapsed, chunk_status, try_number)
                        self.__delete_remote_dir(remote_work_dir, key_path, remote_username, instance_ip)

            # Upload to s3
            with self.iostream_upload:  # Limit concurrent uploads so as not to stall the pipeline.
                command.execute(
                    command_patterns.SingleCommand(
                        cmd="aws",
                        args=[
                            "s3",
                            "cp",
                            "--only-show-errors",
                            multihit_local_outfile,
                            os.path.join(self.chunks_result_dir_s3, "")
                        ]
                    )
                )
            log.write(f"finished alignment for chunk {chunk_id} on {service} server {instance_ip}")

        # Whether lazy or not lazy, we've now got the chunk result locally here.
        return multihit_local_outfile

    def step_description(self, require_docstrings=False):
        if (self.name == "gsnap_out"):
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
        elif (self.name == "rapsearch2_out"):
            return """
                Runs rapsearch remotely.

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
        # If neither, then return default step_description method.
        return super(PipelineStepRunAlignmentRemotely, self).step_description()
