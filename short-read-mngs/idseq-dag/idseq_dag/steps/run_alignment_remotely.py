import os
import threading
import shutil
import random
import traceback
import multiprocessing

from idseq_dag.engine.pipeline_step import PipelineStep

import idseq_dag.util.command as command
import idseq_dag.util.server as server
import idseq_dag.util.log as log
import idseq_dag.util.m8 as m8
from idseq_dag.util.s3 import fetch_from_s3

MAX_CONCURRENT_CHUNK_UPLOADS = 4
DEFAULT_BLACKLIST_S3 = 's3://idseq-database/taxonomy/2018-04-01-utc-1522569777-unixtime__2018-04-04-utc-1522862260-unixtime/taxon_blacklist.txt'

class PipelineStepRunAlignmentRemotely(PipelineStep):
    '''
    Run gsnap/rapsearch2 remotely
    '''

    def __init__(self, *args, **kwrds):
        PipelineStep.__init__(self, *args, **kwrds)
        self.chunks_in_flight = threading.Semaphore(self.additional_attributes['chunks_in_flight'])
        self.chunks_result_dir_local = os.path.join(self.output_dir_local, "chunks")
        self.chunks_result_dir_s3 = os.path.join(self.output_dir_s3, "chunks")
        self.iostream_upload = multiprocessing.Semaphore(MAX_CONCURRENT_CHUNK_UPLOADS)

        command.execute("mkdir -p %s" % self.chunks_result_dir_local)

    def count_reads(self):
        pass

    def run(self):
        ''' Run alignmment remotely '''
        input_fas = self.get_input_fas()
        [output_m8, deduped_output_m8, output_hitsummary, output_counts_json] = self.output_files_local()
        service = self.additional_attributes["service"]
        assert service in ("gsnap", "rapsearch2")

        # TODO: run the alignment remotely and make lazy_chunk=True, revisit this later
        self.run_remotely(input_fas, output_m8, service)

        # get database
        lineage_db = fetch_from_s3(self.additional_files["lineage_db"], self.ref_dir_local, allow_s3mi=True)
        accession2taxid_db = fetch_from_s3(self.additional_files["accession2taxid_db"], self.ref_dir_local, allow_s3mi=True)
        blacklist_s3_file = self.additional_attributes.get('taxon_blacklist', DEFAULT_BLACKLIST_S3)
        taxon_blacklist = fetch_from_s3(blacklist_s3_file, self.ref_dir_local)
        m8.call_hits_m8(output_m8, lineage_db, accession2taxid_db,
                        deduped_output_m8, output_hitsummary, taxon_blacklist)

        # check deuterostome
        deuterostome_db = None
        db_type = 'NT' if service == 'gsnap' else 'NR'
        evalue_type = 'log10' if service == 'rapsearch2' else 'raw'
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = fetch_from_s3(self.additional_files["deuterostome_db"],
                                            self.ref_dir_local, allow_s3mi=True)
        m8.generate_taxon_count_json_from_m8(
            deduped_output_m8, output_hitsummary, evalue_type, db_type,
            lineage_db, deuterostome_db, output_counts_json)


    def run_remotely(self, input_fas, output_m8, service):
        key_path = self.fetch_key(os.environ['KEY_PATH_S3'])
        sample_name = self.output_dir_s3.rstrip('/').replace('s3://', '').replace('/', '-')
        chunk_size = self.additional_attributes["chunk_size"]
        index_dir_suffix = self.additional_attributes.get("index_dir_suffix")
        remote_username = "ec2-user"
        remote_home_dir = "/home/%s" % remote_username
        if service == "gsnap":
            remote_index_dir = "%s/share" % remote_home_dir
        elif service == "rapsearch2":
            remote_index_dir = "%s/references/nr_rapsearch" % remote_home_dir

        if index_dir_suffix:
            remote_index_dir = os.path.join(remote_index_dir, index_dir_suffix)

        remote_work_dir = "%s/batch-pipeline-workdir/%s" % (remote_home_dir, sample_name)

        # Split files into chunks for performance
        part_suffix, input_chunks = self.chunk_input(input_fas, chunk_size)

        # Process chunks
        chunk_output_files = [None] * len(input_chunks)
        chunk_threads = []
        mutex = threading.RLock()
        # Randomize execution order for performance
        randomized = list(enumerate(input_chunks))
        random.shuffle(randomized)

        for n, chunk_input_files in randomized:
            self.chunks_in_flight.acquire()
            self.check_for_errors(mutex, chunk_output_files, input_chunks, service)
            t = threading.Thread(
                target=PipelineStepRunAlignmentRemotely.run_chunk_wrapper,
                args=[
                    self.chunks_in_flight, chunk_output_files, n, mutex, self.run_chunk,
                    [
                        part_suffix, remote_home_dir, remote_index_dir,
                        remote_work_dir, remote_username, chunk_input_files,
                        key_path, service, True
                    ]
                ])
            t.start()
            chunk_threads.append(t)

        # Check chunk completion
        for ct in chunk_threads:
            ct.join()
            self.check_for_errors(mutex, chunk_output_files, input_chunks, service)
        assert None not in chunk_output_files
        # Concatenate the pieces and upload results
        self.concatenate_files(chunk_output_files, output_m8)

    def get_input_fas(self):
        service = self.additional_attributes["service"]
        if len(self.input_files_local[0]) == 1:
            return self.input_files_local[0]
        if len(self.input_files_local[0]) == 3:
            if service == 'gsnap':
                return self.input_files_local[0][0:2]
            if service == 'rapsearch2':
                return [self.input_files_local[0][2]]
        return None

    def fetch_key(self, key_path_s3):
        key_path = fetch_from_s3(key_path_s3, self.output_dir_local)
        command.execute("chmod 400 %s" % key_path)
        return key_path

    def chunk_input(self, input_files, chunksize):
        """Chunk input files into pieces for performance and parallelism."""
        part_lists = []  # Lists of partial files
        known_nlines = None
        part_suffix = ""
        chunk_nlines = chunksize * 2

        for input_file in input_files:
            # Count number of lines in the file
            nlines = int(command.execute_with_output("wc -l %s" % input_file)
                         .strip().split()[0])
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
            command.execute("split -a %d --numeric-suffixes -l %d %s %s" %
                            (ndigits, chunk_nlines, input_file, out_prefix))
            command.execute_with_retries(f"aws s3 sync --only-show-errors {self.chunks_result_dir_local}/ {self.chunks_result_dir_s3}/ --exclude '*' --include '{out_prefix_base}*'")

            # Get the partial file names
            partial_files = []
            paths = command.execute_with_output("ls %s*" % out_prefix).rstrip().split("\n")
            for pf in paths:
                partial_files.append(os.path.basename(pf))

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
        with open(output_m8, 'wb') as outf:
            for f in chunk_output_files:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, outf)

    @staticmethod
    def run_chunk_wrapper(chunks_in_flight, chunk_output_files, n, mutex, target, args):
        result = "error"
        try:
            result = target(*args)
        except:
            with mutex:
                traceback.print_exc()
        finally:
            with mutex:
                chunk_output_files[n] = result
            chunks_in_flight.release()

    def run_chunk(self, part_suffix, remote_home_dir, remote_index_dir, remote_work_dir,
                  remote_username, input_files, key_path, service, lazy_run):
        """Dispatch a chunk to worker machines for distributed GSNAP or RAPSearch
        group machines and handle their execution.
        """
        assert service in ("gsnap", "rapsearch2")

        chunk_id = input_files[0].split(part_suffix)[-1]
        multihit_basename = f"multihit-{service}-out{part_suffix}{chunk_id}.m8"
        multihit_local_outfile = os.path.join(self.chunks_result_dir_local, multihit_basename)
        multihit_remote_outfile = os.path.join(remote_work_dir, multihit_basename)
        multihit_s3_outfile = os.path.join(self.chunks_result_dir_s3, multihit_basename)

        base_str = "aws s3 cp --only-show-errors {s3_path}/{input_fa} {remote_work_dir}/{input_fa} "
        download_input_from_s3 = " ; ".join(
            base_str.format(
                s3_path=self.chunks_result_dir_s3,
                input_fa=input_fa,
                remote_work_dir=remote_work_dir) for input_fa in input_files)

        base_str = "mkdir -p {remote_work_dir} ; {download_input_from_s3} ; "
        environment = self.additional_attributes["environment"]
        if service == "gsnap":
            if environment == "staging":   # Pre-release testing of gsnap 2018-10-26.
                commands = base_str + "{remote_home_dir}/bin/gsnapl -A m8 --batch=0 --use-shared-memory=0 --gmap-mode=none --npaths=100 --ordered -t 36 --max-mismatches=40 -D {remote_index_dir} -d nt_k16 {remote_input_files} > {multihit_remote_outfile}"
            else:
                # max_search argument is still required in prod, until we upgrade to 2018-10-20 gsnap
                commands = base_str + "{remote_home_dir}/bin/gsnapl -A m8 --batch=0 --use-shared-memory=0 --gmap-mode=none --npaths=100 --ordered -t 36 --maxsearch=1000 --max-mismatches=40 -D {remote_index_dir} -d nt_k16 {remote_input_files} > {multihit_remote_outfile}"
        else:
            commands = base_str + "/usr/local/bin/rapsearch -d {remote_index_dir}/nr_rapsearch -e -6 -l 10 -a T -b 0 -v 50 -z 24 -q {remote_input_files} -o {multihit_remote_outfile}"

        commands = commands.format(
            remote_work_dir=remote_work_dir,
            download_input_from_s3=download_input_from_s3,
            remote_home_dir=remote_home_dir,
            remote_index_dir=remote_index_dir,
            remote_input_files=" ".join(
                remote_work_dir + "/" + input_fa for input_fa in input_files),
            multihit_remote_outfile=multihit_remote_outfile
            if service == "gsnap" else multihit_remote_outfile[:-3]
            # Strip the .m8 for RAPSearch as it adds that
        )

        if not lazy_run or not fetch_from_s3(multihit_s3_outfile,
                                             multihit_local_outfile):
            correct_number_of_output_columns = 12
            min_column_number = 0
            max_tries = 2
            try_number = 1
            instance_ip = ""


            def interpret_min_column_number_string(min_column_number_string,
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

            # Check if every row has correct number of columns (12) in the output
            # file on the remote machine
            while min_column_number != correct_number_of_output_columns \
                    and try_number <= max_tries:
                log.write("waiting for {} server for chunk {}".format(
                    service, chunk_id))
                max_concurrent = self.additional_attributes["max_concurrent"]

                with server.ASGInstance(service, key_path,
                                        remote_username, environment,
                                        max_concurrent, chunk_id,
                                        self.additional_attributes.get("max_interval_between_describe_instances") or 900,
                                        self.additional_attributes.get("job_tag_prefix") or "RunningIDseqBatchJob_",
                                        self.additional_attributes.get("job_tag_refresh_seconds") or 600,
                                        self.additional_attributes.get("draining_tag") or "draining") as instance_ip:
                    command.execute(command.remote(commands, key_path, remote_username, instance_ip))

                    if service == "gsnap":
                        verification_command = "cat %s" % multihit_remote_outfile
                    else:
                        # For rapsearch, first remove header lines starting with '#'
                        verification_command = "grep -v '^#' %s" % multihit_remote_outfile
                    verification_command += " | awk '{print NF}' | sort -nu | head -n 1"
                    min_column_number_string = command.execute_with_output(
                        command.remote(verification_command, key_path, remote_username, instance_ip))
                    min_column_number = interpret_min_column_number_string(
                        min_column_number_string, correct_number_of_output_columns,
                        try_number)

                try_number += 1

            # Move output from remote machine to local machine
            msg = "Chunk %s output corrupt; not copying to S3. Re-start pipeline " \
                  "to try again." % chunk_id
            assert min_column_number == correct_number_of_output_columns, msg

            with self.iostream_upload:  # Limit concurrent uploads so as not to stall the pipeline.
                command.execute(
                    command.scp(key_path, remote_username, instance_ip,
                                multihit_remote_outfile, multihit_local_outfile))
                command.execute("aws s3 cp --only-show-errors %s %s/" %
                                (multihit_local_outfile,
                                 self.chunks_result_dir_s3))
            log.write("finished alignment for chunk %s on %s server %s" % (chunk_id, service, instance_ip))
        return multihit_local_outfile
