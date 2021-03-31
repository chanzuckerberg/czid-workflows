import datetime
import json
import os
import threading
import time
import traceback
from abc import abstractmethod
from enum import Enum, IntEnum

import pytz

import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.s3
import idseq_dag.util.count as count

from idseq_dag.util.count import load_duplicate_cluster_sizes
from idseq_dag.exceptions import InvalidInputFileError, InvalidOutputFileError


class StepStatus(IntEnum):
    INITIALIZED = 0
    STARTED = 1  # step.start() called
    FINISHED = 2  # step.run() finished
    UPLOADED = 3  # all results uploaded to s3
    INVALID_INPUT = 4  # an error occurred when validating the input file


class PipelineStep(object):
    ''' Each Pipeline Run Step i.e. run_star, run_bowtie2, etc '''

    def __init__(self, name, input_files, output_files,
                 output_dir_local, output_dir_s3, ref_dir_local,
                 additional_files, additional_attributes):
        ''' Set up all the input_files and output_files here '''
        self.name = name
        self.input_files = input_files  # list of list files
        self.output_files = output_files  # s3 location
        self.output_dir_local = output_dir_local
        self.output_dir_s3 = output_dir_s3.rstrip('/')
        self.ref_dir_local = ref_dir_local
        self.create_local_dirs()

        self.additional_files = additional_files
        self.additional_attributes = additional_attributes

        self.status = StepStatus.INITIALIZED
        self.exec_thread = None
        self.upload_thread = None
        self.input_files_local = []

        # Extra output files for internal use, not visible to users
        #   Automatically uploaded to s3
        self.additional_output_files_hidden = []
        # Extra output files available to users
        #   Automatically uploaded to s3 and
        #   added to status_dict with key "additional_output"
        self.additional_output_files_visible = []
        # Extra output folders, not visible to users
        #   Currently all extra output folders are hidden
        #   Automatically uploaded to s3
        self.additional_output_folders_hidden = []

        # This exists only to check for filename collisions.
        self._files_seen = set()

        self.counts_dict = {}
        self.should_terminate = False
        self.should_count_reads = False

        self.input_file_error = None
        self.upload_results_with_checksum = False

    @abstractmethod
    def run(self):
        ''' implement what is actually being run in this step '''

    @abstractmethod
    def count_reads(self):
        ''' count reads '''

    def stop_waiting(self):
        ''' stop waiting to run '''
        self.should_terminate = True

    def save_counts(self):
        if self.counts_dict:
            count_file_name = "%s/%s.count" % (self.output_dir_local, self.name)
            with open(count_file_name, 'w') as count_file:
                json.dump(self.counts_dict, count_file)
            self.additional_output_files_hidden.append(count_file_name)

    def output_files_local(self):
        ''' Get list of output files on local folder '''
        return [os.path.join(self.output_dir_local, f) for f in self.output_files]

    def create_local_dirs(self):
        ''' make sure proper local directories are created for files with subdirs '''
        for f in self.output_files_local():
            command.make_dirs(os.path.dirname(f))

    def uploading_results(self):
        ''' Upload output files to s3 '''
        files_to_upload = self.output_files_local() + self.additional_output_files_hidden + self.additional_output_files_visible
        for f in files_to_upload:
            if f in self._files_seen:
                raise InvalidOutputFileError({
                    "error": "Filename conflict: {} already exists in {}".format(f, self._files_seen),
                    "step": self.name
                })
            else:
                self._files_seen.add(f)

            # upload to S3 - TODO(Boris): parallelize the following with better calls
            s3_path = self.s3_path(f)
            idseq_dag.util.s3.upload_with_retries(f, s3_path, checksum=self.upload_results_with_checksum)
        for f in self.additional_output_folders_hidden:
            idseq_dag.util.s3.upload_folder_with_retries(f, self.s3_path(f), checksum=self.upload_results_with_checksum)
        self.status = StepStatus.UPLOADED

    @staticmethod
    def done_file(filename):
        ''' get the done file for a particular local file '''
        return "%s.done" % filename

    def wait_for_input_files(self):
        ''' wait for all the input files to be available and update input_files_local '''
        for fl in self.input_files:
            flist = []
            for f in fl:
                if f in self._files_seen:
                    raise InvalidInputFileError({
                        "error": "Filename conflict: {} already exists in {}".format(f, self._files_seen),
                        "step": self.name
                    })
                else:
                    self._files_seen.add(f)

                local_file = os.path.join(self.output_dir_local, f)
                while True:
                    if os.path.exists(local_file) and os.path.exists(self.done_file(local_file)):
                        flist.append(local_file)
                        break
                    else:
                        if self.should_terminate:
                            # If the step is not supposed to be run any more.
                            raise RuntimeError("Step %s being terminated" % self.name)
                        time.sleep(5)
            self.input_files_local.append(flist)

    def validate_input_files(self):
        ''' Validate input files before running the step.
        Should assign any error encountered to self.input_file_error
        '''
        pass

    def save_progress(self):
        ''' save progress after step run '''
        # save stats
        # save other info
        # TO BE IMPLEMENTED
        pass

    def validate(self):
        ''' Make sure all the output files are generated. '''
        for f in self.output_files_local():
            if not os.path.exists(f):
                raise RuntimeError("output file %s should be generated after run" % f)
            # Tag the done files
            done_file = self.done_file(f)
            fmt_now = datetime.datetime.now(tz=pytz.UTC).strftime("%a %b %e %H:%M:%S %Z %Y")
            command.write_text_to_file(fmt_now, done_file)
        self.count_reads()

    def wait_until_finished(self):
        self.exec_thread.join()
        if self.status == StepStatus.INVALID_INPUT:
            raise InvalidInputFileError({
                "error": self.input_file_error.name,
                "step": self.name
            })

        if self.status < StepStatus.FINISHED:
            raise RuntimeError("step %s run failed" % self.name)

    def wait_until_all_done(self):
        self.wait_until_finished()
        # run finished
        self.upload_thread.join()

        if self.status < StepStatus.UPLOADED:
            raise RuntimeError("step %s uploading failed" % self.name)

    def thread_run(self):
        ''' Actually running the step '''
        self.status = StepStatus.STARTED

        v = {"step": self.name}
        with log.log_context("dag_step", v):
            with log.log_context("substep_wait_for_input_files", v):
                self.wait_for_input_files()
            with log.log_context("substep_validate_input_files", v):
                self.validate_input_files()

            # If an input file error was detected, stop execution.
            if self.input_file_error:
                log.write("Invalid input detected for step %s" % self.name)
                self.status = StepStatus.INVALID_INPUT
                return

            with log.log_context("substep_run", v):
                self.run()
            with log.log_context("substep_validate", v):
                self.validate()
            with log.log_context("substep_save_progress", v):
                self.save_progress()
            with log.log_context("substep_save_counts", v):
                self.save_counts()
        self.upload_thread = threading.Thread(target=self.uploading_results)
        self.upload_thread.start()
        self.status = StepStatus.FINISHED

    def start(self):
        ''' function to be called after instantiation to start running the step '''
        self.exec_thread = threading.Thread(target=self.thread_run)
        self.exec_thread.start()

    def relative_path(self, local_path):
        return os.path.relpath(local_path, self.output_dir_local)

    def s3_path(self, local_path):
        return os.path.join(self.output_dir_s3, self.relative_path(local_path))

    def step_description(self, require_docstrings=False):
        ''' Retrieves description for the given step.
        By default, it pulls the docstring of the class but
        should be overridden for more dynamic descriptions
        that depends on the inputs. If no docstring is provided,
        it throws an exception.
        '''
        docstring = self.__doc__ or ""
        if not docstring and require_docstrings:
            raise TypeError(f"No docstring for step {self.name}")
        return docstring.strip()

    def step_resources(self):
        ''' Returns a dictionary of resources in the form of display name => url.

        These will be used on the sidebar of the pipeline visualization on the idseq-web app.
        By default, show link to idseq-dag documentation.
        '''
        return {"IDseq Docs": "https://github.com/chanzuckerberg/idseq-dag/wiki"}


class PipelineCountingStep(PipelineStep):

    """PipelineStep that counts nonunique reads based on back-calculation from cluster sizes
TSV file emitted by `PipelineStepRunIDSeqDedup`. Only steps that follow idseq-dedup are eligible
for this, and not all of them (not all steps count their outputs)."""

    def input_cluster_sizes_path(self):
        # The last last input to PipelineCountingStep is cluster_sizes.tsv
        tsv = self.input_files_local[-1][-1]
        assert tsv.endswith(".tsv"), str(self.input_files_local)
        return tsv

    def _count_reads_work(self, cluster_key, counter_name, fasta_files):
        # Count reads including duplicates (expanding duplicate clusters).
        self.should_count_reads = True
        self.counts_dict[counter_name] = count.reads_in_group(
            file_group=fasta_files,
            cluster_sizes=load_duplicate_cluster_sizes(self.input_cluster_sizes_path()),
            cluster_key=cluster_key)

    def count_reads(self):
        # Steps which decorate their read IDs must override this function to specify
        # a cluster_key that inverses the read ID decorator so that the original
        # cluster sizes can be looked up.
        self._count_reads_work(
            cluster_key=lambda x: x,
            counter_name=self.name,
            fasta_files=self.output_files_local()[0:2]
        )

    @abstractmethod
    def run(self):
        ''' implement what is actually being run in this step '''
