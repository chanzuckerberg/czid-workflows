import json
import sys
import os
import threading
import time
import idseq_dag.util.command as command
from abc import abstractmethod

class PipelineStep(object):
    ''' Each Pipeline Run Step i.e. run_star, run_bowtie, etc '''
    def __init__(self, name, input_files, output_files,
                 output_dir_local, output_dir_s3, ref_dir_local,
                 additional_files, additional_attributes):
        ''' Set up all the input_files and output_files here '''
        self.name = name
        self.input_files = input_files # list of list files
        self.output_files = output_files # s3 location
        self.output_dir_local = output_dir_local
        self.output_dir_s3 = output_dir_s3.rstrip('/')
        self.ref_dir_local = ref_dir_local
        self.additional_files = additional_files
        self.additional_attributes = additional_attributes

        self.started = False
        self.finished = False # done running
        self.all_done = False # done running and uploading
        self.thread = None
        self.input_files_local = []

    @abstractmethod
    def run(self):
        ''' implement what is actually being run in this step '''

    def output_files_local(self):
        ''' Get list of output files on local folder '''
        return [os.path.join(self.output_dir_local, f) for f in self.output_files]

    def start_uploading_results(self):
        ''' Upload output files to s3 '''
        for f in self.output_files_local():
            # upload to S3 - TODO(Boris): parallelize the following with better calls
            command.execute("aws s3 cp %s %s/" % (f, self.output_dir_s3))

    @staticmethod
    def done_file(filename):
        ''' get the done file for a particular local file '''
        return "%s.done" % filename

    @staticmethod
    def wait_for_file_existence(filename):
        ''' wait until file exists '''
        # TODO: (Yunfang) What if the file will never appear because something else failed?
        # It is customary in those cases for an exception to be raised in this loop.
        while True:
            if os.path.exists(filename):
                return
            time.sleep(5) # wait for 5 seconds to check again

    def wait_for_input_files(self):
        ''' wait for all the input files to be avaiable and update input_files_local '''
        for fl in self.input_files:
            flist = []
            for f in fl:
                local_file = os.path.join(self.output_dir_local, f)
                self.wait_for_file_existence(local_file)
                self.wait_for_file_existence(self.done_file(local_file))
                flist.append(local_file)
            self.input_files_local.append(flist)


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
            command.execute("date > %s" % done_file)

    def wait_until_finished(self):
        # TODO(yf): use self.join
        while True:
            if self.finished:
                return
            time.sleep(5)

    def wait_until_all_done(self):
        while True:
            if self.all_done:
                return
            time.sleep(5)

    def thread_run(self):
        ''' Actually running the step '''
        #Timer.start()
        self.started = True
        self.wait_for_input_files()
        self.run()
        self.validate()
        self.finished = True
        self.save_progress()
        #Timer.finalize()
        # TODO(yf): move the uploading to a diffdrent thread
        self.start_uploading_results()
        self.all_done = True

    def start(self):
        ''' function to be called after instanitation to start running the step '''
        self.thread = threading.Thread(target=self.thread_run)
        self.thread.start()
