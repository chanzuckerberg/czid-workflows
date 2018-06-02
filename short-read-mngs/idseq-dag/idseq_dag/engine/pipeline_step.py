import json
import sys
import os
import threading
import subprocess
import time

class PipelineStep:
    def __init__(self, input_files, output_files,
                 output_dir_local, output_dir_s3, ref_dir_local,
                 additional_files, additional_attributes):
      ''' Set up all the input_files and output_files here '''
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
        return [os.path.join(self.output_dir_local,f) for f in self.output_files]

    def start_uploading_results(self):
        for f in self.output_files_local():
            # Tag the done files
            done_file = self.done_file(f)
            subprocess.check_call("date > %s" % done_file, shell=True)
        for f in self.output_files_local():
            # upload to S3 - TODO(Boris): parallelize the following with better calls
            subprocess.check_call("aws s3 cp %s %s/" % (f, self.output_dir_s3), shell=True)

    @staticmethod
    def done_file(filename):
        return "%s.done" % filename

    def wait_for_file_existence(self, filename):
        while True:
            if os.path.exists(filename):
                return
            time.sleep(5) # wait for 5 seconds to check again

    def wait_for_input_files(self):
        for fl in self.input_files:
            flist = []
            for f in fl:
                local_file = os.path.join(self.output_dir_local,f)
                self.wait_for_file_existence(local_file)
                self.wait_for_file_existence(self.done_file(local_file))
                flist.append(f)
            self.input_files_local.append(flist)


    def save_progress(self):
        # save stats
        # save other info
        # TO BE IMPLEMENTED
        pass

    def thread_run(self):
        #Timer.start()
        self.started = True
        self.wait_for_input_files()
        self.run()
        self.finished = True
        self.save_progress()
        #Timer.finalize()
        self.start_uploading_results()
        self.all_done = True

    def start(self):
        ''' Something like the following. It should be non blocking '''
        self.thread = threading.Thread(target=self.thread_run)
        self.thread.start()


