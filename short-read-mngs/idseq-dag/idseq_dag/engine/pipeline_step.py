import json
import sys
import os
import threading

class PipelineStep:
    def __init__(self, input_files, output_files,
                 output_dir_local, output_dir_s3, ref_dir_local,
                 additional_files, additional_attributes, ref_dir_local):
      ''' Set up all the input_files and output_files here '''
      self.finished = False
      self.started = False
      self.input_files = input_files # local location
      self.output_files = output_files # s3 location
      self.additional_files = additional_files
      self.additional_attributes = additional_attributes
      self.thread = None

    @abstractmethod
    def run(self):
      ''' implement what is actually being run in this step '''
    def finished(self):
       ''' Finished Running '''
    def all_done(self):
       ''' Finished running and uploading results '''
    def wait_until_finished(self):
       ''' return true when the run is done (before results are uploaded) '''

    def thread_start(self):
        #Timer.start()
        self.ready_to_run()
        self.run()
        self.save_progress_to_s3()
        Timer.finalize()
        #self.start_uploading_results()

    def start(self):
        ''' Something like the following. It should be non blocking '''
        self.thread = threading.Thread(target=self.thread_start)
        self.thread.start()


