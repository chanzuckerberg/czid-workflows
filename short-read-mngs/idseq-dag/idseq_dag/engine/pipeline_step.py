import json
import sys
import os
import multiprocess

class PipelineStep:
    def __init__(self, lazy_run,
                 input_files_local, output_files_s3, output_dir_s3,
                 additional_files, additional_attributes):
      ''' Set up all the input_files and output_files here '''
      self.finished = False
      self.started = False
      self.lazy_run = lazy_run
      self.input_files = input_files # local location
      self.output_files = output_files # s3 location
      self.additional_files = additional_files
      self.additional_attributes = additional_attributes

    @abstractmethod
    def run(self):
      ''' implement what is actually being run in this step '''
    def need_to_run(self):
       ''' check if output_file exists, signature and stuff '''
       if not lazy_run:
           return True
    def start_downloader(self):
       ''' will start downloading the heavy files '''
    def finished(self):
       ''' Finished Running '''
    def all_done(self):
       ''' Finished running and uploading results '''
    def wait_until_finished(self):
       ''' return true when the run is done (before results are uploaded) '''

    def start(self):
        ''' Something like the following. It should be non blocking '''
        return unless self.need_to_run() or not self.started
        Timer.start()
        self.ready_to_run()
        self.run()
        self.save_progress_to_s3()
        Timer.finalize()
        self.start_uploading_results()


