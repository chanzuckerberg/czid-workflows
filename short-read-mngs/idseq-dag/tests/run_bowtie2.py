import unittest

import os
import subprocess
import time

from idseq_dag.steps.run_bowtie2 import PipelineStepRunBowtie2

class RunBowtie2Test(unittest.TestCase):
    @staticmethod
    def setup_input_files(input_dir_s3, result_dir_local, input_files):
        # We might be able to implement this just in the PipelineFlow
        for fl in self.input_files:
            for f in fl:
                local_file = os.path.join(result_dir_local, f)
                done_file = PipelineStep.done_file(local_file)
                subprocess.check_call("aws s3 cp %s/%s %s/" % (input_dir_s3, f, result_dir_local), shell=True)
                subprocess.check_call("date > %s" % done_file)


    def test_step_paired(self):
        input_dir_s3 = "s3://idseq-samples-prod/tests_samples/1/results"
        output_dir_s3 = os.path.join(input_dir_s3, "testrun_bowtie2_paired_%d" % int(time.time()))
        result_dir_local ="/mnt/idseq/results/bowtie2_paired/%d" % os.getpid()
        input_files = [["lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta",
                        "lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta"]]
        output_files=[["bowtie2.1.fa", "bowtie2.2.fa"]]
        self.setup_input_files(input_dir_s3, result_dir_local, input_files)

        self.runstep = PipelineStepRunBowtie2(
            input_files=input_files,
            output_files=output_files,
            output_dir_local=result_dir_local,
            output_dir_s3=output_dir_s3,
            ref_dir_local='/mnt/idseq/ref',
            additional_files={},
            additional_attributes={}
        )
        self.runstep.start()
        self.runstep.wait_until_succeed()
        # Check results
        # Clean up the folder

    def test_step_single(self):
        input_dir_s3 = "s3://idseq-samples-prod/tests_samples/2/results"
        output_dir_s3 = os.path.join(input_dir_s3, "testrun_bowtie2_single_%d" % int(time.time()))
        result_dir_local ="/mnt/idseq/results/bowtie2-single/%d" % os.getpid()
        input_files = [["lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta"]]
        output_files=[["bowtie2.1.fa"]]
        self.setup_input_files(input_dir_s3, result_dir_local, input_files)
        self.runstep.start()


