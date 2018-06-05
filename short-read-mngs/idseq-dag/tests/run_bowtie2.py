import unittest

import os
import subprocess
import time

from idseq_dag.engine.pipeline_flow import PipelineFlow
from idseq_dag.steps.run_bowtie2 import PipelineStepRunBowtie2

class RunBowtie2Test(unittest.TestCase):

    def test_step_paired(self):
        input_dir_s3 = "s3://idseq-samples-prod/test_samples/1/results"
        output_dir_s3 = os.path.join(input_dir_s3, "testrun_bowtie2_paired_%d" % int(time.time()))
        result_dir_local ="/mnt/idseq/results/bowtie2_paired/%d" % os.getpid()
        input_files = [["lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta",
                        "lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta"]]
        output_files=["bowtie2.1.fa", "bowtie2.2.fa"]
        PipelineFlow.fetch_input_files_from_s3(input_files[0],
                                               input_dir_s3,
                                               result_dir_local)

        runstep = PipelineStepRunBowtie2(
            input_files=input_files,
            output_files=output_files,
            output_dir_local=result_dir_local,
            output_dir_s3=output_dir_s3,
            ref_dir_local='/mnt/idseq/ref',
            additional_files={},
            additional_attributes={}
        )
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder

    def test_step_single(self):
        input_dir_s3 = "s3://idseq-samples-prod/test_samples/2/results"
        output_dir_s3 = os.path.join(input_dir_s3, "testrun_bowtie2_single_%d" % int(time.time()))
        result_dir_local ="/mnt/idseq/results/bowtie2-single/%d" % os.getpid()
        input_files = [["lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta"]]
        output_files=["bowtie2.1.fa"]
        PipelineFlow.fetch_input_files_from_s3(input_files[0],
                                               input_dir_s3,
                                               result_dir_local)
        runstep = PipelineStepRunBowtie2(
            input_files=input_files,
            output_files=output_files,
            output_dir_local=result_dir_local,
            output_dir_s3=output_dir_s3,
            ref_dir_local='/mnt/idseq/ref',
            additional_files={},
            additional_attributes={}
        )

        runstep.start()
        runstep.wait_until_finished()


