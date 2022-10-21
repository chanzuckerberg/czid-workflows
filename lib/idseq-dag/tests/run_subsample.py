import unittest

from idseq_dag.steps.run_subsample import PipelineStepRunSubsample
from .idseq_step_setup import IdseqStepSetup

class RunSubsampleTest(unittest.TestCase):

    def test_step_paired(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunSubsample,
                                                 "subsampled_out",
                                                 paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # TODO
        # Check results
        # Clean up the folder

    def test_step_single(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunSubsample,
                                                 "subsampled_out",
                                                 paired=False)
        runstep.start()
        runstep.wait_until_finished()
        # TODO
        # Check results
        # Clean up the folder
