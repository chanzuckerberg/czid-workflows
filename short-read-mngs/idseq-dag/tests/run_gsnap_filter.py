import unittest

from idseq_dag.steps.run_gsnap_filter import PipelineStepRunGsnapFilter
from .idseq_step_setup import IdseqStepSetup

class RunGsnapFilterTest(unittest.TestCase):

    def test_step_paired(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunGsnapFilter,
                                                 "gsnap_filter_out",
                                                 paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # TODO
        # Check results
        # Clean up the folder

    def test_step_single(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunGsnapFilter,
                                                 "gsnap_filter_out",
                                                 paired=False)
        runstep.start()
        runstep.wait_until_finished()
        # TODO
        # Check results
        # Clean up the folder
