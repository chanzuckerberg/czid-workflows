import unittest

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.steps.run_star import PipelineStepRunStar


class RunStarTest(unittest.TestCase):
    def test_step_paired(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunStar, "star_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # TODO: Check results
        # TODO: Clean up the folder

    def test_step_single(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunStar, "star_out", paired=False)
        runstep.start()
        runstep.wait_until_finished()
        # TODO: Check results
        # TODO: Clean up the folder
