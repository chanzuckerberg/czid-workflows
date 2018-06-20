import unittest

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.steps.generate_alignment_viz import PipelineStepGenerateAlignmentViz


class GenerateAlignmentVizTest(unittest.TestCase):
    def test_step_single(self):
        run_step = IdseqStepSetup.get_step_object(PipelineStepGenerateAlignmentViz, "alignment_viz_out", paired=False)
        run_step.start()
        run_step.wait_until_finished()
        # TODO: Check results
        # TODO: Clean up the folder
