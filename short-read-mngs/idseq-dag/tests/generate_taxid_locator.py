import unittest

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.steps.generate_taxid_locator import PipelineStepGenerateTaxidLocator


class RunGenerateTaxidLocatorTest(unittest.TestCase):
    @staticmethod
    def test_step_single():
        run_step = IdseqStepSetup.get_step_object(PipelineStepGenerateTaxidLocator, "taxid_locator_out", paired=False)
        run_step.start()
        run_step.wait_until_finished()
        # TODO: Check results
        # TODO: Clean up the folder
