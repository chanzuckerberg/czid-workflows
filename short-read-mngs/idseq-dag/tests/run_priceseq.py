import unittest
from .idseq_step_setup import IdseqStepSetup
from idseq_dag.steps.run_priceseq import PipelineStepRunPriceSeq


class RunPriceSeqTest(unittest.TestCase):
    def test_step_paired(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunPriceSeq, "priceseq_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder

    def test_step_single(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunPriceSeq, "priceseq_out", paired=False)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder
