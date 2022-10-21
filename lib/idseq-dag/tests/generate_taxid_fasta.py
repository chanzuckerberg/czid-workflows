import unittest

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.steps.generate_taxid_fasta import PipelineStepGenerateTaxidFasta


class GenerateTaxidFastaTest(unittest.TestCase):
    def test_step_single(self):
        run_step = IdseqStepSetup.get_step_object(PipelineStepGenerateTaxidFasta, "taxid_fasta_out", paired=False)
        run_step.start()
        run_step.wait_until_finished()
        # TODO: Check results
        # TODO: Clean up the folder
