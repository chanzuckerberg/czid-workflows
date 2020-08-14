import unittest

import os
import subprocess
import time

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.engine.pipeline_flow import PipelineFlow
from idseq_dag.steps.generate_annotated_fasta import PipelineStepGenerateAnnotatedFasta

class RunAnnotatedFastaTest(unittest.TestCase):

    def test_step(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepGenerateAnnotatedFasta, "annotated_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder

