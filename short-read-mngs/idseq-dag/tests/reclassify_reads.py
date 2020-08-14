import unittest

import os
import subprocess
import time

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.engine.pipeline_flow import PipelineFlow
from idseq_dag.steps.reclassify_reads import PipelineStepReclassifyReads

class ReclassifyReadsTest(unittest.TestCase):

    def test_step_gsnap(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepReclassifyReads, "refined_gsnap_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder

    def test_step_rapsearch2(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepReclassifyReads, "refined_rapsearch2_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder



