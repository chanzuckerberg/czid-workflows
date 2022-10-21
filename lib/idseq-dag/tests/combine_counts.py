import unittest

import os
import subprocess
import time

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.engine.pipeline_flow import PipelineFlow
from idseq_dag.steps.combine_taxon_counts import PipelineStepCombineTaxonCounts

class RunCombineCountsTest(unittest.TestCase):

    def test_step(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepCombineTaxonCounts, "taxon_count_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder

