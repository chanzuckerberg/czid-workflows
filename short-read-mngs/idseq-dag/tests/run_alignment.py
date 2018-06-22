import unittest

import os
import subprocess
import time

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.engine.pipeline_flow import PipelineFlow
from idseq_dag.steps.run_alignment_remotely import PipelineStepRunAlignmentRemotely

class RunAlignmentTest(unittest.TestCase):

    def test_step_gsnap_paired(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunAlignmentRemotely, "gsnap_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # TODO: Check results
        # Clean up the folder

    def test_step_rapsearch2_paired(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunAlignmentRemotely, "rapsearch2_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # TODO: Check results
        # Clean up the folder

    def test_step_gsnap_single(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunAlignmentRemotely, "gsnap_out", paired=False)
        runstep.start()
        runstep.wait_until_finished()
        # TODO: Check results
        # Clean up the folder

    def test_step_rapsearch2_single(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunAlignmentRemotely, "rapsearch2_out", paired=False)
        runstep.start()
        runstep.wait_until_finished()
        # TODO: Check results
        # Clean up the folder
