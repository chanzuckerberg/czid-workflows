import os
import time
import unittest

from tests.unit.unittest_helpers import relative_file_path

from idseq_dag.steps.generate_lz4 import PipelineStepGenerateLZ4

INPUT_FILE = relative_file_path(__file__, 'doesnotexist')

class TestPipelineStepGenerateLZ4(unittest.TestCase):

    def setUp(self):
        self.step = PipelineStepGenerateLZ4(
            name='test_generate_lz4',
            input_files=[[INPUT_FILE]],
            output_files=[],
            output_dir_local='',
            ref_dir_local='',
            output_dir_s3='',
            additional_files={},
            additional_attributes={},
        )

    def test_get_command(self):
        command = self.step.get_command(INPUT_FILE)
        self.assertEqual('lz4', command.cmd)
        self.assertSequenceEqual(
            ['-9', '-f', INPUT_FILE, INPUT_FILE + '.lz4'],
            command.args
        )
