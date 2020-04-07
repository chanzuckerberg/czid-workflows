import subprocess
import re
import unittest
from unittest.mock import patch, call, ANY
import itertools

from tests.unit.unittest_helpers import MATCH_RE

# Class under test
from idseq_dag.steps.run_alignment_remotely import PipelineStepRunAlignmentRemotely, CORRECT_NUMBER_OF_OUTPUT_COLUMNS, CHUNK_DONT_RETRY_AFTER_SECONDS

FAKE_BAD_IP = "98.76.54.32"
FAKE_GOOD_IP = "12.34.56.78"
FAKE_RUN_CHUNK_ARGS = {"part_suffix": "chunksize-20000-numparts-5-part-",
                       "remote_home_dir": "/tmp/home",
                       "remote_index_dir": "/tmp/index",
                       "remote_work_dir": "/home/ec2-user/batch-pipeline-workdir/idseq-samples-TEST-samples-150-18503-results-3.7",
                       "remote_username": "dummy",
                       "input_files": ["multihit-rapsearch2-out-chunksize-20000-numparts-5-part-3"],
                       "key_path": "/tmp/key_path",
                       "service": "rapsearch2",
                       "lazy_run": False}

def _side_effect_command_execute(cmd):
    if FAKE_BAD_IP in cmd:
        raise subprocess.CalledProcessError(255, cmd)
    return 0


class TestRunAlignmentRemotely(unittest.TestCase):
    '''Tests for `util/trace_rlock.py`'''

    def setUp(self):
        self.step = PipelineStepRunAlignmentRemotely(name="something_out",
                                                     input_files=["host_filter_out"],
                                                     output_files=["rapsearch2.m8", "rapsearch2.deduped.m8",
                                                                   "rapsearch2.hitsummary.tab", "rapsearch2_counts.json"],
                                                     output_dir_local="/tmp/output_local_dir/",
                                                     output_dir_s3="s3://idseq-samples-TEST/samples/150/18503/results/3.7/",
                                                     ref_dir_local="/tmp/ref_dir_local/",
                                                     additional_files=[],
                                                     additional_attributes={
                                                         "service": "rapsearch2",
                                                         "chunks_in_flight": 32,
                                                         "chunk_size": 20000,
                                                         "max_concurrent": 4,
                                                         "environment": "prod",
                                                         "max_interval_between_describe_instances": 900,
                                                         "job_tag_prefix": "RunningIDseqBatchJob_",
                                                         "job_tag_refresh_seconds": 600,
                                                         "draining_tag": "draining",
                                                         "index_dir_suffix": "2018-12-01"
                                                     })


    @patch('idseq_dag.steps.run_alignment_remotely.server.ASGInstance')
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute', side_effect=_side_effect_command_execute)
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute_with_output', return_value=CORRECT_NUMBER_OF_OUTPUT_COLUMNS)
    def test_run_chunk_success(self, mock_command_execute_with_output, mock_command_execute, mock_asg_instance):
        '''Run command remotely with success'''
        mock_asg_instance.return_value.__enter__.return_value = FAKE_GOOD_IP

        self.step.run_chunk(**FAKE_RUN_CHUNK_ARGS)

        mock_command_execute.assert_has_calls([
            call(MATCH_RE(f"^ssh.*{re.escape(FAKE_GOOD_IP)}")),
            call(MATCH_RE(f"^scp.*{re.escape(FAKE_GOOD_IP)}")),
            call(MATCH_RE(f"^aws s3 cp"))
        ])
        mock_command_execute_with_output.assert_called_once_with(MATCH_RE('^ssh.*grep'))


    @patch('idseq_dag.steps.run_alignment_remotely.server.ASGInstance')
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute', side_effect=_side_effect_command_execute)
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute_with_output', return_value=CORRECT_NUMBER_OF_OUTPUT_COLUMNS)
    @patch('idseq_dag.steps.run_alignment_remotely.PipelineStepRunAlignmentRemotely._exponential_backoff')
    def test_run_chunk_with_retry(self, mock_exponential_backoff, mock_command_execute_with_output, mock_command_execute, mock_asg_instance):
        '''Run command remotely succeeded after a retry'''
        mock_asg_instance.return_value.__enter__.side_effect = [FAKE_BAD_IP, FAKE_GOOD_IP]

        self.step.run_chunk(**FAKE_RUN_CHUNK_ARGS)

        mock_command_execute.assert_has_calls([
            call(MATCH_RE(f"^ssh.*{re.escape(FAKE_GOOD_IP)}")),
            call(MATCH_RE(f"^scp.*{re.escape(FAKE_GOOD_IP)}")),
            call(MATCH_RE(f"^aws s3 cp"))
        ])
        mock_command_execute_with_output.assert_called_once_with(MATCH_RE('^ssh.*grep'))
        mock_exponential_backoff.assert_called()


    @patch('idseq_dag.steps.run_alignment_remotely.server.ASGInstance')
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute', side_effect=_side_effect_command_execute)
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute_with_output', return_value=CORRECT_NUMBER_OF_OUTPUT_COLUMNS)
    @patch('idseq_dag.steps.run_alignment_remotely.PipelineStepRunAlignmentRemotely._exponential_backoff')
    def test_run_chunk_fail_after_retries(self, mock_exponential_backoff, mock_command_execute_with_output, mock_command_execute, mock_asg_instance):
        '''Run command remotely failed after retries due to ssh failures'''
        mock_asg_instance.return_value.__enter__.side_effect = itertools.repeat(FAKE_BAD_IP)

        with self.assertRaises(subprocess.CalledProcessError):
            self.step.run_chunk(**FAKE_RUN_CHUNK_ARGS)

        mock_command_execute.assert_has_calls([
            call(MATCH_RE(f"^ssh.*{re.escape(FAKE_BAD_IP)}"))
        ])
        mock_command_execute_with_output.assert_not_called()
        mock_exponential_backoff.assert_called()


    @patch('idseq_dag.steps.run_alignment_remotely.server.ASGInstance')
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute', side_effect=_side_effect_command_execute)
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute_with_output', side_effect=itertools.repeat(-1))
    @patch('idseq_dag.steps.run_alignment_remotely.PipelineStepRunAlignmentRemotely._exponential_backoff')
    def test_run_chunk_fail_chunk_size(self, mock_exponential_backoff, mock_command_execute_with_output, mock_command_execute, mock_asg_instance):
        '''Run command remotely failed after retries due to chunk assertion errors'''
        mock_asg_instance.return_value.__enter__.return_value = FAKE_GOOD_IP

        with self.assertRaises(AssertionError):
            self.step.run_chunk(**FAKE_RUN_CHUNK_ARGS)

        mock_command_execute.assert_has_calls([
            call(MATCH_RE(f"^ssh.*{re.escape(FAKE_GOOD_IP)}"))
        ])
        mock_command_execute_with_output.assert_has_calls([
            call(MATCH_RE(f"^ssh.*{re.escape(FAKE_GOOD_IP)}"))
        ])
        mock_exponential_backoff.assert_called()


    @patch('idseq_dag.steps.run_alignment_remotely.server.ASGInstance')
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute', side_effect=_side_effect_command_execute)
    @patch('idseq_dag.steps.run_alignment_remotely.command.execute_with_output', return_value=CORRECT_NUMBER_OF_OUTPUT_COLUMNS)
    @patch('idseq_dag.steps.run_alignment_remotely.PipelineStepRunAlignmentRemotely._exponential_backoff')
    @patch('idseq_dag.steps.run_alignment_remotely.PipelineStepRunAlignmentRemotely._elapsed_seconds', return_value=CHUNK_DONT_RETRY_AFTER_SECONDS+100)
    def test_run_chunk_do_not_retry(self, _mock_elapsed_seconds, mock_exponential_backoff, _mock_command_execute_with_output, _mock_command_execute, mock_asg_instance):
        '''Run command remotely failed with no retries because it took longer than CHUNK_DONT_RETRY_AFTER_SECONDS threshold'''
        mock_asg_instance.return_value.__enter__.return_value = FAKE_BAD_IP

        with self.assertRaises(subprocess.CalledProcessError):
            self.step.run_chunk(**FAKE_RUN_CHUNK_ARGS)

        mock_exponential_backoff.assert_not_called()
