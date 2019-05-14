import unittest
from unittest.mock import patch, call

# Module under test
import idseq_dag.util.log as log

class TestLog(unittest.TestCase):
    '''Tests for `util/log.py`'''

    @staticmethod
    @patch('idseq_dag.util.log.write')
    def test_log_event(mock_log_write):
        '''Test a successful log event'''
        log.log_event('fake_event_name', values={"name_a": 1, "name_b": "fake_string", "name_c": 3.1415}, warning=True, flush=False, extra_fields={"extra_field_1": "abc"})

        mock_log_write.assert_has_calls([
            call('{"event": "fake_event_name", "extra_field_1": "abc", "values": {"name_a": 1, "name_b": "fake_string", "name_c": 3.1415}}', warning=True, flush=False)
        ])

    @staticmethod
    @patch('idseq_dag.util.log.write')
    def test_log_event_2(mock_log_write):
        '''Test log event with a non serializable field'''
        invalid_field = bytes()

        log.log_event('fake_event_name', values={"name_a": invalid_field, "name_b": "fake_string", "name_c": 3.1415})

        mock_log_write.assert_has_calls([
            call('{"event": "fake_event_name", "values": {"name_a": "<<non-serializable: bytes>>", "name_b": "fake_string", "name_c": 3.1415}}', warning=False, flush=True)
        ])
