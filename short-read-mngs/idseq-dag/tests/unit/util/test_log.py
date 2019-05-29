import datetime
import io
import logging
import unittest
from unittest.mock import patch, call, Mock

# Module under test
import idseq_dag.util.log as log

class TestLog(unittest.TestCase):
    '''Tests for `util/log.py`'''

    @patch('logging.getLogger')
    def test_write_debug_text(self, mock_getLogger):
        '''Test writing a info text message to the logs'''

        log.write("Text message", debug=True)

        mock_getLogger.assert_has_calls([
            call().debug(msg='Text message', extra={"obj_data": None})
        ])


    @patch('logging.getLogger')
    def test_write_info_obj_data(self, mock_getLogger):
        '''Test writing a info object data to the logs'''

        log.write(obj_data={"name_a": 1, "name_b": "fake_string", "name_c": 3.1415})

        mock_getLogger.assert_has_calls([
            call().info(msg=None, extra={"obj_data": {"name_a": 1, "name_b": "fake_string", "name_c": 3.1415}})
        ])


    @patch('logging.getLogger')
    def test_write_info_non_serializable(self, mock_getLogger):
        '''Test writing log entry with a non serializable field'''
        non_serializable_byte_field = bytes()

        log.write(obj_data={"name_a": 1, "name_b": non_serializable_byte_field, "name_c": 3.1415})

        mock_getLogger.assert_has_calls([
            call().info(msg=None, extra={"obj_data": {"name_a": 1, "name_b": non_serializable_byte_field, "name_c": 3.1415}})
        ])


    @staticmethod
    @patch('idseq_dag.util.log.write')
    def test_log_event(mock_log_write):
        '''Test a successful log event'''
        log.log_event('fake_event_name', values={"name_a": 1, "name_b": "fake_string", "name_c": 3.1415}, warning=True, flush=False, extra_fields={"extra_field_1": "abc"})

        mock_log_write.assert_has_calls([
            call(obj_data={"event": "fake_event_name", "extra_field_1": "abc", "values": {"name_a": 1, "name_b": "fake_string", "name_c": 3.1415}}, warning=True, debug=False, flush=False)
        ])


    @staticmethod
    @patch('idseq_dag.util.log.write')
    def test_log_event_3(mock_log_write):
        '''Test a successful log event with duration'''
        start_time = datetime.datetime(2019, 5, 27, 16, 0, 0)
        now = start_time + datetime.timedelta(minutes=1, seconds=18)
        mock_datetime = Mock(datetime.datetime)
        mock_datetime.now.return_value = now

        with patch.object(datetime, 'datetime', mock_datetime):
            log.log_event('fake_event_name', values={"name_a": 1}, warning=False, flush=False, start_time=start_time)

        mock_log_write.assert_has_calls([
            call(obj_data={"event": "fake_event_name", "values": {"name_a": 1}, "duration_ms": 78000}, warning=False, debug=False, flush=False)
        ])


    def test_json_formatter(self):
        '''Test if logger is formatting to json output'''
        output = io.StringIO()
        handler = logging.StreamHandler(output)
        formatter = log.JsonFormatter()
        handler.setFormatter(formatter)
        logger = logging.Logger("my_test_logger")
        logger.addHandler(handler)

        with patch('logging.getLogger', return_value=logger):
            log.write(message='test message', obj_data={"abc": 1, "f2": bytes()})

        out_str = output.getvalue().rstrip()

        LOG_OUTPUT_REGEX = r'^{"timestamp": "\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}.\d{3}", "msg": "test message", "data": {"abc": 1, "f2": "<<non-serializable: bytes>>"}, "thread": "[^"]+", "pid": \d, "level": "INFO"}'
        # sample expected message
        self.assertRegex('{"timestamp": "2019-05-28T21:06:46.728", "msg": "test message", "data": {"abc": 1, "f2": "<<non-serializable: bytes>>"}, "thread": "MainThread", "pid": 6, "level": "INFO"}', LOG_OUTPUT_REGEX)
        # actual message
        self.assertRegex(out_str, LOG_OUTPUT_REGEX)
