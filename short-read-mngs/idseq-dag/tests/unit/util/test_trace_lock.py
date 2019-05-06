import unittest
from unittest.mock import patch, call

# Collaborators to inject
import threading
import time

# Class under test
from idseq_dag.util.trace_lock import TraceLock

class TestTraceLock(unittest.TestCase):
    '''Tests for `util/trace_rlock.py`'''

    @staticmethod
    @patch('idseq_dag.util.log.log_event')
    def test_aquire_free_lock(mock_log_event):
        '''Test lock acquiring when lock is free'''
        lock = TraceLock('lock1')

        with lock:
            pass

        mock_log_event.assert_has_calls([
            call('trace_lock', values={'state': 'acquired', 'lock_name': 'lock1', 'thread_name': 'MainThread'}),
            call('trace_lock', values={'state': 'released', 'lock_name': 'lock1', 'thread_name': 'MainThread'})
        ])

    @staticmethod
    @patch('idseq_dag.util.log.log_event')
    def test_aquire_lock(mock_log_event):
        '''Test waiting to acquire lock'''
        lock = TraceLock('lock2') # object under test
        event = threading.Event()
        def thread_run():
            with lock: # lock acquired
                event.set()
                time.sleep(0.5) # make main thread wait for this lock a bit
        first_thread = threading.Thread(target=thread_run)
        first_thread.name = "Thread-1"
        first_thread.start()

        event.wait()  # wait first_thread acquire the lock

        with lock: # lock will be held for a while we want
            pass

        mock_log_event.assert_has_calls([
            call('trace_lock', values={'lock_name': 'lock2', 'thread_name': 'Thread-1', 'state': 'acquired'}),
            call('trace_lock', values={'lock_name': 'lock2', 'thread_name': 'MainThread', 'state': 'waiting'}),
            call('trace_lock', values={'lock_name': 'lock2', 'thread_name': 'Thread-1', 'state': 'released'}),
            call('trace_lock', values={'lock_name': 'lock2', 'thread_name': 'MainThread', 'state': 'acquired_after_wait'}),
            call('trace_lock', values={'lock_name': 'lock2', 'thread_name': 'MainThread', 'state': 'released'})
        ])
