#!/bin/env python3
import threading
import traceback
import random
import time
from typing import Iterable, List, Any

class ThreadWithResult(threading.Thread):
    """A simple replacement for Python's built-in threading.Thread class,
    adding support for target funcs that return results or raise exceptions.
    Usage example:

        t = ThreadWithResult(target=fetch_something)
        t.start()
        t.join()
        assert t.completed and not t.exception
        print(t.result)
    """
    def __init__(self, target, args=(), kwargs=None):
        super(ThreadWithResult, self).__init__()
        self.args = args
        if kwargs is None:
            # This is a bit arcane, but we want to create a new empty dict
            # here, rather than having kwargs={} above, because in python
            # {} default values are singleton globals.  This matters
            # only in the most rare and exceptional circumstances.
            self.kwargs = {}
        else:
            self.kwargs = kwargs
        self.target = target
        self.exception = None
        self.completed = False
        self.result = None
        self.print_traceback = True

    def run(self) -> None:
        try:
            self.result = self.target(*self.args, **self.kwargs)
            self.exception = False
        except: #pylint: disable=bare-except
            if self.print_traceback:
                traceback.print_exc()
            self.exception = True
        finally:
            self.completed = True

def run_all(threads: Iterable[ThreadWithResult]) -> List[Any]:
    """Run the provided threads.  If all complete without raising an exception,
    return a list of their results."""
    threads = list(threads)
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    for i, t in enumerate(threads):
        assert t.completed and not t.exception, f"Problem in thread {i}."
    return [t.result for t in threads]

def mt_starmap(func, seq: Iterable) -> List[Any]:
    """Like the built-in starmap function, but runs each func invocation in its own thread."""
    return run_all(ThreadWithResult(func, args) for args in seq)

def mt_map(func, seq: Iterable):
    """Like the built-in map function, but runs each func invocation in its own thread."""
    return mt_starmap(func, ((arg,) for arg in seq))

# ONLY TESTS BELOW THIS LINE
# Run this module as a command to trigger the unit test.

def _unit_test_1(r=random.Random(time.time())):
    test = "thread_with_result: exception handling test"
    try:
        threads_to_fail = set([2, 3, 5, 7])
        def target(i: int) -> int:
            time.sleep(r.random())
            if i in threads_to_fail:
                raise RuntimeError()
            return i
        threads = [
            ThreadWithResult(target, (i,))
            for i in range(16)
        ]
        for i, t in enumerate(threads):
            if i in threads_to_fail:
                t.print_traceback = False
        try:
            run_all(threads)
        except AssertionError:
            pass
        except:  #pylint: disable=bare-except
            print("Unexpected exception.")
            raise
        else:
            assert False, "execute_all should have raised an AssertionError"
        for i, t in enumerate(threads):
            assert t.completed
            if i in threads_to_fail:
                assert t.exception, f"Thread {i} should have raised an exception."
                assert t.result is None, f"Thread {i} should not have returned a result."
            else:
                assert not t.exception, f"Thread {i} should not have raised an exception."
                assert t.result == i, f"Thread {i} should have returned result {i}."
        print(f"{test} passed")
    except:
        print(f"{test} failed")
        raise

def _unit_test_2(r=random.Random(time.time())):
    test = "thread_with_result: mt_map and mt_starmap test"
    try:
        numbers_to_modify = set([2, 3, 5, 7])
        def target(i: int) -> int:
            time.sleep(r.random())
            if i in numbers_to_modify:
                return 100 + i
            return i
        try:
            results = mt_map(target, range(16))
        except:  #pylint: disable=bare-except
            print("Unexpected exception.")
            raise
        for i, r in enumerate(results):
            if i in numbers_to_modify:
                assert r == 100 + i
            else:
                assert r == i
        print(f"{test} passed")
    except:
        print(f"{test} failed")
        raise

if __name__ == "__main__":
    _unit_test_1()
    _unit_test_2()
