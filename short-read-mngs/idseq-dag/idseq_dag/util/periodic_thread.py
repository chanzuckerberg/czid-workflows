import threading
import time
import traceback

class PeriodicThread(threading.Thread):
    '''
    Thread that keeps executing its target func periodically until it is told to stop.
    Usage example:
        t = PeriodicThread(target=poll_something, wait_seconds=60, stop_latency_seconds=10, args=(some_argument_1, some_argument_2))
        t.start()
        ...
        t.stop()
        t.join()
    '''

    def __init__(self, target, wait_seconds, stop_latency_seconds, args=(), kwargs=None):
        super(PeriodicThread, self).__init__()
        self.target = target
        self.wait_seconds = wait_seconds
        self.stop_latency_seconds = stop_latency_seconds
        self.args = args
        if kwargs is None:
            self.kwargs = {}
        else:
            self.kwargs = kwargs
        self._stop_event = threading.Event()

    def stop(self):
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()

    def run(self) -> None:
        last_target_run_unixtime = int(time.time()) - self.stop_latency_seconds
        while not self.stopped():
            unixtime_now = int(time.time())
            if unixtime_now - last_target_run_unixtime >= self.wait_seconds:
                try:
                    self.target(*self.args, **self.kwargs)
                except:
                    traceback.print_exc()
                finally:
                    last_target_run_time = unixtime_now  # noqa
            time.sleep(self.stop_latency_seconds)
