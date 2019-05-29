import multiprocessing
import threading
import idseq_dag.util.log as log

class TraceLock():
    r"""
    This class is a wrapper to RLocks that can log states of the lock for each thread.

    TraceLock can be used in a `with` block, like a regular RLock:

    ```
    with TraceLock("my_lock_name"):
        ... something that requires a lock
    ```

    This will generate log entries like:
    ```
        {"event":"trace_lock", values={"lock_name": "lock1", "thread_name": "Thread-1", "state": "acquired"})
        {"event":"trace_lock", values={"lock_name": "lock1", "thread_name": "Thread-2", "state": "waiting"})
        {"event":"trace_lock", values={"lock_name": "lock1", "thread_name": "Thread-1", "state": "released"})
        {"event":"trace_lock", values={"lock_name": "lock1", "thread_name": "Thread-2", "state": "acquired_after_wait"})
        {"event":"trace_lock", values={"lock_name": "lock1", "thread_name": "Thread-2", "state": "released"})
    ```

    State diagram:
    ```
           some thread is trying               
             to acquire a lock                 
                     ||                        
                     /\                        
           +-------- \/ --------+ lock is      
 lock is   |                    | not available
 available |              +-----v----+         
           |              | waiting  |         
           |              +-----|----+         
           |                    | lock is now  
           |                    | available    
     +-----v----+     +---------v---------+    
     | acquired |     |acquired_after_wait|    
     +-----|----+     +---------|---------+    
           |                    |              
           |   thread realeases |              
           |      the lock      |              
           +-------->/\<--------+              
                     \/                        
                     ||                        
                     vv                        
              +--------------+                 
              |   released   |                 
              +--------------+
    ```
    """
    def __init__(self, lock_name, lock=multiprocessing.RLock(), debug=True):
        self._lock = lock
        self._lock_name = lock_name
        self.debug = debug

    def acquire(self):
        v = {"lock_name": self._lock_name, "thread_name": threading.current_thread().name}
        if self._lock.acquire(False):
            log.log_event("trace_lock", values={**v, "state": "acquired"}, debug=self.debug)
        else:
            log.log_event("trace_lock", values={**v, "state": "waiting"}, debug=self.debug)
            self._lock.acquire(True)
            log.log_event("trace_lock", values={**v, "state": "acquired_after_wait"}, debug=self.debug)

    def __enter__(self):
        self.acquire()

    def release(self):
        log.log_event("trace_lock", values={"lock_name": self._lock_name,
                                            "thread_name": threading.current_thread().name,
                                            "state": "released"}, debug=self.debug)
        self._lock.release()

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.release()