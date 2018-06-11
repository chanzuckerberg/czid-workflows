import multiprocessing
import random
import subprocess
import sys
import threading
import time
from functools import wraps
import idseq_dag.util.log as log


class Updater(object):
    """Base for CommandTracker."""

    def __init__(self, update_period, update_function):
        self.update_period = update_period
        self.update_function = update_function
        self.timer_thread = None
        self.exited = False
        self.t_start = time.time()

    def relaunch(self, initial_launch=False):
        if self.exited:
            return
        if self.timer_thread and not initial_launch:
            t_elapsed = time.time() - self.t_start
            self.update_function(t_elapsed)
        self.timer_thread = threading.Timer(self.update_period, self.relaunch)
        self.timer_thread.start()

    def __enter__(self):
        self.relaunch(initial_launch=True)
        return self

    def __exit__(self, *args):
        self.timer_thread.cancel()
        self.exited = True


class CommandTracker(Updater):
    """CommandTracker is for running external and remote commands and
    monitoring their progress with log updates and timeouts.
    """
    lock = multiprocessing.RLock()
    count = multiprocessing.Value('i', 0)

    def __init__(self, update_period=15):
        super(CommandTracker, self).__init__(
            update_period, self.print_update_and_enforce_timeout)
        # User can set the watchdog to a function that takes self.id and
        # t_elapsed as single arg
        self.proc = None  # Value indicates registered subprocess.
        self.timeout = None
        self.t_sigterm_sent = None  # First sigterm, then sigkill.
        self.t_sigkill_sent = None
        self.grace_period = update_period / 2.0
        with CommandTracker.lock:
            self.id = CommandTracker.count.value
            CommandTracker.count.value += 1

    def print_update_and_enforce_timeout(self, t_elapsed):
        """Log an update after every polling period to indicate the command is
        still active.
        """
        with log.print_lock:
            if self.proc is None or self.proc.poll() is None:
                log.write("Command %d still running after %3.1f seconds." %
                      (self.id, t_elapsed))
            else:
                # This should be uncommon, unless there is lengthy python
                # processing following the command in the same CommandTracker
                # "with" block. Note: Not to be confused with post-processing
                # on the data.
                log.write("Command %d still postprocessing after %3.1f seconds." %
                      (self.id, t_elapsed))
            sys.stdout.flush()
        self.enforce_timeout(t_elapsed)

    def enforce_timeout(self, t_elapsed):
        """Check the timeout and send SIGTERM then SIGKILL to end a command's
        execution.
        """
        if self.timeout is None or not self.proc or \
                t_elapsed <= self.timeout or self.proc.poll() is not None:
            # Skip if unregistered subprocess, subprocess not yet timed out,
            # or subprocess already exited.
            pass
        elif not self.t_sigterm_sent:
            # Send SIGTERM first.
            with log.print_lock:
                msg = "Command %d has exceeded timeout of %3.1f seconds. " \
                      "Sending SIGTERM." % (self.id, self.timeout)
                log.write(msg)
                sys.stdout.flush()
            self.t_sigterm_sent = time.time()
            self.proc.terminate()
        elif not self.t_sigkill_sent:
            # Grace_period after SIGTERM, send SIGKILL.
            if time.time() > self.t_sigterm_sent + self.grace_period:
                with log.print_lock:
                    msg = "Command %d still alive %3.1f seconds after " \
                          "SIGTERM. Sending SIGKILL." % (self.id, time.time() - self.t_sigterm_sent)
                    log.write(msg)
                    sys.stdout.flush()
                self.t_sigkill_sent = time.time()
                self.proc.kill()
        else:
            with log.print_lock:
                msg = "Command %d still alive %3.1f seconds after " \
                      "SIGKILL." % (self.id, time.time() - self.t_sigkill_sent)
                log.write(msg)
                sys.stdout.flush()


class ProgressFile(object):
    def __init__(self, progress_file):
        self.progress_file = progress_file
        self.tail_subproc = None

    def __enter__(self):
        # TODO: Do something else here. Tail gets confused if the file
        # pre-exists. Also need to rate-limit.
        if self.progress_file:
            self.tail_subproc = subprocess.Popen(
                "touch {pf} ; tail -f {pf}".format(pf=self.progress_file),
                shell=True)
        return self

    def __exit__(self, *args):
        if self.tail_subproc:
            # TODO: Do we need to join the tail subproc after killing it?
            self.tail_subproc.kill()


def run_in_subprocess(target):
    """
    Decorator that executes a function synchronously in a subprocess.
    Use case:

        thread 1:
             compute_something(x1, y1, z1)

        thread 2:
             compute_something(x2, y2, z2)

        thread 3:
             compute_something(x3, y3, z3)

    If compute_something() is CPU-intensive, the above threads won't really run
    in parallel because of the Python global interpreter lock (GIL). To avoid
    this problem without changing the above code, simply decorate the definition
    of compute_something() like so:

         @run_in_subprocess
         def compute_something(x, y, z):
             ...

    Typical subprocess limitations or caveats apply:
       a. The caller can't see any value returned by the decorated function.
          It should output to a file, a pipe, or a multiprocessing queue.
       b. Changes made to global variables won't be seen by parent process.
       c. Use multiprocessing semaphores/locks/etc, not their threading versions.

    Tip: If input from the same file is needed in all invocations of the
    decorated function, do the I/O before the first call, to avoid accessing
    the file multiple times.
    """

    @wraps(target)
    def wrapper(*args, **kwargs):
        p = multiprocessing.Process(target=target, args=args, kwargs=kwargs)
        p.start()
        p.join()
        if p.exitcode != 0:
            raise RuntimeError("Failed {} on {}, {}".format(
                target.__name__, args, kwargs))

    return wrapper


def retry(operation, randgen=random.Random().random):
    """Retry decorator for external commands."""
    # Note the use of a separate random generator for retries so transient
    # errors won't perturb other random streams used in the application.
    @wraps(operation)
    def wrapped_operation(*args, **kwargs):
        remaining_attempts = 3
        delay = 1.0
        while remaining_attempts > 1:
            try:
                return operation(*args, **kwargs)
            except:
                # Random jitter and exponential delay
                time.sleep(delay * (1.0 + randgen()))
                delay *= 3.0
                remaining_attempts -= 1
        # The last attempt is outside try/catch so caller can handle exception
        return operation(*args, **kwargs)

    return wrapped_operation


def execute(command,
            progress_file=None,
            timeout=None,
            grace_period=None,
            capture_stdout=False,
            merge_stderr=False):
    """Primary way to start external commands in subprocesses and handle
    execution with logging.
    """
    with CommandTracker() as ct:
        with log.print_lock:
            log.write("Command {}: {}".format(ct.id, command))
        with ProgressFile(progress_file):
            if timeout:
                ct.timeout = timeout
            if grace_period:
                ct.grace_period = grace_period
            if capture_stdout:
                # Capture only stdout. Child stderr = parent stderr unless
                # merge_stderr specified. Child input = parent stdin.
                ct.proc = subprocess.Popen(
                    command,
                    shell=True,
                    stdin=sys.stdin.fileno(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT
                    if merge_stderr else sys.stderr.fileno())
                stdout, _ = ct.proc.communicate()
            else:
                # Capture nothing. Child inherits parent stdin/out/err.
                ct.proc = subprocess.Popen(command, shell=True)
                ct.proc.wait()
                stdout = None

            if ct.proc.returncode:
                raise subprocess.CalledProcessError(ct.proc.returncode,
                                                    command, stdout)
            if capture_stdout:
                return stdout


def execute_with_output(command,
                        progress_file=None,
                        timeout=None,
                        grace_period=None,
                         merge_stderr=False):
    return execute(
        command,
        progress_file,
        timeout,
        grace_period,
        capture_stdout=True,
        merge_stderr=merge_stderr)


def execute_remote(base_command, key_path, remote_username, instance_ip):
    # ServerAliveInterval to fix issue with containers keeping open an SSH
    # connection even after worker machines had finished running.
    return 'ssh -o "StrictHostKeyChecking no" -o "ConnectTimeout 15" ' \
           '-o "ServerAliveInterval 60" -i %s %s@%s "%s"' % (
               key_path, remote_username, instance_ip, base_command)


def scp(key_path, remote_username, instance_ip, remote_path, local_path):
    assert " " not in key_path
    assert " " not in remote_path
    assert " " not in local_path
    # ServerAliveInterval to fix issue with containers keeping open an SSH
    # connection even after worker machines had finished running.
    return 'scp -o "StrictHostKeyChecking no" -o "ConnectTimeout 15" ' \
           '-o "ServerAliveInterval 60" -i {key_path} ' \
           '{username}@{ip}:{remote_path} {local_path}'.format(
        key_path=key_path,
        username=remote_username,
        ip=instance_ip,
        remote_path=remote_path,
        local_path=local_path)
