import multiprocessing
import random
import subprocess
import sys
import os
import threading
import time
import glob as _glob
import shutil
import pathlib
from functools import wraps
from typing import Union
import pkg_resources
import idseq_dag.util.log as log
from idseq_dag.util.trace_lock import TraceLock
import idseq_dag.util.command_patterns as command_patterns


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
        self.timer_thread.name = "TimerThread"
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
    lock = TraceLock("CommandTracker", multiprocessing.RLock())
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
        if self.proc is None or self.proc.poll() is None:
            log.write("Command %d still running after %3.1f seconds." %
                      (self.id, t_elapsed))
        else:
            # This should be uncommon, unless there is lengthy python
            # processing following the command in the same CommandTracker
            # "with" block. Note: Not to be confused with post-processing
            # on the data.
            log.write(
                "Command %d still postprocessing after %3.1f seconds." %
                (self.id, t_elapsed))
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
            msg = "Command %d has exceeded timeout of %3.1f seconds. " \
                "Sending SIGTERM." % (self.id, self.timeout)
            log.write(msg)
            self.t_sigterm_sent = time.time()
            self.proc.terminate()
        elif not self.t_sigkill_sent:
            # Grace_period after SIGTERM, send SIGKILL.
            if time.time() > self.t_sigterm_sent + self.grace_period:
                msg = "Command %d still alive %3.1f seconds after " \
                    "SIGTERM. Sending SIGKILL." % (self.id, time.time() - self.t_sigterm_sent)
                log.write(msg)
                self.t_sigkill_sent = time.time()
                self.proc.kill()
        else:
            msg = "Command %d still alive %3.1f seconds after " \
                "SIGKILL." % (self.id, time.time() - self.t_sigkill_sent)
            log.write(msg)


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
                shell=True, executable="/bin/bash")
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
        with log.print_lock:
            p = multiprocessing.Process(target=target, args=args, kwargs=kwargs)
            p.start()
        p.join()
        """
        Occasionally after p.join() returns, the exit code is `None` for a short amount of time.
        As a mitigation, we wait for max 10 seconds for the exit code to return, then check the error code.
        """
        timeout = 0
        while p.exitcode is None and timeout < 10:
            time.sleep(1)
            timeout += 1
        if p.exitcode != 0:
            raise RuntimeError(f"Failed {target.__qualname__} with code {p.exitcode} on {list(args)}, {kwargs}")  # singleton list prints prettier than singleton tuple
    return wrapper


def retry(operation, MAX_TRIES=3):
    """Retry decorator for external commands."""
    # Note the use of a separate random generator for retries so transient
    # errors won't perturb other random streams used in the application.
    invocation = [0]  # a python trick, without which the nested function would not increment a counter defined in the parent

    @wraps(operation)
    def wrapped_operation(*args, **kwargs):
        randgen = None
        remaining_attempts = MAX_TRIES
        delay = 1.0
        while remaining_attempts > 1:
            try:
                t_start = time.time()
                return operation(*args, **kwargs)
            except:
                t_end = time.time()
                if randgen == None:
                    invocation[0] += 1
                    randgen = random.Random(os.getpid() * 10000 + invocation[0]).random  # seed based on pid so subprocesses won't retry in lockstep
                if t_end - t_start > 30:
                    # For longer operations, the delay should increase, so that the backoff will meaningfully reduce load on the failing service.
                    delay = (t_end - t_start) * 0.2
                wait_time = delay * (1.0 + 2.0 * randgen())
                log.write(f"Sleeping {wait_time} seconds before retry {MAX_TRIES - remaining_attempts + 1} of {operation} with {args}, {kwargs}.")
                time.sleep(wait_time)
                delay *= 3.0
                remaining_attempts -= 1
        # The last attempt is outside try/catch so caller can handle exception
        return operation(*args, **kwargs)
    return wrapped_operation


@retry
def execute_with_retries(command,
                         progress_file=None,
                         timeout=None,
                         grace_period=None,
                         merge_stderr=False,
                         log_context_mode=log.LogContextMode.START_END_LOG_EVENTS):
    execute(
        command=command,
        progress_file=progress_file,
        timeout=timeout,
        grace_period=grace_period,
        merge_stderr=merge_stderr,
        log_context_mode=log_context_mode
    )


def execute(command: Union[command_patterns.CommandPattern, str],
            progress_file: str = None,
            timeout: int = None,
            grace_period: int = None,
            capture_stdout: bool = False,
            merge_stderr: bool = False,
            log_context_mode: log.LogContextMode = log.LogContextMode.START_END_LOG_EVENTS) -> Union[str, None]:
    """Primary way to start external commands in subprocesses and handle
    execution with logging.
    """
    if not isinstance(command, command_patterns.CommandPattern):
        # log warning if using legacy format
        log.write(warning=True, message="Command parameter is using legacy type str. Use idseq_dag.util.command_patterns.", obj_data={"cmd": command, "type": type(command)})
        cmd = command_patterns.ShellScriptCommand(script=command, args=[])
    else:
        cmd = command

    with CommandTracker() as ct:
        log_values = {"cid": f"Command {ct.id}", "command": cmd.as_dict()}
        with log.log_context('command_execute', values=log_values, log_context_mode=log_context_mode) as lctx:
            with ProgressFile(progress_file):
                if timeout:
                    ct.timeout = timeout
                if grace_period:
                    ct.grace_period = grace_period
                if capture_stdout:
                    # Capture only stdout. Child stderr = parent stderr unless
                    # merge_stderr specified. Child input = parent stdin.
                    ct.proc = cmd.open(
                        stdin=sys.stdin.fileno(),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT if merge_stderr else sys.stderr.fileno()
                    )
                    stdout, _ = ct.proc.communicate()
                else:
                    # Capture nothing. Child inherits parent stdin/out/err.
                    ct.proc = cmd.open()
                    ct.proc.wait()
                    stdout = None

                lctx.values.update({"returncode": ct.proc.returncode})

                if ct.proc.returncode:
                    raise subprocess.CalledProcessError(ct.proc.returncode,
                                                        str(command), stdout)
                if capture_stdout:
                    return stdout


def execute_with_output(command: Union[command_patterns.CommandPattern, str],
                        progress_file: str = None,
                        timeout: int = None,
                        grace_period: int = None,
                        merge_stderr: bool = False,
                        log_context_mode: log.LogContextMode = log.LogContextMode.START_END_LOG_EVENTS):
    return execute(
        command=command,
        progress_file=progress_file,
        timeout=timeout,
        grace_period=grace_period,
        capture_stdout=True,
        merge_stderr=merge_stderr,
        log_context_mode=log_context_mode
    ).decode('utf-8')


def make_dirs(path: str):
    if not os.path.isdir(path):
        with log.log_context(context_name="command.make_dirs", values={'path': path}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
            os.makedirs(path, exist_ok=True)


def write_text_to_file(text: str, file_path: str):
    with log.log_context(context_name='command.write_text_to_file', values={'path': file_path, 'text': text}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        with open(file_path, "w") as f:
            print(text, file=f)


def copy_file(src: str, dst: str):
    with log.log_context(context_name='command.copy_file', values={'src': src, 'dest': dst}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        shutil.copy(src, dst)


def move_file(src: str, dst: str):
    with log.log_context(context_name='command.move_file', values={'src': src, 'dest': dst}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        shutil.move(src, dst)


def rename(src: str, dst: str):
    with log.log_context(context_name='command.rename', values={'src': src, 'dest': dst}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        os.rename(src, dst)


def touch(path, exist_ok=True):
    with log.log_context(context_name='command.touch', values={'path': path}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        pathlib.Path(path).touch(exist_ok=exist_ok)


def remove_file(file_path: str):
    with log.log_context(context_name='command.remove_file', values={'path': file_path}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        os.remove(file_path)


def remove_rf(path: str):
    '''Mimics behavior of rm -rf linux command.'''
    def _remove_entry(path_entry):
        if os.path.isdir(path_entry) and not os.path.islink(path_entry):
            shutil.rmtree(path_entry)
        elif os.path.exists(path_entry):
            os.remove(path_entry)
    with log.log_context(context_name='command.remove_rf', values={'path': path}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        path_list = _glob.glob(path)
        if len(path_list) == 1 and path_list[0] == path:
            _remove_entry(path)
        else:
            for path_entry in path_list:
                with log.log_context(context_name='command.remove_rf._remove_entry', values={'path_entry': path_entry}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
                    _remove_entry(path_entry)


def chmod(path: str, mode: int):
    '''Execute a chmod operation.
       Parameter 'mode' must be in octal format. Ex: chmod('/tmp/test.txt', 0o400)'''
    with log.log_context(context_name='command.chmod', values={'path': path, 'mode': oct(mode)}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        os.chmod(path, mode)


def glob(glob_pattern: str, strip_folder_names: bool = False, max_results: int = 0):
    '''
        Execute a glob pattern to local file system.
            Parameters:
                glob_pattern(str): A glob pattern. Ex: /tmp/*.gz
                max_results(int): Limit the number of results to be returned. Zero means not limit is set.
                strip_folder_names(bool): Return only the file names without folder information.
                                          Ex: "/tmp/123.txt" is returned as "123.txt"
            Returns:
                Array of strings containing the files found. Empty array if none is found.
    '''
    values = {'glob_pattern': glob_pattern, 'strip_folder_names': strip_folder_names, 'max_results': max_results}
    with log.log_context(context_name='command.glob', values=values, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT):
        results = _glob.glob(glob_pattern)
        results.sort()
        if max_results > 0:
            results = results[:max_results]
        if strip_folder_names:
            results = list(map(os.path.basename, results))
        values["results"] = results
        return results


def scp(key_path, remote_username, instance_ip, remote_path, local_path):
    assert " " not in key_path
    assert " " not in remote_path
    assert " " not in local_path
    # ServerAliveInterval to fix issue with containers keeping open an SSH
    # connection even after worker machines had finished running.
    return command_patterns.SingleCommand(
        cmd="scp",
        args=[
            "-o", "StrictHostKeyChecking no",
            "-o", "ConnectTimeout 15",
            "-o", "ServerAliveInterval 60",
            "-i", key_path,
            f"{remote_username}@{instance_ip}:{remote_path}",
            local_path
        ]
    )


def remote(base_command, key_path, remote_username, instance_ip):
    # ServerAliveInterval to fix issue with containers keeping open an SSH
    # connection even after worker machines had finished running.
    return command_patterns.SingleCommand(
        cmd="ssh",
        args=[
            "-o", "StrictHostKeyChecking no",
            "-o", "ConnectTimeout 15",
            "-o", "ServerAliveInterval 60",
            "-i", key_path,
            f"{remote_username}@{instance_ip}",
            base_command
        ]
    )


def get_resource_filename(root_relative_path, package='idseq_dag'):
    '''
        Given a file path relative to the root of the package, it returns an absolute path.

        Example:
            command.get_resource_filename("scripts/fastq-fasta-line-validation.awk")

        will return a string containing:
            /app/idseq_dag/scripts/fastq-fasta-line-validation.awk
    '''
    return pkg_resources.resource_filename(package, root_relative_path)


class LongRunningCodeSection(Updater):
    """
    Make sure we print something periodically while a long section of code is running.
    """
    lock = TraceLock("LongRunningCodeSection", multiprocessing.RLock())
    count = multiprocessing.Value('i', 0)

    def __init__(self, name, update_period=15):
        super(LongRunningCodeSection, self).__init__(
            update_period, self.print_update)
        with LongRunningCodeSection.lock:
            self.id = LongRunningCodeSection.count.value
            LongRunningCodeSection.count.value += 1
        self.name = name

    def print_update(self, t_elapsed):
        """Log an update after every polling period to indicate the code section is
        still active.
        """
        log.write("LongRunningCodeSection %d (%s) still running after %3.1f seconds." %
                  (self.id, self.name, t_elapsed))
