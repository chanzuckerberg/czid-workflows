"""
Helper routines for idseq-workflows unit tests: boilerplate for invoking miniwdl locally,
collecting its results, etc.

Most test cases access these via the 'util' fixture (defined in ../conftest.py)
"""

import json
import os
import subprocess
import tempfile
import unittest
from subprocess import check_output
from pathlib import Path


def repo_dir():
    return Path(os.path.dirname(os.path.dirname(__file__)))


_miniwdl_run_dir = None


def miniwdl_run(*args, returncode=0):
    """
    miniwdl_run(*args, returncode) passes the command-line arguments list (usually the WDL filename, followed by
    inputs) through to `miniwdl run`, checks its exit status == returncode, and returns its parsed JSON stdout.
    """
    global _miniwdl_run_dir
    assert _miniwdl_run_dir

    args = (
        ["miniwdl", "run"]
        + list(str(arg) for arg in args)
        + ["--verbose", "--error-json", "--dir", _miniwdl_run_dir]
    )
    rslt = subprocess.run(args, stdout=subprocess.PIPE, check=False)
    assert (
        rslt.returncode == returncode
    ), f"miniwdl run returncode {rslt.returncode} (expected {returncode})"
    return json.loads(rslt.stdout)


def miniwdl_inputs_outputs(call_dir):
    """
    miniwdl_inputs_outputs(call_dir) parses the input & output JSON from a completed run directory. Use this to load
    one task's call directory from a previous end-to-end workflow test case, and manipulate the inputs to form the
    basis of new unit tests to exercise different code paths in that task.
    """

    with open(os.path.join(call_dir, "inputs.json")) as infile:
        inputs = json.load(infile)
    with open(os.path.join(call_dir, "outputs.json")) as infile:
        outputs = json.load(infile)
    return (inputs, outputs)


def wdl_error_message(runfailed):
    """
    given RunFailed info (parsed error JSON), attempt to retrieve the error message written as JSON on the last line of
    the failed task's standard error stream
    """
    ans = None
    if runfailed.get("error") == "RunFailed" and runfailed["cause"]["error"] == "CommandFailed":
        stderr_file = runfailed["cause"]["stderr_file"]
        tail = subprocess.run(["tail", "-n", "1", stderr_file], stdout=subprocess.PIPE, check=True)
        try:
            ans = json.loads(tail.stdout)
        except Exception:
            pass
    return ans


class WDLTestCase(unittest.TestCase):
    def run_miniwdl(self, args=[], task=None, task_input=None, docker_image_id=os.environ["DOCKER_IMAGE_ID"]):
        cmd = ["miniwdl", "run", "--verbose", self.wdl] + args + [f"docker_image_id={docker_image_id}"]
        if task:
            cmd += ["--task", task]
            if task_input:
                tf = tempfile.NamedTemporaryFile(prefix="idseq-workflows-test-", suffix=".json", mode="wt")
                json.dump(task_input, tf)
                tf.flush()
                cmd += ["--input", tf.name]
        else:
            tf = tempfile.NamedTemporaryFile(prefix="idseq-workflows-test-", suffix=".json", mode="wt")
            json.dump({k: v for k, v in self.common_inputs.items() if k+'=' not in ''.join(cmd)}, tf)
            tf.flush()
            cmd += ["--input", tf.name]
        td = tempfile.TemporaryDirectory(prefix="idseq-workflows-test-").name
        cmd += ["--verbose", "--error-json", "--dir", td]
        print(cmd)
        res = check_output(cmd)
        return json.loads(res)

    def assertRunFailed(self, ecm, task, error, cause):
        miniwdl_error = json.loads(ecm.exception.output)
        self.assertEqual(miniwdl_error["error"], "RunFailed")
        self.assertEqual(miniwdl_error["cause"]["error"], "CommandFailed")
        self.assertIn(miniwdl_error["cause"]["run"], {task, f"call-{task}"})
        with open(miniwdl_error["cause"]["stderr_file"]) as fh:
            last_line = fh.read().splitlines()[-1]
            idseq_error = json.loads(last_line)
        self.assertEqual(idseq_error["wdl_error_message"], True)
        self.assertEqual(idseq_error["error"], error)
        self.assertEqual(idseq_error["cause"], cause)
