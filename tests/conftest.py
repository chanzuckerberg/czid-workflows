import os
import json
import subprocess
import pytest


@pytest.fixture(scope="session")
def repo_dir():
    return os.path.dirname(os.path.dirname(__file__))


@pytest.fixture(scope="session")
def miniwdl_run(tmpdir_factory):
    """
    miniwdl_run(*args, returncode) passes the command-line arguments list (usually the WDL filename, followed by
    inputs) through to `miniwdl run`, checks its exit status == returncode, and returns its parsed JSON stdout.
    """
    tmpdir = tmpdir_factory.mktemp("miniwdl_run_")

    def kappa(*args, returncode=0):
        args = ["miniwdl", "run"] + list(args) + ["--verbose", "--error-json", "--dir", tmpdir]
        rslt = subprocess.run(args, stdout=subprocess.PIPE, check=False)
        assert rslt.returncode == returncode
        return json.loads(rslt.stdout)

    return kappa


@pytest.fixture(scope="session")
def miniwdl_inputs_outputs():
    """
    miniwdl_inputs_outputs(call_dir) parses the input & output JSON from a completed run directory. Use this to load
    one task's call directory from a previous end-to-end workflow test case, and manipulate the inputs to form the
    basis of new unit tests to exercise different code paths in that task.
    """

    def kappa(call_dir):
        with open(os.path.join(call_dir, "inputs.json")) as infile:
            inputs = json.load(infile)
        with open(os.path.join(call_dir, "outputs.json")) as infile:
            outputs = json.load(infile)
        return (inputs, outputs)

    return kappa


@pytest.fixture(scope="session")
def RunFailed_stderr_msg():
    """
    given RunFailed info (parsed error JSON), attempt to retrieve the error message written as JSON on the last line of
    the failed task's standard error stream
    """

    def kappa(info):
        assert info["error"] == "RunFailed"
        ans = None
        if info["cause"]["error"] == "CommandFailed":
            stderr_file = info["cause"]["stderr_file"]
            tail = subprocess.run(
                ["tail", "-n", "1", stderr_file], stdout=subprocess.PIPE, check=True
            )
            try:
                ans = json.loads(tail.stdout)
            except Exception:
                pass
        return ans

    return kappa
