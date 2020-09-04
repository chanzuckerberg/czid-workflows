import sys
import os
import json
import logging
import subprocess
import pytest
import WDL

# In any test case inputs.json, the following keys will be overridden with the corresponding value
# IF PRESENT. So the inputs.json can just put null expecting this to fill it.
INPUT_OVERRIDES = {
    "docker_image_id": os.environ.get(
        "DOCKER_IMAGE_ID",
        "docker.pkg.github.com/chanzuckerberg/idseq-workflows/idseq-main-public:3db222b",
    ),
    "s3_wd_uri": "s3://DUMMY_URI/",
}

@pytest.fixture
def repo_dir():
    return os.path.dirname(os.path.dirname(__file__))


@pytest.fixture(scope="session")
def miniwdl_run_cfg():
    "miniwdl run session-wide initialization"
    logging.basicConfig(level=logging.DEBUG)
    with WDL._util.configure_logger():
        logger = logging.getLogger(__name__)
        cfg = WDL.runtime.config.Loader(logger)
        return cfg


@pytest.fixture(scope="session")
def load_task():
    "function to load the WDL task object"
    return lambda wdl_path, task_name: next(
        task for task in WDL.load(wdl_path).tasks if task.name == task_name
    )


@pytest.fixture(scope="session")
def load_inputs_outputs():
    "function to load inputs.json and expected_outputs.json for a test case"

    def kappa(exe, call_dir):
        with open(os.path.join(call_dir, "inputs.json")) as infile:
            inputs = json.load(infile)
        for key in INPUT_OVERRIDES:
            if key in inputs:
                inputs[key] = INPUT_OVERRIDES[key]
        # absolutify paths to make cwd irrelevant
        inputs = WDL.values_from_json(inputs, exe.available_inputs, exe.required_inputs)
        inputs = WDL.Value.rewrite_env_files(
            inputs, lambda fn: os.path.join(call_dir, fn) if "://" not in fn else fn
        )
        inputs = WDL.values_to_json(inputs)
        with open(os.path.join(call_dir, "expected_outputs.json")) as infile:
            expected_outputs = WDL.values_from_json(
                json.load(infile), exe.effective_outputs, namespace=exe.name
            )
        expected_outputs = WDL.Value.rewrite_env_files(
            expected_outputs, lambda fn: os.path.join(call_dir, fn)
        )
        expected_outputs = WDL.values_to_json(expected_outputs, exe.name)
        return (inputs, expected_outputs)

    return kappa


@pytest.fixture
def miniwdl_run(miniwdl_run_cfg, tmpdir):
    "function to run the exe on given inputs (not necessarily the saved ones, e.g. to test error cases)"

    def kappa(exe, inputs, **kwargs):
        inputs = WDL.values_from_json(inputs, exe.available_inputs, exe.required_inputs)
        run_dir, outputs = WDL.runtime.run(
            miniwdl_run_cfg, exe, inputs, run_dir=str(tmpdir), **kwargs
        )
        return (run_dir, WDL.values_to_json(outputs, exe.name))

    return kappa


@pytest.fixture(scope="session")
def compare_outputs():
    "function to perform structured comparison of WDL outputs"

    def kappa(exe, actual, expected, file_digest=None, float_places=None):
        if file_digest is None:
            # default: represent each file just by its basename and size
            file_digest = (
                lambda fn: f"{os.path.basename(fn)}/{os.path.getsize(fn) if os.path.isfile(fn) else None}"
            )
        actual = WDL.values_from_json(actual, exe.effective_outputs, namespace=exe.name)
        actual = WDL.Value.rewrite_env_files(actual, file_digest)
        actual = _round_floats(WDL.values_to_json(actual, exe.name), float_places)
        expected = WDL.values_from_json(expected, exe.effective_outputs, namespace=exe.name)
        expected = WDL.Value.rewrite_env_files(expected, file_digest)
        expected = _round_floats(WDL.values_to_json(expected, exe.name), float_places)
        try:
            _compare_outputs(actual, expected)
        except AssertionError:
            print(json.dumps({"actual": actual, "expected": expected}, indent=2), file=sys.stderr)
            raise

    return kappa


def _round_floats(j, places):
    if places is not None:
        if isinstance(j, list):
            return [_round_floats(elt, places) for elt in j]
        if isinstance(j, dict):
            ans = {}
            for key in j:
                ans[key] = _round_floats(j[key], places)
            return ans
        if isinstance(j, float):
            return round(j, places)
    return j


def _compare_outputs(actual, expected):
    if isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected)
        for i in range(len(expected)):
            _compare_outputs(actual[i], expected[i])
    elif isinstance(expected, dict):
        assert isinstance(actual, dict) and not (set(expected) - set(actual))
        for key in expected:
            _compare_outputs(actual[key], expected[key])
    elif isinstance(expected, (str, int, float)):
        assert type(actual) is type(expected) and actual == expected
    elif expected is None:
        assert actual is None
    else:
        assert False


@pytest.fixture(scope="session")
def RunFailed_stderr_json():
    "given WDL.runtime.RunFailed exception, attempt to JSON-parse the last line of failing task's standard error"

    def kappa(exc):
        ans = None
        if isinstance(exc, WDL.runtime.RunFailed):
            info = WDL.runtime.error_json(exc)
            if info["cause"]["error"] == "CommandFailed":
                # read last line of failed task's standard error
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
