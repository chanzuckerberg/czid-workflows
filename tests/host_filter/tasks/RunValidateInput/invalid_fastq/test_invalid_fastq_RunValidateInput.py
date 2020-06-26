import os
import pytest
import WDL


@pytest.fixture(scope="session")
def inputs_outputs(exe, load_inputs_outputs):
    "(test_inputs, expected_outputs) for this case, as JSON-like dicts"
    return load_inputs_outputs(exe, os.path.dirname(__file__))


def test_invalid_fastq_RunValidateInput(exe, inputs_outputs, miniwdl_run, RunFailed_stderr_json):
    # for this case ./expected_outputs.json is {} because we expect failure
    (inputs, expected_outputs) = inputs_outputs
    # run on inputs and catch the resulting RunFailed exception
    with pytest.raises(WDL.runtime.RunFailed) as exc:
        miniwdl_run(exe, inputs)
    # RunFailed_stderr_json helper extracts the last stderr line of the failed task and parses it
    # as JSON if possible, otherwise returns None.
    info = RunFailed_stderr_json(exc.value)
    # check for expected error message
    assert "File does not follow fastq format" in str(info)
