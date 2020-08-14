import os
import pytest
import WDL


@pytest.fixture()
def inputs_outputs(exe, load_inputs_outputs):
    "(test_inputs, expected_outputs) for this case, as JSON-like dicts"
    return load_inputs_outputs(exe, os.path.dirname(__file__))


def test_bench3_ercc_db_RunStar(exe, inputs_outputs, miniwdl_run, compare_outputs):
    """
    Runs the bench3 inputs using the 1.6GB ERCC index instead of the 27GB human index
    """
    # Load the test inputs & expected outputs
    (inputs, expected_outputs) = inputs_outputs
    # Run & get the actual outputs
    run_dir, actual_outputs = miniwdl_run(exe, inputs)
    # Check actual outputs against the expected outputs
    # Optional:
    # - del keys from expected_outputs that needn't be checked
    # - custom assertions on actual/expected
    compare_outputs(exe, actual_outputs, expected_outputs)


# Optional: add more test cases, perhaps with perturbed inputs
