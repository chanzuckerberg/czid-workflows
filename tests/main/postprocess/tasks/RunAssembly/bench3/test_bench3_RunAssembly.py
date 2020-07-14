
import os
import math
import pytest
import WDL


@pytest.fixture()
def inputs_outputs(exe, load_inputs_outputs):
    "(test_inputs, expected_outputs) for this case, as JSON-like dicts"
    return load_inputs_outputs(exe, os.path.dirname(__file__))


def test_bench3_RunAssembly(exe, inputs_outputs, miniwdl_run, compare_outputs):
    # Load the test inputs & expected outputs
    (inputs, expected_outputs) = inputs_outputs
    # Run & get the actual outputs
    run_dir, actual_outputs = miniwdl_run(exe, inputs)
    # Check actual outputs against the expected outputs
    # There seems to be slight variability in the output size.
    def file_digest(fn):
        return f"{os.path.basename(fn)}/{math.floor(os.path.getsize(fn)/1024) if os.path.isfile(fn) else None}"

    compare_outputs(exe, actual_outputs, expected_outputs, file_digest=file_digest)

# Optional: add more test cases, perhaps with perturbed inputs