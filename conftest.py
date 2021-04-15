"""
This pytest fixture "util" provides an easy way for test cases in WORKFLOWNAME/test directories to
import the test_util/ helper functions, without manipulating PYTHONPATH or such.
"""

import sys
import os
import pytest

sys.path.append(os.path.dirname(__file__))
import test_util  # noqa


@pytest.fixture(scope="session")
def util(tmpdir_factory):
    if not test_util._miniwdl_run_dir:
        test_util._miniwdl_run_dir = str(tmpdir_factory.getbasetemp() / "miniwdl_run")
    return test_util
