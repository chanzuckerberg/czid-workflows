
import os
import pytest
import WDL

@pytest.fixture(scope="session")
def exe(load_task):
    "load the task to be tested"
    return load_task(
        os.path.join(os.path.dirname(__file__),
        "../../../../main/host_filter.wdl"),
        "RunLZW"
    )
