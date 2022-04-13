import os
from test_util import WDLTestCase


class TestIndexGeneration(WDLTestCase):
    """Tests the RunValidateInput function"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "index_generation.wdl")
    common_inputs = {
        "ncbi_server": "https://idseq-samples-test.s3-us-west-2.amazonaws.com/index-generation/inputs",
    }

    def testIndexGeneration(self):
        res = self.run_miniwdl(["index_name=2020-04-20"])
        assert res
