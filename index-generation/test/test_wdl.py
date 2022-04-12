import os
from test_util import WDLTestCase


class TestIndexGeneration(WDLTestCase):
    """Tests the RunValidateInput function"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "index_generation.wdl")

    @classmethod
    def setUpClass(self):
        args = ["index_name=2020-04-20"]
        self.rv_args = args

    def testIndexGeneration(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "windows1.fastq.gz")
        args = self.rv_args + [f"fastqs={fastqs_0}"]
        res = self.run_miniwdl(args, task="RunValidateInput")
        assert res
