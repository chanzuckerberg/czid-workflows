import os
from test_util import WDLTestCase


class TestIndexGeneration(WDLTestCase):
    """Tests the RunValidateInput function"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    common_inputs = {
        "input_fasta": os.path.join(os.path.dirname(__file__), "fixtures/input.fasta"),
        "host_name": "test",
        "ercc_fasta": os.path.join(os.path.dirname(__file__), "fixtures/ERCC.fasta"),
        "ercc_gtf": os.path.join(os.path.dirname(__file__), "fixtures/ERCC.gtf"),
    }

    def testIndexGeneration(self):
        res = self.run_miniwdl()
        assert res, res
