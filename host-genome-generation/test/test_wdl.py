import os
from test_util import WDLTestCase


class TestIndexGeneration(WDLTestCase):
    """Tests the RunValidateInput function"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "host_genome_generation.wdl")
    common_inputs = {
        "input_fasta": "fixtures/input.fasta",
        "host_name": "test",
        "ercc_fasta": "fixtures/ERCC.fasta",
        "ercc_gtf": "fixtures/ERCC.gtf",
    }

    def testIndexGeneration(self):
        res = self.run_miniwdl()
        assert res, res
