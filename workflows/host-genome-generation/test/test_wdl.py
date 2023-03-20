import os
from test_util import WDLTestCase


class TestIndexGeneration(WDLTestCase):
    """Tests the RunValidateInput function"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "host_genome_generation.wdl")
    common_inputs = {
        "genome_name": "test",
        "genome_fasta_gz": os.path.join(os.path.dirname(__file__), "fixtures/input.fa.gz"),
        "ERCC_fasta_gz": os.path.join(os.path.dirname(__file__), "fixtures/ERCC.fa.gz"),
    }

    def testIndexGeneration(self):
        res = self.run_miniwdl()
        assert res, res
