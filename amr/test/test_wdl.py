import os
from test_util import WDLTestCase


class TestAMR(WDLTestCase):
    """Tests AMR"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")

    def testRgiMain(self):
        inputs = {
            "contigs": os.path.join(os.path.dirname(__file__), "contigs.fasta"),
            "card_json": "s3://czid-public-references/test/AMRv2/card.json",
        }
        res = self.run_miniwdl(task="RunRgiMain", task_input=inputs)
        with open(res["outputs"]["RunRgiMain.main_amr_results"]) as main_results:
            lastline = main_results.read().splitlines()[-1]
        main = lastline.split("\t")
        self.assertEqual(main[2], '3')
        self.assertEqual(main[3], '281')
    

    def testRgiMain(self):
        inputs = {
            "non_host_reads": [os.path.join(os.path.dirname(__file__), "gsnap_filter_1.fa"),
                               os.path.join(os.path.dirname(__file__), "gsnap_filter_2.fa")],
            "card_json": "s3://czid-public-references/test/AMRv2/card.json",
        }
        res = self.run_miniwdl(task="RunRgiBwtKma", task_input=inputs)
        with open(res["outputs"]["RunRgiBwtKma.kma_amr_results"]) as kma_results:
            lastline = kma_results.read().splitlines()[-1] 
        kma = lastline.split("\t")
        self.assertEqual(kma[8], 'Elizabethkingia anophelis')
        
