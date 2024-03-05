import csv
import gzip
import os
from test_util import WDLTestCase

import marisa_trie


class TestIndexGeneration(WDLTestCase):
    """Tests the RunValidateInput function"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "index-generation.wdl")
    common_inputs = {
        "index_name": "2020-04-20",
        "ncbi_server": "https://idseq-samples-test.s3-us-west-2.amazonaws.com/index-generation/inputs",
        "provided_nt": "https://idseq-samples-test.s3-us-west-2.amazonaws.com/index-generation/inputs/blast/nt",
        "provided_nr": "https://idseq-samples-test.s3-us-west-2.amazonaws.com/index-generation/inputs/blast/nr",
    }

    def testIndexGeneration(self):
        # res = self.run_miniwdl(["index_name=2020-04-20"])
        res = self.run_miniwdl()
        outputs = res["outputs"]
        with gzip.open(outputs["index_generation.versioned_taxid_lineages_csv"], "rt") as f:
            for row in csv.DictReader(f):
                self.assertEqual(row["version_start"], "2020-04-20")
                self.assertEqual(row["version_end"], "2020-04-20")
        self.assertTrue(outputs["index_generation.nt_loc_db"].endswith(
            "marisa"), os.path.basename(outputs["index_generation.nt_loc_db"]))

        trie = marisa_trie.RecordTrie("256pI")
        trie.load(outputs["index_generation.nt_info_db"])
        self.assertEqual(list(trie.items()), [
            ("kraken:taxid|1434323|NC_027349.1", (b'Escherichia phage HY01, complete genome', 11200)),
            ("kraken:taxid|2557553|NC_048759.1", (b'Serratia phage MTx, partial genome', 68621))
        ])
