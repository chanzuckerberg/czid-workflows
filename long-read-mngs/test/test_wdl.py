import csv
import gzip
import os
from test_util import WDLTestCase


class TestLongReadMNGS(WDLTestCase):
    """Tests the czid_long_read_mngs workflow"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    common_inputs = {
        "s3_wd_uri": "",
        "input_fastq": os.path.join(os.path.dirname(__file__), "test_files/dummy.fastq"),
        "minimap_host_db": os.path.join(os.path.dirname(__file__), "test_files/human.fasta"),
        "minimap_human_db": os.path.join(os.path.dirname(__file__), "test_files/human.fasta"),
        "minimap2_local_db_path": os.path.join(os.path.dirname(__file__), "test_files/viral_nt_minimap2_index"),
        "accession2taxid_db": "s3://czid-public-references/mini-database/alignment_indexes/2020-08-20-viral/viral_accessions2taxid.marisa",
        "lineage_db": "s3://czid-public-references/taxonomy/2021-01-22/taxid-lineages.db",
        "taxon_blacklist": "s3://czid-public-references/taxonomy/2021-01-22/taxon_blacklist.txt",
        "deuterostome_db": "s3://czid-public-references/taxonomy/2021-01-22/deuterostome_taxids.txt",
    }

    def testLongReadMNGS(self):
        res = self.run_miniwdl([])
        outputs = res["outputs"]
        with gzip.open(outputs["czid_long_read_mngs.foo"], "rt") as f:
            pass
