import csv
import gzip
import os
from test_util import WDLTestCase


class TestLongReadMNGS(WDLTestCase):
    """Tests the czid_long_read_mngs workflow"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    common_inputs = {
        "s3_wd_uri": "",
        "input_fastq": os.path.join(os.path.dirname(__file__), "test_files/test.fastq"),
        "minimap_host_db": "s3://czid-public-references/host_filter/ercc/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/ERCC.fasta",
        "minimap_human_db": "s3://czid-public-references/host_filter/ercc/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/ERCC.fasta",
        "minimap2_local_db_path": "s3://czid-public-references/test/viral-alignment-indexes/viral_nt",
        "diamond_local_db_path": "s3://czid-public-references/test/viral-alignment-indexes/viral_nr",
        "accession2taxid_db": "s3://czid-public-references/mini-database/alignment_indexes/2020-08-20-viral/viral_accessions2taxid.marisa",
        "lineage_db": "s3://czid-public-references/taxonomy/2021-01-22/taxid-lineages.db",
        "taxon_blacklist": "s3://czid-public-references/taxonomy/2021-01-22/taxon_blacklist.txt",
        "deuterostome_db": "s3://czid-public-references/taxonomy/2021-01-22/deuterostome_taxids.txt",
    }

    def testLongReadMNGS(self):
        res = self.run_miniwdl([])
        outputs = res["outputs"]
        self.assertIn("czid_long_read_mngs.nt_deduped_out_m8", outputs)
