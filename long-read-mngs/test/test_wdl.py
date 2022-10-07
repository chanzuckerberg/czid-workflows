import csv
import os
from test_util import WDLTestCase


class TestLongReadMNGS(WDLTestCase):
    """Tests the czid_long_read_mngs workflow"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    ref_bucket = "s3://czid-public-references"
    ercc_prefix = "host_filter/ercc/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime"
    alignment_indexes_prefix = "mini-database/alignment_indexes/2020-08-20-viral"
    common_inputs = {
        "s3_wd_uri": "",
        "input_fastq": os.path.join(os.path.dirname(__file__), "test_files/test.fastq"),
        "minimap_host_db": os.path.join(ref_bucket, ercc_prefix, "ERCC.fasta"),
        "minimap_human_db": os.path.join(ref_bucket, ercc_prefix, "ERCC.fasta"),
        "minimap2_local_db_path": os.path.join(ref_bucket, "test/viral-alignment-indexes/viral_nt"),
        "diamond_local_db_path": os.path.join(ref_bucket, "test/viral-alignment-indexes/viral_nr"),
        "accession2taxid_db": os.path.join(ref_bucket, alignment_indexes_prefix, "viral_accessions2taxid.marisa"),
        "lineage_db": os.path.join(ref_bucket, "taxonomy/2021-01-22/taxid-lineages.db"),
        "taxon_blacklist": os.path.join(ref_bucket, "taxonomy/2021-01-22/taxon_blacklist.txt"),
        "deuterostome_db": os.path.join(ref_bucket, "taxonomy/2021-01-22/deuterostome_taxids.txt"),
    }

    def testLongReadMNGS(self):
        res = self.run_miniwdl([])
        outputs = res["outputs"]
        self.assertIn("czid_long_read_mngs.nt_deduped_out_m8", outputs)
        with open(outputs["czid_long_read_mngs.tallied_hits"]) as f:
            rows = list(csv.reader(f))
            self.assertEqual(rows[0], ["final_taxid", "aln_len", "seq_len"])
            self.assertGreater(len(rows), 1)

        with open(outputs["czid_long_read_mngs.genus_tallied_hits"]) as f:
            rows = list(csv.reader(f))
            self.assertEqual(rows[0], ["genus", "aln_len", "seq_len"])
            self.assertGreater(len(rows), 1)
