import csv
import os
from typing import Dict
from test_util import WDLTestCase


class TestLongReadMNGS(WDLTestCase):
    """Tests the czid_long_read_mngs workflow"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    ref_bucket = "s3://czid-public-references"
    ercc_prefix = "host_filter/ercc/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime"
    alignment_indexes_prefix = "mini-database/alignment_indexes/2020-08-20-viral"
    common_inputs = {
        # this mode cuts down on memory usage for testin
        "guppy_basecaller_setting": "super",
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

    def _tallied_hits_assertions(self, outputs: Dict[str, str], name: str):
        with open(outputs[f"czid_long_read_mngs.{name}"]) as f:
            rows = list(csv.reader(f))
            self.assertEqual(rows[0], ["taxid", "level", "total_sequence_length", "total_alignment_length"])
            self.assertGreater(len(rows), 1)
            prev = None
            for row in rows[1:]:
                self.assertRegex(row[0], r"\d+")
                self.assertIn(row[1], ["genus", "species"])
                self.assertRegex(row[2], r"\d+")
                self.assertRegex(row[3], r"\d+")
                if prev:
                    self.assertGreaterEqual(prev, int(row[3]))
                prev = int(row[3])

    def _unmapped_reads_assertions(self, outputs: Dict[str, str]):
        self.assertIn("czid_long_read_mngs.unmapped_reads", outputs)
        with open(outputs["czid_long_read_mngs.unmapped_reads"]) as f:
            unmapped = f.readlines()
            self.assertEqual(len(unmapped), 44)

    def testLongReadMNGS(self):
        res = self.run_miniwdl([])
        outputs = res["outputs"]
        self.assertIn("czid_long_read_mngs.nt_deduped_out_m8", outputs)
        self.assertIn("czid_long_read_mngs.nr_deduped_out_m8", outputs)

        # test tally hits
        self._tallied_hits_assertions(outputs, "nt_tallied_hits")
        self._tallied_hits_assertions(outputs, "nr_tallied_hits")

        # test unmapped reads
        self._unmapped_reads_assertions(outputs)
