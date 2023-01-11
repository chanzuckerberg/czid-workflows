import csv
import gzip
import json
import os
import tempfile
from typing import Dict
from test_util import WDLTestCase


class TestLongReadMNGS(WDLTestCase):
    """Tests the czid_long_read_mngs workflow"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    ref_bucket = "s3://czid-public-references"
    ercc_prefix = "host_filter/ercc/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime"
    alignment_indexes_prefix = "mini-database/alignment_indexes/2020-08-20-viral"
    common_inputs = {
        # this mode cuts down on memory usage for testing
        "guppy_basecaller_setting": "super",
        "minimap_host_db": os.path.join(ref_bucket, ercc_prefix, "ERCC.fasta"),
        "minimap_human_db": os.path.join(ref_bucket, ercc_prefix, "ERCC.fasta"),
        "minimap2_local_db_path": os.path.join(ref_bucket, "test/viral-alignment-indexes/viral_nt"),
        "diamond_local_db_path": os.path.join(ref_bucket, "test/viral-alignment-indexes/viral_nr"),
        "accession2taxid_db": os.path.join(ref_bucket, alignment_indexes_prefix, "viral_accessions2taxid.marisa"),
        "lineage_db": os.path.join(ref_bucket, "taxonomy/2021-01-22/taxid-lineages.db"),
        "taxon_blacklist": os.path.join(ref_bucket, "taxonomy/2021-01-22/taxon_blacklist.txt"),
        "deuterostome_db": os.path.join(ref_bucket, "taxonomy/2021-01-22/deuterostome_taxids.txt"),
        "nt_info_db": os.path.join(ref_bucket, "test/viral-alignment-indexes/viral_nt_info.marisa"),
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

    def _read_length_metrics_assertions(self, outputs: Dict[str, str]):
        with open(outputs["czid_long_read_mngs.read_length_metrics"]) as f:
            metrics = json.load(f)
            self.assertLessEqual(metrics["read_length_min"], metrics["read_length_max"])
            self.assertLessEqual(metrics["read_length_mean"], metrics["read_length_max"])
            self.assertLessEqual(metrics["read_length_median"], metrics["read_length_max"])

    def _unmapped_reads_assertions(self, outputs: Dict[str, str]):
        self.assertIn("czid_long_read_mngs.unmapped_reads", outputs)
        with open(outputs["czid_long_read_mngs.unmapped_reads"]) as f:
            unmapped = f.readlines()
            self.assertGreaterEqual(len(unmapped), 40)  # TODO: figure out why this changes

    def testLongReadMNGSZipped(self):
        res = self.run_miniwdl([f"input_fastq={os.path.join(os.path.dirname(__file__), 'test_files/test.fastq.gz')}"])
        outputs = res["outputs"]
        self.assertIn("czid_long_read_mngs.nt_deduped_m8", outputs)
        self.assertIn("czid_long_read_mngs.nr_deduped_m8", outputs)

        # test tally hits
        self._tallied_hits_assertions(outputs, "nt_tallied_hits")
        self._tallied_hits_assertions(outputs, "nr_tallied_hits")

        # test unmapped reads
        self._unmapped_reads_assertions(outputs)

    def testLongReadMNGS(self):
        input_path = os.path.join(os.path.dirname(__file__), 'test_files/test.fastq.gz')
        with tempfile.NamedTemporaryFile(suffix=".fastq") as f, gzip.open(input_path) as zipped_f:
            for line in zipped_f:
                f.write(line)
            f.seek(0)

            res = self.run_miniwdl([f"input_fastq={f.name}"])
            outputs = res["outputs"]
            self.assertIn("czid_long_read_mngs.nt_deduped_m8", outputs)
            self.assertIn("czid_long_read_mngs.nr_deduped_m8", outputs)

            # test tally hits
            self._tallied_hits_assertions(outputs, "nt_tallied_hits")
            self._tallied_hits_assertions(outputs, "nr_tallied_hits")

            # test unmapped reads
            self._unmapped_reads_assertions(outputs)
