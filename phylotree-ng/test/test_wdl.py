import os
# import sys
# import tempfile
# import json

from test_util import WDLTestCase

class TestPhylotree(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    samples_dir = os.path.join(os.path.dirname(__file__), "samples")
    common_inputs = {
        "samples": [
            {
                "sample_name": "test_sample1",
                "workflow_run_id": 12345,
                "contig_fasta": os.path.join(samples_dir, "test_sample1", "contigs.fasta"),
                "combined_contig_summary": os.path.join(samples_dir, "test_sample1", "combined_contig_summary.json"),
            },
            {
                "sample_name": "test_sample2",
                "workflow_run_id": 23456,
                "contig_fasta": os.path.join(samples_dir, "test_sample2", "contigs.fasta"),
                "combined_contig_summary": os.path.join(samples_dir, "test_sample2", "combined_contig_summary.json"),
            },
        ],
        "reference": {"taxon_id": 64320},
        "additional_references": [
            {"accession_id": "MH157204.1"},
            {"accession_id": "MH900227.1"}
        ],
        "superkingdom_name": "fixme"
    }

    def test_phylotree(self):
        res = self.run_miniwdl(args=[])
        self.assertEqual(res, {})
