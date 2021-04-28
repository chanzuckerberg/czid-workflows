import os
# import sys
# import tempfile
# import json

from test_util import WDLTestCase

class TestPhylotree(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    common_inputs = {
        "samples": [
            {"sample_name": "test_sample", "workflow_run_id": 12345}
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
