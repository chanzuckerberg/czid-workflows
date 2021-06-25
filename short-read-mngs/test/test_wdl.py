import os
import json
import hashlib
import yaml
from test_util import WDLTestCase


class TestSTAR(WDLTestCase):
    """Tests the RunSTAR function
    the inputs are minimal, with only 100 reads
    should only add ~1 min to testing time
    """

    wdl = os.path.join(os.path.dirname(__file__), "..", "host_filter.wdl")
    with open(os.path.join(os.path.dirname(__file__), "local_test.yml")) as fh:
        common_inputs = yaml.safe_load(fh)
    star_args = None

    @classmethod
    def setUpClass(self):
        fastqs_0 = os.path.join(
            os.path.dirname(__file__),
            "host_filter",
            "star_inputs",
            "valid_input1.fastq",
        )
        fastqs_1 = os.path.join(
            os.path.dirname(__file__),
            "host_filter",
            "star_inputs",
            "valid_input2.fastq",
        )
        summary_json = os.path.join(
            os.path.dirname(__file__),
            "host_filter",
            "star_inputs",
            "validate_input_summary.json",
        )
        args = [
            "s3_wd_uri=''",
            f"validate_input_summary_json={summary_json}",
            f"valid_input_fastq={fastqs_0}",
            f"valid_input_fastq={fastqs_1}",
            "star_genome=s3://idseq-public-references/host_filter/ercc"
            "/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/STAR_genome.tar",
        ]
        self.star_args = args

    def test_star(self):
        """test the basic star parameters"""
        args = self.star_args + ["nucleotide_type=DNA", "host_genome=human"]
        res = self.run_miniwdl(args, task="RunStar")
        with open(res["outputs"]["RunStar.output_read_count"]) as f:
            count = json.load(f)

        self.assertEqual(count["star_out"], 100)
        with open(res["outputs"]["RunStar.unmapped1_fastq"]) as f:
            hash = hashlib.md5(f.read().encode("utf-8")).hexdigest()
            self.assertEqual(hash, "c4d71e1b9b01734f7c3d300a7eac327a")
        with open(res["outputs"]["RunStar.unmapped2_fastq"]) as f:
            hash = hashlib.md5(f.read().encode("utf-8")).hexdigest()
            self.assertEqual(hash, "6b46fe79bf089c8b3f6377fab34b9744")

    def test_star_rna(self):
        """test the nucleotide_type of RNA works, should run STAR with TranscriptomeSAM"""
        args = self.star_args + ["nucleotide_type=RNA", "host_genome=human"]
        res = self.run_miniwdl(args, task="RunStar")
        with open(res["outputs"]["RunStar.output_read_count"]) as f:
            count = json.load(f)
        self.assertEqual(count["star_out"], 100)
        self.assertIn("TranscriptomeSAM", res["outputs"]["RunStar.step_description_md"])

    def test_star_nonhuman(self):
        """test that there is no output BAM file if the host is non-human"""
        args = self.star_args + ["nucleotide_type=DNA", "host_genome=pig"]
        res = self.run_miniwdl(args, task="RunStar")

        with open(res["outputs"]["RunStar.output_read_count"]) as f:
            count = json.load(f)
        self.assertEqual(count["star_out"], 100)
        self.assertIsNone(res["outputs"]["RunStar.aligned_file"])

    def test_starlong(self):
        """tests that STARLong runs if # of reads with length > 500 is >1
        the validation input has been modified, but there are no actual long reads
        """

        args = self.star_args + ["nucleotide_type=DNA", "host_genome=human"]
        args[1] = args[1].replace(".json", "_long.json")
        res = self.run_miniwdl(args, task="RunStar")
        with open(res["outputs"]["RunStar.output_read_count"]) as f:
            count = json.load(f)
        self.assertEqual(count["star_out"], 100)
