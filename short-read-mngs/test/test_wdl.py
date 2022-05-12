import os
import json
import hashlib
import yaml
import csv
from test_util import WDLTestCase
from subprocess import CalledProcessError


class TestRunValidate(WDLTestCase):
    """Tests the RunValidateInput function"""

    wdl = os.path.join(os.path.dirname(__file__), "..", "host_filter.wdl")

    @classmethod
    def setUpClass(self):
        args = ["max_input_fragments=1", "file_ext=fastq", "s3_wd_uri=''"]
        self.rv_args = args

    def testValidateWindows(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "windows1.fastq.gz")
        args = self.rv_args + [f"fastqs={fastqs_0}"]
        res = self.run_miniwdl(args, task="RunValidateInput")
        with open(res["outputs"]["RunValidateInput.valid_input1_fastq"]) as f:
            hash = hashlib.md5(f.read().encode("utf-8")).hexdigest()
        self.assertEqual(hash, "a410dd184a01187d9c7c1823f5fc353e")

    def testInvalidInput(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "host_filter", "test_RunValidateInput_invalid_char.fastq")
        args = self.rv_args + [f"fastqs={fastqs_0}"]

        with self.assertRaises(CalledProcessError) as ecm:
            self.run_miniwdl(args, task="RunValidateInput")
        miniwdl_error = json.loads(ecm.exception.output)
        with open(miniwdl_error["cause"]["stderr_file"]) as stderr:
            error_json = stderr.readlines()[-1]
            cause = json.loads(error_json.strip())["cause"]
            self.assertEqual(cause, "PARSE ERROR: not an ascii file. Line 4 contains non-ascii characters.")


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
            "star_genome=s3://czid-public-references/host_filter/ercc"
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


class TestAlign(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "non_host_alignment.wdl")
    with open(os.path.join(os.path.dirname(__file__), "local_test.yml")) as fh:
        common_inputs = yaml.safe_load(fh)

    @classmethod
    def setUpClass(self):
        self.m8_file = os.path.join(
            os.path.dirname(__file__),
            "non_host_alignment",
            "call_hits_inputs",
            "gsnap.m8",
        )
        self.common_args = {
            "lineage_db": "s3://czid-public-references/taxonomy/2021-01-22/taxid-lineages.db",
            "accession2taxid": "s3://czid-public-references/mini-database/alignment_indexes/"
            "2020-08-20-viral/viral_accessions2taxid.db",
            "taxon_blacklist": "s3://czid-public-references/taxonomy/2021-01-22/taxon_blacklist.txt",
            "deuterostome_db": "s3://czid-public-references/taxonomy/2021-01-22/deuterostome_taxids.txt",
            "duplicate_cluster_size": os.path.join(
                os.path.dirname(__file__), "duplicate_cluster_sizes.tsv"
            ),
            "s3_wd_uri": "",
        }

    def testMinimap2CallHits(self):
        args = dict(self.common_args)
        args.update({"m8_file": self.m8_file, "prefix": "minimap2"})
        res = self.run_miniwdl(task="RunCallHitsMinimap2", task_input=args)
        m8_file = res["outputs"]["RunCallHitsMinimap2.deduped_out_m8"]
        deduped = {}
        with open(m8_file) as f:
            rd = csv.reader(f, delimiter="\t")
            for row in rd:
                deduped[row[0]] = row[1]

        self.assertEqual(len(deduped.keys()), 22)
        self.assertEqual(
            deduped[
                "NC_007795.1_64__benchmark_lineage_93061_1280_1279_90964__s0000000966"
            ],
            "NC_004615.1",
        )

    def testDiamondCallHits(self):
        args = dict(self.common_args)
        args.update(
            {
                "m8_file": self.m8_file,
                "prefix": "diamond",
                "min_read_length": 0,
                "count_type": "NR",
            }
        )
        res = self.run_miniwdl(task="RunCallHitsDiamond", task_input=args)
        m8_file = res["outputs"]["RunCallHitsDiamond.deduped_out_m8"]
        deduped = {}
        with open(m8_file) as f:
            rd = csv.reader(f, delimiter="\t")
            for row in rd:
                deduped[row[0]] = row[1]
        self.assertEqual(len(deduped.keys()), 22)
