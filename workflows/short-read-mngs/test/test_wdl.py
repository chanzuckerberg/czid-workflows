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
            "lineage_db": "s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/"
            "taxid-lineages.marisa",
            "accession2taxid": "s3://czid-public-references/mini-database/alignment_indexes/"
            "2020-08-20-viral/viral_accessions2taxid.marisa",
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
