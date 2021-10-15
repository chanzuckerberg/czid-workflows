import os
import json
import yaml
from test_util import WDLTestCase
import gzip
import tempfile
from subprocess import CalledProcessError


class TestConsensusGenomes(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    with open(os.path.join(os.path.dirname(__file__), "local_test.yml")) as fh:
        common_inputs = yaml.safe_load(fh)
    sc2_ref_fasta = "s3://idseq-public-references/consensus-genome/MN908947.3.fa"

    def test_sars_cov2_illumina_cg(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "sample_sars-cov-2_paired_r1.fastq.gz")
        fastqs_1 = os.path.join(os.path.dirname(__file__), "sample_sars-cov-2_paired_r2.fastq.gz")
        args = ["sample=test_sample", f"fastqs_0={fastqs_0}", f"fastqs_1={fastqs_1}", "technology=Illumina",
                f"ref_fasta={self.sc2_ref_fasta}"]
        res = self.run_miniwdl(args)
        outputs = res["outputs"]
        with open(outputs["consensus_genome.compute_stats_out_output_stats"]) as fh:
            output_stats = json.load(fh)
        self.assertEqual(output_stats["sample_name"], "test_sample")
        # TODO: track non-determinism
        self.assertGreater(output_stats["depth_avg"], 222)
        self.assertEqual(output_stats["total_reads"], 47108)
        self.assertEqual(output_stats["mapped_reads"], 47054)
        self.assertEqual(output_stats["mapped_paired"], 47035)
        self.assertEqual(output_stats["ercc_mapped_reads"], 0)
        self.assertEqual(output_stats["ref_snps"], 7)
        self.assertEqual(output_stats["ref_mnps"], 0)
        for output_name, output in outputs.items():
            if output_name in {"consensus_genome.minion_log", "consensus_genome.vadr_errors"}:
                continue
            if not isinstance(output, list):
                output = [output]
            for filename in output:
                self.assertGreater(os.path.getsize(filename), 0)

    def test_sars_cov2_ont_cg(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "Ct20K.fastq.gz")
        args = ["sample=test_sample", f"fastqs_0={fastqs_0}", "technology=ONT", f"ref_fasta={self.sc2_ref_fasta}"]
        res = self.run_miniwdl(args)
        outputs = res["outputs"]
        with open(outputs["consensus_genome.compute_stats_out_output_stats"]) as fh:
            output_stats = json.load(fh)
        self.assertEqual(output_stats["sample_name"], "test_sample")
        self.assertGreater(output_stats["depth_avg"], 6)
        self.assertGreater(output_stats["depth_frac_above_10x"], 0.2)
        self.assertEqual(output_stats["depth_frac_above_100x"], 0)
        self.assertEqual(output_stats["total_reads"], 912)
        self.assertEqual(output_stats["mapped_reads"], 689)
        self.assertEqual(output_stats["mapped_paired"], 0)
        self.assertNotIn("ercc_mapped_reads", output_stats)
        self.assertEqual(output_stats["ref_snps"], 2)
        self.assertEqual(output_stats["ref_mnps"], 0)
        self.assertEqual(output_stats["n_actg"], 2461)
        self.assertEqual(output_stats["n_missing"], 27442)
        self.assertEqual(output_stats["n_gap"], 0)
        self.assertEqual(output_stats["n_ambiguous"], 0)
        for output_name, output in outputs.items():
            if output_name in {"consensus_genome.quantify_erccs_out_ercc_out",
                               "consensus_genome.filter_reads_out_filtered_fastqs",
                               "consensus_genome.trim_reads_out_trimmed_fastqs",
                               "consensus_genome.vadr_errors"}:
                continue
            if not isinstance(output, list):
                output = [output]
            for filename in output:
                self.assertIsNotNone(filename, output_name)
                self.assertGreater(os.path.getsize(filename), 0)

    def test_sars_cov2_ont_cg_no_length_filter(self):
        """
        Ensures that the apply_length_filter=false option has an effect by truncating every other input read
        and asserting lower coverage in the output
        """
        fastqs_0 = os.path.join(os.path.dirname(__file__), "Ct20K.fastq.gz")
        tf = tempfile.NamedTemporaryFile(prefix=__name__, suffix=".fastq.gz")
        with gzip.open(fastqs_0, mode="r") as in_fh, gzip.open(tf.name, mode="w") as out_fh:
            lineno = 0
            while True:
                header, seq, sep, qual = in_fh.readline(), in_fh.readline(), in_fh.readline(), in_fh.readline()
                if not header:
                    break
                if lineno % 2 == 0:
                    seq, qual = seq[:100] + b"\n", qual[:100] + b"\n"
                out_fh.write(header + seq + sep + qual)
                lineno += 1
        args = ["sample=test_sample", f"fastqs_0={tf.name}", "technology=ONT", f"ref_fasta={self.sc2_ref_fasta}",
                "apply_length_filter=false"]
        res = self.run_miniwdl(args)
        outputs = res["outputs"]
        with open(outputs["consensus_genome.compute_stats_out_output_stats"]) as fh:
            output_stats = json.load(fh)
        self.assertEqual(output_stats["sample_name"], "test_sample")
        self.assertGreater(output_stats["depth_avg"], 3)
        self.assertLess(output_stats["depth_avg"], 4)

    def test_general_cg(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "SRR11741455_65054_nh_R1.fastq.gz")
        fastqs_1 = os.path.join(os.path.dirname(__file__), "SRR11741455_65054_nh_R2.fastq.gz")
        args = ["sample=test_sample", f"fastqs_0={fastqs_0}", f"fastqs_1={fastqs_1}", "technology=Illumina",
                "filter_reads=false", "ref_accession_id=MF965207.1",
                "primer_bed=s3://idseq-public-references/consensus-genome/na_primers.bed"]
        res = self.run_miniwdl(args)
        for output_name, output in res["outputs"].items():
            if isinstance(output, str):
                self.assertGreater(os.path.getsize(output), 0)
            elif output:
                for path in output:
                    self.assertGreater(os.path.getsize(path), 0)
        with open(res["outputs"]["consensus_genome.compute_stats_out_output_stats"]) as fh:
            output_stats = json.load(fh)
        self.assertEqual(output_stats["sample_name"], "test_sample")
        self.assertGreater(output_stats["depth_avg"], 510)
        self.assertGreaterEqual(output_stats["depth_frac_above_10x"], 0.99)
        self.assertGreaterEqual(output_stats["depth_frac_above_100x"], 0.99)
        self.assertEqual(output_stats["total_reads"], 100000)
        self.assertEqual(output_stats["mapped_reads"], 55688)
        self.assertEqual(output_stats["mapped_paired"], 55680)
        self.assertEqual(output_stats["ercc_mapped_reads"], 0)
        self.assertEqual(output_stats["ref_snps"], 0)
        self.assertEqual(output_stats["ref_mnps"], 0)
        # TODO: address this non-determinism
        self.assertIn(output_stats["n_actg"], [15317])
        self.assertIn(output_stats["n_missing"], [0, 1])
        self.assertEqual(output_stats["n_gap"], 0)
        self.assertEqual(output_stats["n_ambiguous"], 1)

        self.assertEqual(res["outputs"]["consensus_genome.vadr_alerts_out"], None)
        self.assertEqual(res["outputs"]["consensus_genome.vadr_quality_out"], None)
        self.assertEqual(res["outputs"]["consensus_genome.vadr_errors"], None)

        args.append(f"ref_fasta={self.sc2_ref_fasta}")
        with self.assertRaises(CalledProcessError) as ecm:
            self.run_miniwdl(args)
        self.assertRunFailed(ecm, task="MakeConsensus", error="InsufficientReadsError",
                             cause="No reads after MakeConsensus")

    def test_length_filter_midnight_primers(self):
        """
        Test that the length filters are properly set for midnight primers
        """
        fastqs_0 = os.path.join(os.path.dirname(__file__), "blank.fastq.gz")
        args = [
            "sample=test_sample",
            f"fastqs_0={fastqs_0}",
            "technology=ONT",
            f"ref_fasta={self.sc2_ref_fasta}",
            "primer_set=nCoV-2019/V1200",
        ]
        with self.assertRaises(CalledProcessError) as ecm:
            self.run_miniwdl(args)
        miniwdl_error = json.loads(ecm.exception.output)
        with open(
            os.path.join(miniwdl_error["dir"], "call-ApplyLengthFilter", "inputs.json")
        ) as f:
            apply_length_filter_inputs = json.load(f)
        self.assertEqual(apply_length_filter_inputs["min_length"], 250)
        self.assertEqual(apply_length_filter_inputs["max_length"], 1500)

    def test_sars_cov2_ont_cg_no_reads(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "blank.fastq.gz")
        args = ["sample=test_sample", f"fastqs_0={fastqs_0}", "technology=ONT", f"ref_fasta={self.sc2_ref_fasta}"]
        with self.assertRaises(CalledProcessError) as ecm:
            self.run_miniwdl(args)
        self.assertRunFailed(
            ecm,
            task="RemoveHost",
            error="InsufficientReadsError",
            cause="No reads after RemoveHost"
        )
