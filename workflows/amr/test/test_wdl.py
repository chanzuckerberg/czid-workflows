import csv
import os
from test_util import WDLTestCase


def relpath(*args):
    """ helper to get a filepath relative to the current file """
    return os.path.join(os.path.dirname(__file__), *args)


class TestAMR(WDLTestCase):
    """Tests AMR"""

    wdl = relpath("..", "run.wdl")
    CARD = "s3://czid-public-references/test/AMRv2/card.json"

    def testRgiMain(self):
        inputs = {"contigs_fa": relpath("contigs.fasta"), "card_json": self.CARD}
        res = self.run_miniwdl(task="RunRgiMain", task_input=inputs)
        with open(res["outputs"]["RunRgiMain.main_amr_results"]) as main_results:
            lastline = main_results.read().splitlines()[-1]
        main = lastline.split("\t")
        self.assertEqual(main[2], "75")
        self.assertEqual(main[3], "260")

    def testRgiEmptyMain(self):
        inputs = {"contigs_fa": relpath("contigs_failed.fasta"), "card_json": self.CARD}
        res = self.run_miniwdl(task="RunRgiMain", task_input=inputs)
        with open(res["outputs"]["RunRgiMain.main_amr_results"]) as main_results:
            self.assertEqual(len(main_results.read().splitlines()), 1)

    def testRgiBwtKma(self):
        inputs = {
            "non_host_reads_fa": [
                relpath("gsnap_filter_1.fa"),
                relpath("gsnap_filter_2.fa"),
            ],
            "card_json": self.CARD,
        }
        res = self.run_miniwdl(task="RunRgiBwtKma", task_input=inputs)
        with open(res["outputs"]["RunRgiBwtKma.kma_amr_results"]) as kma_results:
            lastline = kma_results.read().splitlines()[-1]
        kma = lastline.split("\t")
        self.assertEqual(kma[8], "Elizabethkingia anophelis")

    def testRunResultsPerSample(self):
        inputs = {
            "sample_name": "sample",
            "gene_coverage": relpath("RunResultsPerSample", "gene_coverage.tsv"),
            "card_ontology": relpath("RunResultsPerSample", "ontology.json"),
            "kma_species_output": relpath(
                "RunResultsPerSample", "sr_species_report_61mer_analysis.gene.txt"
            ),
            "kma_output": relpath(
                "RunResultsPerSample", "sr_amr_report.allele_mapping_data.txt"
            ),
            "main_species_output": relpath(
                "RunResultsPerSample",
                "contig_species_report_61mer_analysis_rgi_summary.txt",
            ),
            "main_output": relpath("RunResultsPerSample", "contig_amr_report.txt"),
        }
        res = self.run_miniwdl(task="RunResultsPerSample", task_input=inputs)
        with open(
            res["outputs"]["RunResultsPerSample.synthesized_report"]
        ) as synthesized_report:
            primary_amr_report = synthesized_report.read().splitlines()

        self.assertEqual(len(primary_amr_report), 6)
        gene_names = [i.split("\t")[0] for i in primary_amr_report]
        for expected_gene_name in [
            "ArnT",
            "Staphylococcus aureus LmrS",
            "fosA5",
            "norC",
            "smeB",
        ]:
            self.assertTrue(expected_gene_name in gene_names)

    def testRunRedup(self):
        inputs = {
            "host_filtered_reads": [
                relpath("RunRedup", "host_filter_1.fastq"),
                relpath("RunRedup", "host_filter_2.fastq")
            ],
            "subsampled_reads": [relpath("RunRedup", "subsampled_1.fa"), relpath("RunRedup", "subsampled_2.fa")],
            "clusters": relpath("RunRedup", "clusters.csv"),
            "cluster_sizes": relpath("RunRedup", "duplicate_cluster_sizes.tsv"),
        }

        # Get all input read ids
        with open(relpath("RunRedup", "clusters.csv"), 'r') as clusters_csv:
            tsvreader = csv.reader(clusters_csv)
            cluster_read_ids = [row for row in tsvreader]

        # Get subsampled read ids
        subsampled_read_ids = set()
        with open(relpath("RunRedup", "subsampled_1.fa"), 'r') as subsampled_1_fa:
            for line in subsampled_1_fa.readlines():
                if line.startswith(">"):
                    subsampled_read_ids.add(line.lstrip(">").rstrip("\n"))

        target_read_ids = set()
        excluded_read_ids = set()
        for rep_read_id, read_id in cluster_read_ids:
            if rep_read_id in subsampled_read_ids:
                target_read_ids.add(read_id)
            else:
                excluded_read_ids.add(read_id)

        # Run task and collect output ids
        res = self.run_miniwdl(task="RunRedup", task_input=inputs)

        output_1 = dict()
        with open(res["outputs"]["RunRedup.redups_fa"][0]) as redups_1_fa:
            reads = []
            this_read = []
            for line in redups_1_fa.readlines():
                if line.startswith(">"):
                    reads.append(this_read)
                    this_read = [line.lstrip(">").rstrip("\n")]
                else:
                    this_read.append(line.rstrip("\n"))
            reads.append(this_read)
            output_1 = {read[0]: "".join(read[1:]) for read in reads if len(read) > 0}

        output_2 = dict()
        with open(res["outputs"]["RunRedup.redups_fa"][1]) as redups_2_fa:
            reads = []
            this_read = []
            for line in redups_2_fa.readlines():
                if line.startswith(">"):
                    reads.append(this_read)
                    this_read = [line.lstrip(">").rstrip("\n")]
                else:
                    this_read.append(line.rstrip("\n"))
            reads.append(this_read)
            output_2 = {read[0]: "".join(read[1:]) for read in reads if len(read) > 0}

        # Check ids
        output_ids_1 = output_1.keys()
        assert len(set(output_ids_1)) == len(output_ids_1)
        output_ids_1 = set(output_ids_1)

        output_ids_2 = output_2.keys()
        assert len(set(output_ids_2)) == len(output_ids_2)
        output_ids_2 = set(output_ids_2)

        assert output_ids_1 == output_ids_2
        output_ids = output_ids_1 & output_ids_2

        for read_id in target_read_ids:
            assert read_id in output_ids
        for read_id in excluded_read_ids:
            assert read_id not in output_ids

        # Check sequences
        host_filter_1 = dict()
        with open(relpath("RunRedup", "host_filter_1.fastq"), "r") as host_filter_1_fq:
            lines = host_filter_1_fq.readlines()
            for index in range(0, len(lines), 4):
                read_id = lines[index].lstrip("@").rstrip("\n")
                read_sequence = lines[index + 1].rstrip("\n")
                host_filter_1[read_id] = read_sequence

        host_filter_2 = dict()
        with open(relpath("RunRedup", "host_filter_2.fastq"), "r") as host_filter_2_fq:
            lines = host_filter_2_fq.readlines()
            for index in range(0, len(lines), 4):
                read_id = lines[index].lstrip("@").rstrip("\n")
                read_sequence = lines[index + 1].rstrip("\n")
                host_filter_2[read_id] = read_sequence

        for read_id, read_sequence in output_1.items():
            assert read_id in host_filter_1
            assert read_sequence == host_filter_1[read_id]
        for read_id, read_sequence in output_2.items():
            assert read_id in host_filter_2
            assert read_sequence == host_filter_2[read_id]

    # def testRunSpadesAssemblyFailedNonUniformCoverage(self):
    #     inputs = {
    #         "reduplicated_reads": [
    #             relpath("RunSpades", "subsampled_1.fa"),
    #             relpath("RunSpades", "subsampled_2.fa")
    #             ],
    #         "min_contig_length": 100,
    #     }

    #     res = self.run_miniwdl(
    #         task="RunSpades",
    #         task_input=inputs,
    #         docker_image_id="ghcr.io/chanzuckerberg/czid-workflows/czid-short-read-mngs-public"
    #     )

    #     with open(res["outputs"]["RunSpades.contigs"]) as contigs_fa:
    #         lines = contigs_fa.readlines()
    #         assert lines[0].startswith(";ASSEMBLY FAILED")

    # def testRunSpadesAssemblyFailedExitCodeZero(self):
    #     inputs = {
    #         "reduplicated_reads": [
    #             relpath("RunSpades", "redups_1.fa"),
    #             relpath("RunSpades", "redups_2.fa")
    #             ],
    #         "min_contig_length": 100,
    #     }

    #     res = self.run_miniwdl(
    #         task="RunSpades",
    #         task_input=inputs,
    #         docker_image_id="ghcr.io/chanzuckerberg/czid-workflows/czid-short-read-mngs-public"
    #     )

    #     with open(res["outputs"]["RunSpades.contigs"]) as contigs_fa:
    #         lines = contigs_fa.readlines()
    #         assert lines[0].startswith(";ASSEMBLY FAILED")
