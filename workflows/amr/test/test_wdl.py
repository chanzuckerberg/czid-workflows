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
