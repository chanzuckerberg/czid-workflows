import unittest
import re
import os
import subprocess
from unittest.mock import patch
from tests.unit.unittest_helpers import relative_file_path, file_contents
from idseq_dag.steps.generate_phylo_tree import PipelineStepGeneratePhyloTree
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

ASSEMBLY_SUMMARY_FILE = relative_file_path(__file__, "../../../examples/fixtures/assembly_summary.txt")
EXAMPLE_VCF_FILE = relative_file_path(__file__, "../../../examples/fixtures/example.vcf")
TMP_VCF_OUT_FILE = "/tmp/tmp_generatephylotree_testcase.vcf"


class GeneratePyhloTreeTestCase(unittest.TestCase):
    '''Tests for idseq_dag/steps/generate_pyhlo_tree.py module'''

    def test_get_taxid_genomes(self):
        results = PipelineStepGeneratePhyloTree.get_taxid_genomes(ASSEMBLY_SUMMARY_FILE, 10298, 2)

        self.assertEqual(
            results,
            [
                'GCA_000859985.2\t10298\t10298\tHuman alphaherpesvirus 1\tftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/859/985/GCA_000859985.2_ViralProj15217',
                'GCA_003052245.1\t10298\t10298\tHuman alphaherpesvirus 1\tftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/052/245/GCA_003052245.1_ASM305224v1'
            ]
        )

    @staticmethod
    def _ncbi_output_stub(cmd):
        if isinstance(cmd, command_patterns.CommandPattern):
            if cmd.cmd == "curl" and cmd.args == ["https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=AY821803.1&usehistory=y"]:
                return r'''<?xml version="1.0" encoding="UTF-8" ?> <!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd"> <eSearchResult><Count>1</Count><RetMax>1</RetMax><RetStart>0</RetStart><QueryKey>1</QueryKey><WebEnv>NCID_1_252904980_130.14.22.76_9001_1564705514_734685590_0MetA0_S_MegaStore</WebEnv><IdList> <Id>60686968</Id> </IdList><TranslationSet/><QueryTranslation/></eSearchResult>'''
            if cmd.cmd == "curl" and cmd.args == ["https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=1&WebEnv=NCID_1_252904980_130.14.22.76_9001_1564705514_734685590_0MetA0_S_MegaStore&rettype=gb&retmode=xml"]:
                return '<?xml version="1.0" encoding="UTF-8"  ?>\n<!DOCTYPE GBSet PUBLIC "-//NCBI//NCBI GBSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd">\n<GBSet>\n<GBSeq_locus>AY821803</GBSeq_locus>\n<!-- fake content -->\n</GBSet>\n'
            raise Exception("Not stubbed", cmd.as_dict())
        raise Exception("Not stubbed", cmd)


    def test_fetch_ncbi(self):
        with patch("idseq_dag.steps.generate_phylo_tree.command.execute_with_output", side_effect=self._ncbi_output_stub):
            result = PipelineStepGeneratePhyloTree.fetch_ncbi("AY821803.1")

            expected_string_in_result = '<GBSeq_locus>AY821803</GBSeq_locus>'
            self.assertEqual(result['search_url'], "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=AY821803.1&usehistory=y")
            self.assertTrue(expected_string_in_result in result['genbank_xml'], f"result should contain {expected_string_in_result}, but it contains {result['genbank_xml'][:500]}...") #, "^" + re.escape(expected_result_prefix) + ".*")


    def test_name_samples_vcf(self):
        with self.subTest("_vcf_new_column_description"):
            new_column_description = PipelineStepGeneratePhyloTree._vcf_new_column_description(
                input_file=EXAMPLE_VCF_FILE,
                sample_names_by_run_ids={'5678': r'sample_/name&\_1'}
            )

            self.assertEqual(new_column_description, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpipeline_run_1234\tnon_numeric_column\tsample_/name&\\_1")

        with self.subTest("_vcf_replace_column_description"):
            # setup
            os.system(f"rm -rf {TMP_VCF_OUT_FILE}")

            # test
            PipelineStepGeneratePhyloTree._vcf_replace_column_description(
                input_file=EXAMPLE_VCF_FILE,
                output_file=TMP_VCF_OUT_FILE,
                new_column_description="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpipeline_run_1234\tnon_numeric_column\tsample_/name&\\_1"
            )

            # assert
            self.assertRegex(file_contents(TMP_VCF_OUT_FILE), "#CHROM.*pipeline_run_1234\tnon_numeric_column\tsample_/name&\\\\_1.*")

            # teardown
            os.system(f"rm -rf {TMP_VCF_OUT_FILE}")
