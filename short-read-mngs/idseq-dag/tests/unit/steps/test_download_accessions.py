import re
import unittest
import tempfile
from multiprocessing import RLock


# Class under test
from idseq_dag.steps.download_accessions import PipelineStepDownloadAccessions

FAKE_SEQUENCE_DATA = (
    'GGREKKKLLKSQKDGSAKDRHNPRAFSVANIVRTQRNVQRNLDRAQKKEYVPLSDRRAARVEEGPPSLVAVVGPPGVGKS\n'
    'TLIRSLVKLYTNHNLTNPTGPITVCTSQTKRITFLECPNTPTAMLDVAKIADLVLLCVDAKFGFEMETFEFLNMMQTHGF\n'
    'KLYAPLSNVGAVSFDKDAVYIDIGRANYTKRENLDLPREEKEEGEGSDDASSLQDVKEGVDQKMQYSSLRLFKGSKAVKA\n'
)
class TestDownloadAccessions(unittest.TestCase):
    def setUp(self):
        self.step = PipelineStepDownloadAccessions(
            name="rapsearch2_out",
            input_files=[
                [
                    "rapsearch2.m8",
                    "rapsearch2.deduped.m8",
                    "rapsearch2.hitsummary.tab",
                    "rapsearch2_counts.json"
                ]
            ],
            output_files=["assembly/nr.refseq.fasta"],
            output_dir_local=tempfile.TemporaryDirectory().name,
            output_dir_s3="s3://dummy_bucket",
            ref_dir_local=tempfile.TemporaryDirectory().name,
            additional_files={
                "lineage_db": "s3://idseq-public-references/taxonomy/2018-12-01/taxid-lineages.sqlite3",
                "loc_db": "s3://idseq-public-references/alignment_data/2018-12-01/nr_loc.sqlite3"
            },
            additional_attributes={
                "db":  "s3://idseq-public-references/alignment_data/2018-12-01/nr",
                "db_type": "nr"
            },
        )


    def test_fix_ncbi_record(self):
        '''Run command remotely with success'''
        parameterized_tests = [
            [
                'WHEN one or more description from a multi-header start with a comma, THEN return the header with comma and spaces from the beginning of each description removed',
                '>XP_002289390.1 , partial [Thalassiosira pseudonana CCMP1335]\x01EED92927.1 ,  Conserved Hypothetical Protein, partial [Thalassiosira pseudonana CCMP1335]\n' + FAKE_SEQUENCE_DATA,
                '>XP_002289390.1 partial [Thalassiosira pseudonana CCMP1335]\x01EED92927.1 Conserved Hypothetical Protein, partial [Thalassiosira pseudonana CCMP1335]\n' + FAKE_SEQUENCE_DATA
            ],
            [
                'WHEN entry is a single-header and the description starts with a comma THEN return the header with the comma and spaces from the beginning of the description removed',
                '>OGG82825.1 , replicative DNA helicase [Candidatus Kaiserbacteria bacterium RIFCSPLOWO2_02_FULL_54_13]\n' + FAKE_SEQUENCE_DATA,
                '>OGG82825.1 replicative DNA helicase [Candidatus Kaiserbacteria bacterium RIFCSPLOWO2_02_FULL_54_13]\n' + FAKE_SEQUENCE_DATA,
            ],
            [
                "WHEN none of the descriptions from a multi-header start with a comma THEN return the original headers",
                '>XP_002289390.1 partial [Thalassiosira pseudonana CCMP1335]\x01EED92927.1 Conserved Hypothetical Protein, partial [Thalassiosira pseudonana CCMP1335]\n' + FAKE_SEQUENCE_DATA,
                '>XP_002289390.1 partial [Thalassiosira pseudonana CCMP1335]\x01EED92927.1 Conserved Hypothetical Protein, partial [Thalassiosira pseudonana CCMP1335]\n' + FAKE_SEQUENCE_DATA
            ],
            [
                "WHEN entry is a single-header and description doesn't start with a comma THEN return the original header",
                '>OGG82825.1 replicative DNA helicase [Candidatus Kaiserbacteria bacterium RIFCSPLOWO2_02_FULL_54_13]\n' + FAKE_SEQUENCE_DATA,
                '>OGG82825.1 replicative DNA helicase [Candidatus Kaiserbacteria bacterium RIFCSPLOWO2_02_FULL_54_13]\n' + FAKE_SEQUENCE_DATA,
            ],
            [
                "WHEN entry is a single-header that has no description THEN return the original header",
                '>OGG82825.1\n' + FAKE_SEQUENCE_DATA,
                '>OGG82825.1\n' + FAKE_SEQUENCE_DATA,
            ],
        ]

        for test_name, accession_data, expected_result in parameterized_tests:
            with self.subTest(test_name):
                result = self.step._fix_ncbi_record(accession_data)

                self.assertEqual(result, expected_result)
        
