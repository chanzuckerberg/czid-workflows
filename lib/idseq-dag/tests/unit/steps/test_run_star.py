import unittest

# Class under test
from idseq_dag.steps.run_star import PipelineStepRunStar

class TestPipelineStepRunStar(unittest.TestCase):
    '''Tests for `util/run_star.py`'''

    # @patch('idseq_dag.steps.run_alignment_remotely.server.ASGInstance')
    def test_extract_rid(self):
        '''Extract rid from header string'''
        param_list = [
            ("@K11111:222:XXXXXXXXX:3:4444:5555:00000/2\t00", "@K11111:222:XXXXXXXXX:3:4444:5555:00000"),
            ("@K11111:222:XXXXXXXXX:3:4444:5555:00000\t00", "@K11111:222:XXXXXXXXX:3:4444:5555:00000"),
            (">K11111:222:XXXXXXXXX:3:4444:5555:00000\t00", ">K11111:222:XXXXXXXXX:3:4444:5555:00000")
        ]
        for line, expected in param_list:
            with self.subTest(line=line, expected=expected):
                rid = PipelineStepRunStar.extract_rid(line)

                self.assertEqual(rid, expected)
