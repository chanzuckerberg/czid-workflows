import unittest
from tempfile import TemporaryFile

from idseq_dag.util.parsing import BlastnOutput6Reader, BlastnOutput6Writer

from tests.unit.unittest_helpers import relative_file_path


class TestBlastn6Reader(unittest.TestCase):
    def test_reading(self):
        blastn_input_6 = relative_file_path(__file__, "m8-test/sample.m8")
        with open(blastn_input_6) as blastn_input_6_f:
            rows = list(BlastnOutput6Reader(blastn_input_6_f))
            self.assertEqual(len(rows), 50)
            self.assertEqual(rows[0]["pident"], 100.0)

    def test_writing(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6Writer(f)
            writer.writerow({
                "qseqid": "1",
                "sseqid": "2",
                "pident": 1.0,
                "length": 5,
                "mismatch": 10,
                "gapopen": 15,
                "qstart": 20,
                "qend": 25,
                "sstart": 30,
                "send": 35,
                "evalue": 0.1,
                "bitscore": 1.0,
            })

    def test_writing_error(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6Writer(f)
            self.assertRaises(Exception, lambda : writer.writerow({"bad key": 1}))

    def test_filtration(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5",
            "2	MK468611.1	135.0	126	0	0	1	126	8433	8308	1.1e-74	290.5",  # pident too high
            "3 	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1.1e-74	290.5",  # pident too low
            "4	MK468611.1	-0.25	126	0	0	1	126	8433	8308	NaN	290.5",  # NaN error
            "5	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1	290.5",  # error too high
        ]

        rows = list(BlastnOutput6Reader(blastn_input_6, filter_invalid=True))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["qseqid"], "1")

    def test_read_error(self):
        blastn_input_6 = [""]
        self.assertRaises(Exception, lambda : next(BlastnOutput6Reader(blastn_input_6)))
