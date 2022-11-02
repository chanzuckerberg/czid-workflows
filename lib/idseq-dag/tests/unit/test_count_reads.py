import os
import sys
import tempfile
import unittest

from idseq_dag.util.count import count_reads
from idseq_dag.exceptions import InvalidInputFileError

class TestCountReads(unittest.TestCase):
    def test_count_reads(self):
        expect_reads = {
            os.path.join(os.path.dirname(__file__), "fixtures", "reads.fasta"): 402,
            os.path.join(os.path.dirname(__file__), "fixtures", "reads.fasta.gz"): 402,
            os.path.join(os.path.dirname(__file__), "fixtures", "reads.fastq"): 100,
            os.path.join(os.path.dirname(__file__), "fixtures", "reads.fastq.gz"): 100
        }
        for filename in expect_reads:
            self.assertEqual(count_reads(filename), expect_reads[filename])

        with tempfile.NamedTemporaryFile() as tf:
            self.assertEqual(count_reads(tf.name), 0)

            tf.write(b"test")
            tf.flush()
            with self.assertRaises(InvalidInputFileError):
                count_reads(tf.name)

        with tempfile.NamedTemporaryFile() as tf:
            with open(os.path.join(os.path.dirname(__file__), "fixtures", "reads.fastq"), "rb") as fh:
                tf.write(fh.read(90))
                tf.flush()
            with self.assertRaises(InvalidInputFileError):
                count_reads(tf.name)
