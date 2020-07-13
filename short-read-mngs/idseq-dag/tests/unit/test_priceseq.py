import os
import sys
import tempfile
import unittest

from idseq_dag.steps.run_priceseq import PipelineStepRunPriceSeq
from idseq_dag.exceptions import InvalidInputFileError

fasta_filename = os.path.join(os.path.dirname(__file__), "fixtures", "reads.fasta")
fastq_filename = os.path.join(os.path.dirname(__file__), "fixtures", "reads.fastq")

@unittest.skipIf(os.uname().sysname != "Linux", "Skipping test on incompatible platform")
class TestPriceSeq(unittest.TestCase):
    def run_priceseqfilter(self, **args):
        with tempfile.TemporaryDirectory() as td, tempfile.NamedTemporaryFile(suffix="." + args["file_type"]) as tf:
            args["out_files"] = [tf.name] * len(args["in_files"])
            try:
                wd = os.getcwd()
                os.chdir(td)
                PipelineStepRunPriceSeq.run_priceseqfilter(None, **args)
            finally:
                os.chdir(wd)

    def test_run_priceseqfilter(self):
        for args in [dict(in_files=[fasta_filename], is_paired=False, file_type="fasta"),
                     dict(in_files=[fastq_filename], is_paired=False, file_type="fastq"),
                     dict(in_files=[fasta_filename, fasta_filename], is_paired=True, file_type="fasta"),
                     dict(in_files=[fastq_filename, fastq_filename], is_paired=True, file_type="fastq")]:
            self.run_priceseqfilter(**args)
        for args in [dict(in_files=[fasta_filename], is_paired=True, file_type="fastq"),
                     dict(in_files=[fastq_filename], is_paired=True, file_type="fasta"),
                     dict(in_files=[fasta_filename, fasta_filename], is_paired=False, file_type="fasta"),
                     dict(in_files=[fastq_filename, fastq_filename], is_paired=False, file_type="fastq")]:
            with self.assertRaises(Exception):
                self.run_priceseqfilter(**args)
