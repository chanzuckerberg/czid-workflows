import os
import unittest
import shutil
from os.path import dirname, join

from idseq_dag.util.fasta import sort_fastx_by_entry_id, multilinefa2singlelinefa

class TestFasta(unittest.TestCase):
    def test_sort_fastx_by_entry_id_fastq(self):
        try:
            fastq = join(dirname(__file__), "..", "fixtures", "reads.fastq")
            tmp_fastq_path = fastq + ".tmp-sort-fastx-by-entry-id-fastq"
            shutil.copyfile(fastq, tmp_fastq_path)
            sort_fastx_by_entry_id(tmp_fastq_path)
            prev = None
            line_count = 0
            with open(tmp_fastq_path) as result:
                for i, line in enumerate(result):
                    if i % 4 == 0:
                        if prev:
                            assert line > prev, f"Expected '{line}' to come before '{prev}'"
                        prev = line
                    line_count = i + 1
            assert line_count % 4 == 0, f"Expected number of lines in fastq to be a factor of 4 but it was {line_count}"
        finally:
            os.remove(tmp_fastq_path)

    def test_sort_fastx_by_entry_id_fasta(self):
        try:
            fasta = join(dirname(__file__), "..", "fixtures", "reads.fasta")
            tmp_fasta_path = fasta + ".tmp-sort-fastx-by-entry-id-fasta"
            # The sort function only supports single line fasta so we must
            #   convert the fixture
            multilinefa2singlelinefa(fasta, tmp_fasta_path)
            sort_fastx_by_entry_id(tmp_fasta_path)
            prev = None
            with open(tmp_fasta_path) as result:
                for i, line in enumerate(result):
                    if i % 4 == 0:
                        if prev:
                            assert line > prev, f"Expected '{line}' to come before '{prev}'"
                        prev = line
                    line_count = i + 1
            assert line_count % 2 == 0, f"Expected number of lines in fasta to be a factor of 2 but it was {line_count}"
        finally:
            os.remove(tmp_fasta_path)
