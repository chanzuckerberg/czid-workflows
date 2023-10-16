
import gzip
import csv
import tempfile
import unittest
from datetime import datetime
from os.path import dirname, join

import importlib

lineage_csvs = importlib.import_module("workflows.index-generation.generate_lineage_csvs")

test_files_dir = join(dirname(__file__), "test_files")

class TestIndexGeneration(unittest.TestCase):
    def test_version_taxon_lineages(self):
        with tempfile.NamedTemporaryFile('wb', suffix="csv.gz") as previous_lineages, \
            tempfile.NamedTemporaryFile('wb', suffix="csv.gz") as lineages, \
            tempfile.NamedTemporaryFile('r', suffix="csv.gz") as output:

            with open(join(test_files_dir, "lineages.csv"), "rb") as lineages_raw:
                lineages.write(gzip.compress(lineages_raw.read()))
                lineages.seek(0)

            with open(join(test_files_dir, "previous.csv"), "rb") as previous_lineages_raw:
                previous_lineages.write(gzip.compress(previous_lineages_raw.read()))
                previous_lineages.seek(0)

            version = datetime.now().strftime("%Y-%m-%d")
            lineage_csvs.version_taxon_lineages(
                previous_lineages_filename=previous_lineages.name,
                lineages_filename=lineages.name,
                version=version,
                output_filename=output.name,
            )

            output_unzipped = gzip.open(output.name, "rt")
            result = list(csv.DictReader(output_unzipped))
            self.assertEqual(len(result), 9, result)
