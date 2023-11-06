import gzip
import csv
import tempfile
import unittest
from datetime import datetime
from os.path import dirname, join

import importlib

lineage_csvs = importlib.import_module(
    "workflows.index-generation.generate_lineage_csvs"
)

test_files_dir = join(dirname(__file__), "test_files")


class TestIndexGeneration(unittest.TestCase):
    unchanged_taxon_id = 2  # in both old and new csv -> 1 row
    unchanged_taxon_id_with_multiple_rows = (
        3  # unchanged in this update, but had past updates -> 2 total rows
    )
    multiple_row_taxid = 4  # has 2 rows in old csv and a new row in current csv -> 3 rows in versioned csv
    taxid_with_updated_genus = 5  # genus changed between old csv and current csv
    taxid_missing_genus = 6  # has -200 for genus, is assigned to new family in new CSV -> 2 rows in versioned csv
    deprecated_taxid = 7  # taxon present in old csv, not present in current csv
    deprecated_taxid_multiple_rows = (
        8  # taxon present with multiple rows in old csv, not present in current csv
    )
    taxid_phage_to_nonphage = 9  # taxon that changes from phage to non-phage
    new_taxid = 10  # taxon not present in old csv, present in new csv

    version = datetime.now().strftime("%Y-%m-%d")

    def _find_entries_for_taxid(self, results, taxid, sort=True):
        entries = [x for x in results if int(x["taxid"]) == taxid]
        if sort:
            entries.sort(key=lambda x: x["version_end"])
        return entries

    def _unchanged_taxa_assertions(self, results):
        unchanged_taxon_single_row = self._find_entries_for_taxid(
            results, self.unchanged_taxon_id
        )
        self.assertEqual(len(unchanged_taxon_single_row), 1)
        self.assertEqual(unchanged_taxon_single_row[0]["version_end"], self.version)
        unchanged_taxon_multiple_rows = self._find_entries_for_taxid(
            results, self.unchanged_taxon_id_with_multiple_rows
        )
        self.assertEqual(len(unchanged_taxon_multiple_rows), 2)
        # second entry (more recent one) should have its version_end updated
        self.assertNotEqual(
            unchanged_taxon_multiple_rows[0]["version_end"], self.version
        )
        self.assertEqual(unchanged_taxon_multiple_rows[1]["version_end"], self.version)

    def _deprecated_taxa_assertions(self, results):
        deprecated_taxon_rows = self._find_entries_for_taxid(
            results, self.deprecated_taxid
        )
        # deprecated row should still be present in new lineages csv
        self.assertEqual(len(deprecated_taxon_rows), 1)
        # deprecated row should not have its version_end updated
        self.assertNotEqual(deprecated_taxon_rows[0]["version_end"], self.version)

    def _deprecated_taxa_multiple_rows_assertions(self, results):
        deprecated_taxon_multiple_rows = self._find_entries_for_taxid(
            results, self.deprecated_taxid_multiple_rows
        )
        self.assertEqual(len(deprecated_taxon_multiple_rows), 2)
        self.assertNotEqual(
            deprecated_taxon_multiple_rows[0]["version_end"], self.version
        )
        self.assertNotEqual(
            deprecated_taxon_multiple_rows[1]["version_end"], self.version
        )

    def _updated_taxa_assertions(self, results):
        updated_taxa_rows = self._find_entries_for_taxid(
            results, self.taxid_with_updated_genus
        )
        self.assertEqual(len(updated_taxa_rows), 2, updated_taxa_rows)
        # old row should not have its version_end updated
        self.assertNotEqual(updated_taxa_rows[0]["version_end"], self.version)
        # new row should have start and end as the new version
        self.assertEqual(updated_taxa_rows[1]["version_start"], self.version)
        self.assertEqual(updated_taxa_rows[1]["version_end"], self.version)

    def _multiple_taxa_rows_assertions(self, results):
        taxon_with_multiple_rows = self._find_entries_for_taxid(
            results, self.multiple_row_taxid
        )
        self.assertEqual(len(taxon_with_multiple_rows), 3)
        # old rows should not have their version_end field updated
        self.assertNotEqual(taxon_with_multiple_rows[0]["version_end"], self.version)
        self.assertNotEqual(taxon_with_multiple_rows[1]["version_end"], self.version)
        # last row should have its start and end as new version
        self.assertEqual(taxon_with_multiple_rows[2]["version_start"], self.version)
        self.assertEqual(taxon_with_multiple_rows[2]["version_end"], self.version)

    def _updated_taxa_missing_genus_assertions(self, results):
        taxon_rows_no_genus_updated_family = self._find_entries_for_taxid(
            results, self.taxid_missing_genus
        )
        self.assertEqual(len(taxon_rows_no_genus_updated_family), 2)
        # old row should not have its version_end updated
        old_row = taxon_rows_no_genus_updated_family[0]
        self.assertNotEqual(old_row["version_end"], self.version)
        # new row should have its start and end as new version
        new_row = taxon_rows_no_genus_updated_family[1]
        self.assertEqual(new_row["version_start"], self.version)
        self.assertEqual(new_row["version_end"], self.version)

        self.assertNotEqual(old_row["family_taxid"], new_row["family_taxid"])

    def _new_taxid_assertions(self, results):
        new_taxon_rows = self._find_entries_for_taxid(results, self.new_taxid)
        self.assertEqual(len(new_taxon_rows), 1)
        self.assertEqual(new_taxon_rows[0]["version_start"], self.version)
        self.assertEqual(new_taxon_rows[0]["version_end"], self.version)

    def _phage_to_nonphage_assertions(self, results):
        phage_rows = self._find_entries_for_taxid(results, self.taxid_phage_to_nonphage)
        self.assertEqual(len(phage_rows), 2)
        self.assertNotEqual(phage_rows[0]["version_end"], self.version)
        self.assertEqual(phage_rows[1]["version_start"], self.version)
        self.assertEqual(phage_rows[1]["version_end"], self.version)

    def test_version_taxon_lineages(self):
        with tempfile.NamedTemporaryFile(
            "wb", suffix="csv.gz"
        ) as previous_lineages, tempfile.NamedTemporaryFile(
            "wb", suffix="csv.gz"
        ) as lineages, tempfile.NamedTemporaryFile(
            "r", suffix="csv.gz"
        ) as output:
            with open(join(test_files_dir, "lineages.csv"), "rb") as lineages_raw:
                lineages.write(gzip.compress(lineages_raw.read()))
                lineages.seek(0)

            with open(
                join(test_files_dir, "previous.csv"), "rb"
            ) as previous_lineages_raw:
                previous_lineages.write(gzip.compress(previous_lineages_raw.read()))
                previous_lineages.seek(0)

            lineage_csvs.version_taxon_lineages(
                previous_lineages_filename=previous_lineages.name,
                lineages_filename=lineages.name,
                version=self.version,
                output_filename=output.name,
            )

            output_unzipped = gzip.open(output.name, "rt")
            result = list(csv.DictReader(output_unzipped))

            self._unchanged_taxa_assertions(result)
            self._deprecated_taxa_assertions(result)
            self._deprecated_taxa_multiple_rows_assertions(result)
            self._updated_taxa_assertions(result)
            self._multiple_taxa_rows_assertions(result)
            self._updated_taxa_missing_genus_assertions(result)
            self._new_taxid_assertions(result)
            self._phage_to_nonphage_assertions(result)
            self.assertEqual(len(result), 17, result)

    def _is_phage_assertions(self, results):
        phage_taxid = 373406  # Mycobacterium virus Cooper
        non_phage_taxid = 1031976  # Influenza A
        phage_taxid_without_family = 28368 # Corynebacterium diphtheriae phage

        phage_row = self._find_entries_for_taxid(results, phage_taxid, sort=False)[0]
        self.assertEqual(int(phage_row["is_phage"]), 1)
        non_phage_row = self._find_entries_for_taxid(
            results, non_phage_taxid, sort=False
        )[0]
        self.assertEqual(int(non_phage_row["is_phage"]), 0)

        phage_row_without_family = self._find_entries_for_taxid(results, phage_taxid_without_family, sort=False)[0]
        self.assertEqual(int(phage_row_without_family["is_phage"]), 1)

    def test_generate_taxon_lineage_names(self):
        with tempfile.NamedTemporaryFile(
            "wb", suffix="csv.gz"
        ) as names_file, tempfile.NamedTemporaryFile(
            "wb", suffix="csv.gz"
        ) as lineages_taxid_file, tempfile.NamedTemporaryFile(
            "r", suffix="csv.gz"
        ) as output:
            with open(join(test_files_dir, "names_file_test.csv"), "rb") as names_raw:
                names_file.write(gzip.compress(names_raw.read()))
                names_file.seek(0)

            with open(
                join(test_files_dir, "lineages_file_test.csv"), "rb"
            ) as lineages_raw:
                lineages_taxid_file.write(gzip.compress(lineages_raw.read()))
                lineages_taxid_file.seek(0)

            lineage_csvs.generate_taxon_lineage_names(
                names_filename=names_file.name,
                lineages_filename=lineages_taxid_file.name,
                output_filename=output.name,
            )

            output_unzipped = gzip.open(output.name, "rt")
            result = list(csv.DictReader(output_unzipped))
            self._is_phage_assertions(result)
