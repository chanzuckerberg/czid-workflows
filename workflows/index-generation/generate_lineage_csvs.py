import csv
import gzip
import sys
import logging

from datetime import datetime

from typing import Dict, Union

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s|%(levelname)s|%(message)s")

_taxon_levels = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

_fieldnames = [
    "taxid",
    "tax_name",
    "is_phage",
] + [
    f"{level}_{label}"
    for level in _taxon_levels
    for label in ["taxid", "name", "common_name"]
]

_versioning_fieldnames = [
    "version_start",
    "version_end",
    "created_at",
    "updated_at",
]

PHAGE_FAMILIES_NAMES = {
    "Myoviridae",
    "Siphoviridae",
    "Podoviridae",
    "Lipothrixviridae",
    "Rudiviridae",
    "Ampullaviridae",
    "Bicaudaviridae",
    "Clavaviridae",
    "Corticoviridae",
    "Cystoviridae",
    "Fuselloviridae",
    "Globuloviridae",
    "Guttaviridae",
    "Inoviridae",
    "Leviviridae",
    "Microviridae",
    "Plasmaviridae",
    "Tectiviridae",
}


def generate_taxon_lineage_names(
    names_filename: str,
    lineages_filename: str,
    output_filename: str,
):
    """
    Combines the outputs from https://github.com/chanzuckerberg/ncbitax2lin to create
    a single CSV with the lineage of each taxon, including IDs, names, and common names.
    Inputs:
    - names_filename: name of a gzipped CSV file with the following columns: tax_id,name_txt,name_txt_common
    - lineages_filename: name of a gzipped CSV file with the following columns:
        tax_id,superkingdom,kingdom,phylum,class,order,family,genus,species,no_rank,no_rank1,no_rank2,no_rank3,no_rank4
      In this file, taxonomic levels such as kingdom, or phylum are populated by IDs rather than names
    Outputs a gzipped CSV with the following columns:
    - taxid
    - tax_name
    - is_phage
    - superkingdom_taxid
    - phylum_taxid
    - class_taxid
    - order_taxid
    - family_taxid
    - genus_taxid
    - species_taxid
    - superkingdom_name
    - phylum_name
    - class_name
    - order_name
    - family_name
    - genus_name
    - species_name
    - superkingdom_common_name
    - phylum_common_name
    - class_common_name
    - order_common_name
    - family_common_name
    - genus_common_name
    - species_common_name
    Columns other than tax_id may be blank. This CSV is supposed to correspond to the taxon_lineages table in the
    CZID web application without any versioning or date information.
    """

    names = {}
    with gzip.open(names_filename, "rt") as f:
        for row in csv.DictReader(f):
            names[row["tax_id"]] = (row["name_txt"], row["name_txt_common"])

    with gzip.open(lineages_filename, "rt") as rf, gzip.open(
        output_filename, "wt"
    ) as wf:
        writer = csv.DictWriter(wf, fieldnames=_fieldnames)
        writer.writeheader()

        for row in csv.DictReader(rf):
            tax_name = names.get(row["tax_id"], ("", ""))[0]
            new_row = {
                "taxid": row["tax_id"],
                "tax_name": tax_name,
                "is_phage": 1 if tax_name in PHAGE_FAMILIES_NAMES else 0,
            }

            for level in _taxon_levels:
                new_row[f"{level}_taxid"] = row[level]
                name, common_name = names.get(row[level], ("", ""))
                new_row[f"{level}_name"], new_row[f"{level}_common_name"] = (
                    name,
                    common_name,
                )
            writer.writerow(new_row)


def _equals(previous_row: Dict[str, str], row: Dict[str, str]):
    """
    Checks if fields in row equals a previous row, ignoring versioning info
    previous_row will have versioning information already
    row has no versioning information so this just
    checks for equality in fields that are not related to versioning
    """
    for fieldname in _fieldnames:
        if previous_row[fieldname] != row[fieldname]:
            return False

    return True

def _find_highest_taxonomy_changed(previous_row: Dict[str, str], row: Dict[str, str]):
    fieldnames_to_search = [f'{level}_name' for level in _taxon_levels]
    # searches in order of taxid (should be the same), is_phage, superkingdom -> kingdom -> ... -> species
    for fieldname in fieldnames_to_search:
        if previous_row[fieldname] != row[fieldname]:
            return (fieldname, previous_row[fieldname], row[fieldname])


def version_taxon_lineages(
    previous_lineages_filename: Union[str, None],
    lineages_filename: str,
    version: str,
    output_filename: str,
):
    """
    Combines a non-versioned CSV output from generate_taxon_lineage_names with a previous
    output from this function to create a versioned lineage CSV. A versioned lineage CSV
    has historical versions of lineages and the version ranges for which those are valid.
    For example, let's say a lineage is just genus_taxid and species_taxid, the table may
    look something like:
    taxid,version_start,version_end,is_phage,genus_taxid,species_taxid
      123,   2020-01-01, 2020-02-01,       0,         10,          123
      123,   2020-03-01, 2020-03-01,       1,         20,          123
      456,   2020-01-01, 2020-01-01,       0,         10,          456
      789,   2020-03-01, 2020-03-01,       1,         20,          789
    In this table:
    - The 123 species had it's genus changed between the 2020-02-01 version and
      the 2020-03-01 version, so we see two lineages for it. The old one, valid from 2020-01-01
      to 2020-02-01, and the new one, valid from 2020-03-01 to 2020-03-01.
    - The 456 species has been removed between the 2020-01-01 and 2020-02-01 versions, it is
      still in the table, but there are no entries for it under any versions after 2020-01-01
    - The 789 species has been added in the 2020-03-01 version, there is an entry for it in
      that new version but there are no entries for previous versions
    The output of this function is intended to be loaded raw into the CZID web application's
    taxon_lineages table.
    Inputs:
    - previous_lineages_filename: name of a previous output from this file, a gzipped
      CSV file with the following colums (note the inclusion of version_start, version_end,
      created_at, updated_at):
        - taxid
        - tax_name
        - is_phage
        - superkingdom_taxid
        - phylum_taxid
        - class_taxid
        - order_taxid
        - family_taxid
        - genus_taxid
        - species_taxid
        - superkingdom_name
        - phylum_name
        - class_name
        - order_name
        - family_name
        - genus_name
        - species_name
        - superkingdom_common_name
        - phylum_common_name
        - class_common_name
        - order_common_name
        - family_common_name
        - genus_common_name
        - species_common_name
        - version_start
        - version_end
        - created_at
        - updated_at
    - lineages_filename: name the output from generate_taxon_lineage_names
    - version: a version string for this lineage version, a date string in the form: 2020-04-20
    Outputs a gzipped CSV with the new lineage information merged into the old using the versioning
    scheme described above
    """
    previous_lineages = {}
    previous_lineages_version = None
    if previous_lineages_filename:
        with gzip.open(previous_lineages_filename, "rt") as f:
            for row in csv.DictReader(f):
                previous_lineages[(row["taxid"], row["version_end"])] = row
                if (
                    not previous_lineages_version
                    or previous_lineages_version < row["version_end"]
                ):
                    previous_lineages_version = row["version_end"]

    num_existing_rows = len(previous_lineages)
    logging.info(
        f"Number of rows in existing taxon lineages table: {num_existing_rows}"
    )

    with gzip.open(output_filename, "wt") as wf, open("changelog.csv", 'w') as log, open("deleted_log.csv", 'w') as deleted_log, open("new_taxa", 'w') as new_taxa_log:
        writer = csv.DictWriter(wf, fieldnames=_fieldnames + _versioning_fieldnames)
        writer.writeheader()
        # keep track of counts of different types of taxa
        # this allows us to spot check results without loading output file into memory
        num_unchanged_rows = 0
        num_new_taxa_rows = 0
        num_updated_lineage_rows = 0
        num_total_new_rows = 0
        num_deprecated_rows = 0

        log_writer = csv.writer(log)
        log_writer.writerow(["taxid", "taxonomy_level", "old_value", "new_value"])

        with gzip.open(lineages_filename, "rt") as rf:
            for row in csv.DictReader(rf):
                previous_row = previous_lineages.pop(
                    (row["taxid"], previous_lineages_version), None
                )

                if previous_row and _equals(row, previous_row):
                    # We already have this exact lineage, update its version_end
                    #   to keep it from expiring
                    previous_row["version_end"] = version
                    previous_row["updated_at"] = str(datetime.now())
                    writer.writerow(previous_row)
                    num_unchanged_rows += 1
                    num_total_new_rows += 1
                else:
                    # This is either a brand new lineage, or an updated
                    #   lineage. Create a new lineage, and don't update
                    #   the version_end of the possible old lineage so it will
                    #   expire
                    row["version_start"] = version
                    row["version_end"] = version
                    row["created_at"] = str(datetime.now())
                    row["updated_at"] = str(datetime.now())
                    writer.writerow(row)
                    num_total_new_rows += 1

                    if previous_row:
                        (taxonomy_level, old_val, new_val) = _find_highest_taxonomy_changed(previous_row, row)
                        print(taxonomy_level, old_val, new_val)
                        log_writer.writerow([taxonomy_level, old_val, new_val])
                        # this is an updated lineage
                        writer.writerow(previous_row)
                        num_total_new_rows += 1
                        num_updated_lineage_rows += 1
                        num_deprecated_rows += 1
                    else:
                        # this is a new lineage
                        num_new_taxa_rows += 1

            for previous_row in previous_lineages.values():
                # All rows left in previous_lineages are for taxons that have
                #   been removed or outdated lineages for existing taxa. We still need to
                #   write them to the new output file so we have them for older versions,
                #   they just won't have their version updated so they will be considered expired.
                writer.writerow(previous_row)
                num_deprecated_rows += 1
                num_total_new_rows += 1

        summary_counts = (
            f"Number of taxa with unchanged lineages: {num_unchanged_rows}\n"
            f"Number of taxa with updated lineages: {num_updated_lineage_rows}\n"
            f"Number of new taxa: {num_new_taxa_rows}\n"
            f"Number of deprecated lineage rows (outdated lineage or deprecated taxa): {num_deprecated_rows}\n"
            f"Number of total rows written to new table: {num_total_new_rows}"
        )
        logging.info(summary_counts)

        # Assert that the correct number of rows have been written to the new table
        # and that we've correctly calculated the number of taxa in each category
        expected_existing_num_rows = num_unchanged_rows + num_deprecated_rows
        if not expected_existing_num_rows == num_existing_rows:
            logging.warning(
                """
                Number of expected existing rows (deprecated lineages and unchanged lineages)
                 %s does not match number of rows in taxon lineages table %s
                """
                % (expected_existing_num_rows, num_existing_rows)
            )

        expected_total_new_rows = (
            num_existing_rows + num_updated_lineage_rows + num_new_taxa_rows
        )
        if not expected_total_new_rows == num_total_new_rows:
            logging.warning(
                """
                Expected number of rows in new table (length of old table + updated rows + new rows)
                 %s does not match number of rows written %s
                """
                % (expected_total_new_rows, num_total_new_rows)
            )


if __name__ == "__main__":
    names_filename = sys.argv[1]
    lineages_filename = sys.argv[2]
    version = sys.argv[3]
    raw_lineage_output_filename = sys.argv[4]
    versioned_lineage_output_filename = sys.argv[5]
    previous_lineages_filename = sys.argv[6] if len(sys.argv) > 6 else None

    generate_taxon_lineage_names(
        names_filename,
        lineages_filename,
        raw_lineage_output_filename,
    )

    version_taxon_lineages(
        previous_lineages_filename,
        raw_lineage_output_filename,
        version,
        versioned_lineage_output_filename,
    )
