import csv
import gzip
import sys


# We label as 'phage' all of the prokaryotic (bacterial and archaeal) virus families
# listed here: https://en.wikipedia.org/wiki/Bacteriophage
PHAGE_FAMILIES_NAMES = {'Myoviridae', 'Siphoviridae', 'Podoviridae', 'Lipothrixviridae',
                        'Rudiviridae', 'Ampullaviridae', 'Bicaudaviridae', 'Clavaviridae',
                        'Corticoviridae', 'Cystoviridae', 'Fuselloviridae', 'Globuloviridae',
                        'Guttaviridae', 'Inoviridae', 'Leviviridae', 'Microviridae',
                        'Plasmaviridae', 'Tectiviridae'}

phages = {}


def generage_phage_list(versioned_lineages_csv, output_filename):
    with gzip.open(versioned_lineages_csv, "rt") as f:
        for row in csv.DictReader(f):
            taxid = row["taxid"]
            if row["family_name"] in PHAGE_FAMILIES_NAMES:
                entry = phages.get(taxid)
                if not entry:
                    phages[taxid] = {
                        "version_start": row["version_start"],
                        "version_end": row["version_end"],
                    }
                else:
                    version_start = entry["version_start"]
                    if not version_start and row["version_start"]:
                        version_start = row["version_start"]
                    elif version_start:
                        version_start = min(version_start, row["version_start"])

                    version_end = entry["version_end"]
                    if not version_end and row["version_end"]:
                        version_end = row["version_end"]
                    elif version_end:
                        version_end = min(version_end, row["version_end"])

                    phages[taxid] = {
                        "version_start": version_start,
                        "version_end": version_end,
                    }

    with gzip.open(output_filename, "wt") as f:
        writer = csv.DictWriter(f, ["version_start", "version_end", "taxid"])
        for taxid, entry in phages.items():
            writer.writerow(dict(taxid=taxid, **entry))


if __name__ == "__main__":
    generage_phage_list(sys.argv[1], sys.argv[2])
