"""
Accession2Taxid
After alignment, IDseq uses the NCBI accession2taxid database to map accessions to taxonomic IDs.
The full database contains billions of entries, but only ~15% of those are found in either NR or NT databases.
Therefore, the full NCBI accession2taxid database is subsetted to include only the relevant entries and then used
to map accessions to taxonomic IDs.
1. Download NT/NR
2. Extract accessions
3. Subsetting accessions to the ones appearing in NT/NR
"""
import argparse
import logging
from typing import List
import subprocess

import marisa_trie

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


def accession_names(*source_files: str):
    """
    Extracts the accession names from NR/NT
    Source files must be in fasta format where the first word of the header is the accession name
    """
    accession_count = 0
    for source_file in source_files:
        # use grep and cut to extract the accession names from the fasta headers, it is much faster than using Biopython
        proc = subprocess.Popen(f"grep '^>' {source_file} | cut -d' ' -f1", shell=True, stdout=subprocess.PIPE)
        while proc.stdout and (line := proc.stdout.readline()):
            accession_count += 1
            if accession_count % 1_000_000 == 0:
                logger.info(f"\t\textracted {accession_count // 1_000_000}M accession names")
            yield line[1:-1].split(b".", 1)[0]
    logger.info(f"extracted {accession_count} accession names")


def accession_id_to_taxid(mapping_files: List[str], accessions_trie):
    accession_count = 0
    for mapping_file in mapping_files:
        with open(mapping_file, 'rb') as f:
            logger.info(f"\textracting mappings from {mapping_file}")
            for line in f:
                fields = line.split(b'\t')
                accession = fields[0]
                accession_no_version = accession.split(b'.', 1)[0]

                # Only output mappings if the accession is in the source files

                # If using the prot.accession2taxid.FULL file
                if len(fields) < 3 and accession_no_version in accessions_trie:
                    # Remove the version number and the taxid will be at index 1
                    accession_count += 1
                    yield (accession_no_version.decode(), (int(fields[1]),))
                elif accession in accessions_trie:
                    # Otherwise there is a versionless accession ID at index 0 and the taxid is at index 2
                    accession_count += 1
                    yield (accession.decode(), (int(fields[2]),))

                if accession_count % 1_000_000 == 0 and accession_count > 0:
                    logger.info(f"\t\t{accession_count // 1_000_000}M accessions mapped")
    logger.info(f"mapped {accession_count} accessions")


def generate_accession2taxid(accession_mapping_files: List[str], source_files: List[str], accession2taxid_db: str):
    logger.info("starting generating accession2taxid")
    logger.info("extracting accession names from source files")
    # build a marisa-trie of accession names, this is a highly efficient data structure that
    #   enables fast lookups and takes up far less memory than a python set.
    # this is used to filter the accession2taxid database to only include accessions that are in the source files
    accessions_trie = marisa_trie.BinaryTrie(accession_names(*source_files))
    logger.info("extracting mappings from accession mapping files")
    marisa_trie.RecordTrie(
        "L",
        accession_id_to_taxid(accession_mapping_files, accessions_trie),
    ).save(accession2taxid_db)
    logger.info("finished generating accession2taxid")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate accession2taxid.')
    parser.add_argument('accession_mapping_files', nargs=4)
    parser.add_argument('--nt_file')
    parser.add_argument('--nr_file')
    parser.add_argument('--accession2taxid_db')
    args = parser.parse_args()

    source_files = [f for f in [args.nt_file, args.nr_file] if f is not None]

    generate_accession2taxid(
        accession_mapping_files=args.accession_mapping_files,
        source_files=source_files,
        accession2taxid_db=args.accession2taxid_db,
    )
