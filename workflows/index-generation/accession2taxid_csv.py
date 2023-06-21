import argparse
import csv
from typing import Set

import marisa_trie


def marisa_to_csv(accession2taxid_path: str, taxids_to_drop: Set[int], csv_path: str):
    with open(csv_path, "w") as f:
        trie = marisa_trie.RecordTrie("L").load(accession2taxid_path)
        writer = csv.writer(f)
        writer.writerows([accession, taxid] for accession, (taxid,) in trie.items() if taxid not in taxids_to_drop)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession2taxid-path")
    parser.add_argument("--taxids-to-drop", nargs="+", type=int, default=[])
    parser.add_argument("--csv-filepath")
    args = parser.parse_args()

    marisa_to_csv(
        accession2taxid_path=args.accession2taxid_path,
        taxids_to_drop=set(args.taxids_to_drop),
        csv_path=args.csv_filepath,
    )
