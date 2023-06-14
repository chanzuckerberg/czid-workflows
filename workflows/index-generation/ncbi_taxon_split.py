import argparse
import logging
import os
from typing import DefaultDict, Set
from collections import defaultdict


import marisa_trie


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


class TaxidWriteBuffer:
    def __init__(self, output_dir: str, max_size: int) -> None:
        self.output_dir = output_dir
        self.max_size = max_size
        self.bytes_by_taxid: DefaultDict[int, bytes] = defaultdict(bytes)
        self.size = 0

    def write(self, taxid: int, data: bytes) -> None:
        if self.size + len(data) > self.max_size:
            self.flush()
        self.bytes_by_taxid[taxid] += data
        self.size += len(data)

    def flush(self) -> None:
        for taxid, data in self.bytes_by_taxid.items():
            with open(f"{self.output_dir}/{taxid}.fasta", "ab") as f:
                f.write(data)
        self.bytes_by_taxid = defaultdict(bytes)
        self.size = 0


def ncbi_taxon_split(
    nt_filepath: str,
    accession2taxid_path: str,
    taxids_to_drop: Set[int],
    output_dir: str,
    buffer_size: int,
):
    with open(nt_filepath, "rb") as f:
        taxon_count = accession_count = 0
        accession2taxid = marisa_trie.RecordTrie("L").load(accession2taxid_path)
        write_buffer = TaxidWriteBuffer(output_dir, max_size=buffer_size)
        os.makedirs(output_dir, exist_ok=True)
        taxid = None
        for line in f:
            if line.startswith(b">"):
                accession_count += 1
                if accession_count % 1_000_000 == 0:
                    logger.info(f"\t{taxon_count / 1_000_000}M accessions processed")

                non_versioned_accession_id = line.split(b".", 1)[0].decode()
                if non_versioned_accession_id not in accession2taxid:
                    taxid = None
                    continue

                taxid = accession2taxid[non_versioned_accession_id][0][0]
                if taxid in taxids_to_drop:
                    continue

            if taxid:
                write_buffer.write(taxid, line)
        write_buffer.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nt-filepath")
    parser.add_argument("--accession2taxid-path")
    parser.add_argument("--taxids-to-drop", nargs="+", type=int, default=[])
    parser.add_argument("--output-dir")
    parser.add_argument("--buffer-size", type=int, default=1_000_000_000)
    args = parser.parse_args()

    ncbi_taxon_split(
        nt_filepath=args.nt_filepath,
        accession2taxid_path=args.accession2taxid_path,
        taxids_to_drop=set(args.taxids_to_drop),
        output_dir=args.output_dir,
        buffer_size=args.buffer_size,
    )
