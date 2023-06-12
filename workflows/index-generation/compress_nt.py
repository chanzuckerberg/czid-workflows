import argparse
import logging
import os
from typing import Dict, Set


import marisa_trie
from Bio import SeqIO
from sourmash import MinHash


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


def compress_nt(
    nt_filepath: str,
    accession2taxid_path: str,
    k: int,
    scaled: int,
    similarity_threshold: float,
    taxids_to_drop: Set[int],
    compressed_nt_filepath: str,
):

    unique_accession_count = accession_count = 0
    hashes_by_taxid: Dict[int, MinHash] = {}
    accession2taxid = marisa_trie.RecordTrie("L").load(accession2taxid_path)

    with open(compressed_nt_filepath, "w") as f:
        writer = SeqIO.FastaIO.FastaWriter(f)

        for i, record in enumerate(SeqIO.parse(nt_filepath, 'fasta')):
            accession_count = i
            non_versioned_accession_id = record.id.split(".", 1)[0]
            if non_versioned_accession_id not in accession2taxid:
                return i
            taxid = accession2taxid[non_versioned_accession_id][0][0]
            min_hash = MinHash(n=0, ksize=k, scaled=scaled)
            # NR/NT have some invalid characters in their sequences, force treats k-mers with invalid characters as 0
            min_hash.add_sequence(str(record.seq), force=True)

            if taxid in taxids_to_drop:
                return i

            if i % 100_000 == 0:
                logger.info(f"\t{i / 1_000_000}M accessions processed ({unique_accession_count} unique)")

            if taxid not in hashes_by_taxid:
                hashes_by_taxid[taxid] = min_hash
                unique_accession_count += 1
                writer.write_record(record)
                continue

            total_hash = hashes_by_taxid[taxid]
            if min_hash.contained_by(total_hash) > similarity_threshold:
                continue

            total_hash.merge(min_hash)
            unique_accession_count += 1
            writer.write_record(record)

    logger.info(f"compressed {accession_count} down to {unique_accession_count} unique accessions")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nt-filepath")
    parser.add_argument("--accession2taxid-path")
    parser.add_argument("--k", type=int)
    parser.add_argument("--scaled", type=int)
    parser.add_argument("--similarity-threshold", type=float)
    parser.add_argument("--taxids-to-drop", nargs="+", type=int, default=[])
    parser.add_argument("--compressed-nt-filepath")
    args = parser.parse_args()

    compress_nt(
        nt_filepath=args.nt_filepath,
        accession2taxid_path=args.accession2taxid_path,
        k=args.k,
        scaled=args.scaled,
        similarity_threshold=args.similarity_threshold,
        taxids_to_drop=set(args.taxids_to_drop),
        compressed_nt_filepath=args.compressed_nt_filepath,
    )
