import argparse
import logging
import os
import random
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Generator, Iterable, List, Set, TypeVar


import marisa_trie
from Bio import SeqIO
from sourmash import MinHash


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


taxid_lock = threading.Lock()
locked_taxids: Set[int] = set()

class TaxidLock:
    def __init__(self, taxid: int):
        self.taxid = taxid

    def __enter__(self):
        while True:
            with taxid_lock:
                if self.taxid not in locked_taxids:
                    locked_taxids.add(self.taxid)
                    return
            time.sleep(random.random() * 0.1 + 0.05)


    def __exit__(self, exc_type, exc_val, exc_tb):
        with taxid_lock:
            locked_taxids.remove(self.taxid)



T = TypeVar('T')
def chunked(iterable: Iterable[T], n: int) -> Generator[List[T], None, None]:
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == n:
            yield chunk
            chunk = []


def compress_nt(
    nt_filepath: str,
    accession2taxid_path: str,
    k: int,
    scaled: int,
    similarity_threshold: float,
    taxids_to_drop: Set[int],
    parallelism: int,
    compressed_nt_filepath: str,
):
    file_lock = threading.Lock()

    unique_accession_count = 0
    hashes_by_taxid: Dict[int, MinHash] = {}
    accession2taxid = marisa_trie.RecordTrie("L").load(accession2taxid_path)

    with open(compressed_nt_filepath, "w") as f:
        writer = SeqIO.FastaIO.FastaWriter(f)

        def process_record(i: int, record: SeqIO.SeqRecord) -> int:
            # The heavy operations here are executed by sourmash in rust, so we don't need to worry about the GIL
            nonlocal unique_accession_count

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

            with TaxidLock(taxid):
                if taxid not in hashes_by_taxid:
                    hashes_by_taxid[taxid] = min_hash
                    with file_lock:
                        unique_accession_count += 1
                        writer.write_record(record)
                    return i

                total_hash = hashes_by_taxid[taxid]
                # potentially heavy operation
                if min_hash.contained_by(total_hash) > similarity_threshold:
                    return i

                # potentially heavy operation
                total_hash.merge(min_hash)

            with file_lock:
                unique_accession_count += 1
                writer.write_record(record)

            return i

        accession_count = 0
        with ThreadPoolExecutor(max_workers=parallelism) as executor:
            chunks = chunked(enumerate(SeqIO.parse(nt_filepath, 'fasta')), parallelism)
            future_chunks = ({ executor.submit(process_record, i, item) for i, item in chunk} for chunk in chunks)
            # break down into chunks of 1000 to avoid memory issues using itertools chunking
            for futures in future_chunks:
                for future in as_completed(futures):
                    accession_count = max(accession_count, future.result())

    logger.info(f"compressed {accession_count} down to {unique_accession_count} unique accessions")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nt-filepath")
    parser.add_argument("--accession2taxid-path")
    parser.add_argument("--k", type=int)
    parser.add_argument("--scaled", type=int)
    parser.add_argument("--similarity-threshold", type=float)
    parser.add_argument("--taxids-to-drop", nargs="+", type=int, default=[])
    parser.add_argument("--parallelism", type=int, default=os.cpu_count())
    parser.add_argument("--compressed-nt-filepath")
    args = parser.parse_args()

    compress_nt(
        nt_filepath=args.nt_filepath,
        accession2taxid_path=args.accession2taxid_path,
        k=args.k,
        scaled=args.scaled,
        similarity_threshold=args.similarity_threshold,
        taxids_to_drop=set(args.taxids_to_drop),
        parallelism=args.parallelism,
        compressed_nt_filepath=args.compressed_nt_filepath,
    )
