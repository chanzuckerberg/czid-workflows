#!/usr/bin/env python3
from typing import Iterator, List, Tuple, NamedTuple
import sys

class Read(NamedTuple):
    header: str
    sequence: str

def iterator(fasta_file: str) -> Iterator[Read]:
    """Iterate through fasta_file, yielding one Read tuple at a time."""
    # TODO: Support full fasta format, where sequences may be split over multiple lines.
    # Perf: 47 million (unpaired) reads per minute on a high end 2018 laptop.
    with open(fasta_file, 'r', encoding='utf-8') as f:
        while True:
            header = f.readline().rstrip()
            sequence = f.readline().rstrip()
            if not header or not sequence:
                break
            # the performance penalty for these asserts is only 8 percent
            assert header[0] == '>'
            assert sequence[0] != '>'
            # the performance penalty for constructing a Read tuple is 40 percent
            yield Read(header, sequence)

def synchronized_iterator(fasta_files: List[str]) -> Iterator[Tuple[Read, ...]]:
    """Iterate through one or more fasta files in lockstep, yielding tuples of
    matching reads.  When the given list fasta_files has length 1, yield
    1-element tuples.  This facilitates uniform processing of either
    unpaired or paired-end reads."""
    return zip(*map(iterator, fasta_files))

def count_reads(fasta_files: List[str]) -> int:
    return sum(1 for _ in synchronized_iterator(fasta_files))

if __name__ == "__main__":
    # Count reads.  Run with fasta filenames as args.  Just for testing.
    print(count_reads(sys.argv[1:]))
