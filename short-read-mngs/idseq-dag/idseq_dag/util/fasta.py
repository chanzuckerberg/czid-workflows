#!/usr/bin/env python3
from typing import Iterator, List, Tuple, NamedTuple
import sys
import idseq_dag.util.command as command

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

def input_file_type(input_file):
    ''' Check input file type based on first line of file. file needs to be uncompressed '''
    with open(input_file, 'r') as f:
        first_line = f.readline()
    if first_line[0] == '@':
        return 'fastq'
    elif first_line[0] == '>':
        return 'fasta'
    return

def fq2fa(input_fastq, output_fasta):
    ''' FASTQ to FASTA conversion '''
    cmd = f"sed -n '1~4s/^@/>/p;2~4p' <{input_fastq} >{output_fasta}"
    command.execute(cmd)

def multilinefa2singlelinefa(input_fasta, output_fasta):
    ''' Multi-line FASTA to Single-line FASTA conversion '''
    cmd = f"awk 'NR==1 {{print $0}} NR>1 && /^>/ {{printf(\"\\n%s\\n\",$0);next; }} NR>1 {{ printf(\"%s\",$0);}}  END {{printf(\"\\n\");}}' <{input_fasta} > {output_fasta}"
    command.execute(cmd)

if __name__ == "__main__":
    # Count reads.  Run with fasta filenames as args.  Just for testing.
    print(count_reads(sys.argv[1:]))
