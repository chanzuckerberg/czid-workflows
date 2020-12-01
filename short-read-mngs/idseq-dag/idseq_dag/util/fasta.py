#!/usr/bin/env python3
from typing import Iterator, List, Tuple, NamedTuple
import sys
import os
from subprocess import run
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

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

def _count_reads(fasta_files: List[str]) -> int:
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
    command.execute(
        command_patterns.ShellScriptCommand(
            script=r'''sed -n '1~4s/^@/>/p;2~4p' <"${input_fastq}" > "${output_fasta}";''',
            named_args={
                'input_fastq': input_fastq,
                'output_fasta': output_fasta
            }
        )
    )


def multilinefa2singlelinefa(input_fasta, output_fasta):
    ''' Multi-line FASTA to Single-line FASTA conversion '''
    command.execute(
        command_patterns.ShellScriptCommand(
            script=r'''awk 'NR==1 {print $0} NR>1 && /^>/ {printf("\n%s\n",$0);next; } NR>1 { printf("%s",$0);}  END {printf("\n");}' <"${input_fasta}" > "${output_fasta}";''',
            named_args={
                'input_fasta': input_fasta,
                'output_fasta': output_fasta
            }
        )
    )


def sort_fastx_by_entry_id(fastq_path):
    tmp_sorted_path = fastq_path + ".sorted"
    with open(fastq_path, 'rb') as in_file:
        with open(tmp_sorted_path, 'wb') as out_file:
            # Command based on this https://www.biostars.org/p/15011/#103041
            if input_file_type(fastq_path) == 'fastq':
                cmd = "paste -d '`' - - - - | sort -k1,1 -S 3G | tr '`' '\n'"
            else:
                # WARNING: does not support multiline fasta
                cmd = "paste -d '`' - - | sort -k1,1 -S 3G | tr '`' '\n'"
            # By default the sort utility uses a locale-based sort, this is significantly
            #   slower than a simple byte comparison. It also produces a different
            #   order than python's default string comparisons would which makes testing
            #   a bit less convenient. All we care about is producing a consistent order
            #   every time, the order itself is irrelevant, so we set LC_ALL=C to do a
            #   simple byte comparison instead of a locale-based sort which is faster,
            #   produces a consistent result regardless of locale, and produces the same
            #   order python's default string comparison would.
            run(cmd, env={'LC_ALL': 'C'}, stdin=in_file, stdout=out_file, check=True, shell=True)
    os.rename(tmp_sorted_path, fastq_path)


if __name__ == "__main__":
    # Count reads.  Run with fasta filenames as args.  Just for testing.
    print(_count_reads(sys.argv[1:]))
