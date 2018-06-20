import multiprocessing
import random
import subprocess
import sys
import threading
import time
from functools import wraps
import idseq_dag.util.log as log

# TODO(cdebourcy): use samtools to convert the sam files to fasta instead

def generate_unmapped_singles_from_sam(sam_file, output_fa):
    """Output a single file containing every unmapped read after bowtie2.

    SAM file alignments:
    - See: https://en.wikipedia.org/wiki/SAM_(file_format)
         - https://broadinstitute.github.io/picard/explain-flags.html
    - part[0] = query template name
    - part[1] = bitwise flag
    - part[9] = segment sequence
    """
    with open(output_fa, 'w') as output_read:
        with open(sam_file, 'r', encoding='utf-8') as sam_f:
            # Skip headers
            read = sam_f.readline()
            while read and read[0] == '@':
                read = sam_f.readline()
            while read:
                part = read.split("\t")
                # Read unmapped
                if part[1] == "4":
                    # Do NOT append /1 to read id
                    output_read.write(">%s\n%s\n" % (part[0], part[9]))
                read = sam_f.readline()

def generate_unmapped_pairs_from_sam(sam_file, out_fas):
    """
    Output out_fas[0] and out_fas[1] containing the unmapped pairs from bowtie2.
    Optionl: out_fas[2] multiplex read ids by appending /1 and /2.
    SAM file alignments:
    - See: https://en.wikipedia.org/wiki/SAM_(file_format)
         - https://broadinstitute.github.io/picard/explain-flags.html
    - part[0] = query template name
    - part[1] = bitwise flag
    - part[9] = segment sequence
    """

    assert(len(out_fas) == 2 or len(out_fas) == 3)
    out_fa_1 = open(out_fas[0], 'w')
    out_fa_2 = open(out_fas[1], 'w')
    out_fa_merged = open(out_fas[2], 'w') if len(out_fas) == 3 else None
    sam_f = open(sam_file, 'r', encoding='utf-8')
    # Skip headers
    read1 = sam_f.readline()
    while read1 and read1[0] == '@':
        read1 = sam_f.readline()
    read2 = sam_f.readline()

    while read1 and read2:
        part1 = read1.split("\t")
        part2 = read2.split("\t")
        if part1[1] == "77" and part2[1] == "141":  # Both parts unmapped
            out_fa_1.write(">%s\n%s\n" % (part1[0], part1[9]))
            out_fa_2.write(">%s\n%s\n" % (part2[0], part2[9]))
            if out_fa_merged:
                # Append /1 to read id
                out_fa_merged.write(">%s/1\n%s\n" % (part1[0], part1[9]))
                # Append /2 to read id
                out_fa_merged.write(">%s/2\n%s\n" % (part2[0], part2[9]))
        read1 = sam_f.readline()
        read2 = sam_f.readline()
    out_fa_1.close()
    out_fa_2.close()
    if out_fa_merged:
        out_fa_merged.close()
