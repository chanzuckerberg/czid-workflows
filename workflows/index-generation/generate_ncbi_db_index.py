"""
Generate loc db
Generate Loc DB for NT/NR
"""
import sys
import logging
from collections import namedtuple

import marisa_trie


MILLION = 1_000_000


def _accession_id_and_name(header_line):
    # Sometimes multiple accessions will be mapped to a single sequence.
    # In this case, they will be separated by the \x01 char.
    # To get the accession name, just use the string until the first \x01.
    parts = header_line.split('\x01', 1)[0].split(' ', 1)
    accession_id = parts[0][1:]  # Remove the '>'
    accession_name = parts[1] if len(parts) > 1 else ''
    return accession_id, accession_name


Accession = namedtuple('Accession', [
    'accession_id',
    'accession_name',
    'seq_offset',
    'header_len',
    'seq_len',
    'seq_bp_len',
])


def _extract_accessions(db_file: str):
    with open(db_file) as dbf:
        seq_offset = 0
        seq_len = 0
        seq_bp_len = 0
        header_len = 0
        lines = 0
        accession_id = ""
        accession_name = ""
        for line in dbf:
            lines += 1
            if lines % MILLION == 0:
                logging.info(f"{lines/MILLION}M lines")
            if line[0] == '>':  # header line
                if seq_len > 0 and len(accession_id) > 0:
                    yield Accession(accession_id, accession_name, seq_offset, header_len, seq_len, seq_bp_len)

                seq_offset = seq_offset + header_len + seq_len
                header_len = len(line)
                seq_len = 0
                seq_bp_len = 0
                accession_id, accession_name = _accession_id_and_name(line)
            else:
                seq_len += len(line)
                seq_bp_len += len(line.strip())
        if seq_len > 0 and len(accession_id) > 0:
            yield Accession(accession_id, accession_name, seq_offset, header_len, seq_len, seq_bp_len)


def generate_loc_db(db_file, loc_db_file):
    loc_dict = marisa_trie.RecordTrie(
        # unsigned long (8 bytes) followed by 2 unsigned integers (4 bytes)
        # https://docs.python.org/3/library/struct.html#format-characters
        "QII",
        ((a.accession_id, (a.seq_offset, a.header_len, a.seq_len)) for a in _extract_accessions(db_file)),
    )
    loc_dict.save(loc_db_file)


def generate_info_db(db_file, info_db_file):
    info_dict = marisa_trie.RecordTrie(
        # 256 character string followed by an unsigned integer (4 bytes)
        # https://docs.python.org/3/library/struct.html#format-characters
        "256pI",
        ((a.accession_id, (a.accession_name.rstrip().encode(), a.seq_bp_len)) for a in _extract_accessions(db_file)),
    )
    info_dict.save(info_db_file)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    command = sys.argv[1]
    db_file = sys.argv[2]
    index_file = sys.argv[3]

    if command == 'loc':
        generate_loc_db(db_file, index_file)
    elif command == 'info':
        generate_info_db(db_file, index_file)
    else:
        print(f"Unknown command {command}, expected " + " or ".join(['loc', 'info']))
        exit(1)
