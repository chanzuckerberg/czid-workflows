"""
Generate loc db
Generate Loc DB for NT/NR
"""
import dbm
import re
import shelve
import sys
import logging


def generate_loc_db(db_file, loc_db_file, info_db_file):
    # Logic copied from generate_loc_db_work
    #   slightly changed for writing to shelve format
    loc_dict = shelve.Shelf(dbm.ndbm.open(loc_db_file.replace(".db", ""), 'c'))  # type: ignore
    info_dict = shelve.Shelf(dbm.ndbm.open(info_db_file.replace(".db", ""), 'c'))  # type: ignore
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
            if lines % 100000 == 0:
                logging.info(f"{lines/1000000.0}M lines")
            if line[0] == '>':  # header line
                if seq_len > 0 and len(accession_id) > 0:
                    loc_dict[accession_id] = [seq_offset, header_len, seq_len]
                if seq_bp_len > 0 and len(accession_name) > 0:
                    info_dict[accession_id] = [accession_name, seq_bp_len]

                seq_offset = seq_offset + header_len + seq_len
                header_len = len(line)
                seq_len = 0
                seq_bp_len = 0
                accession_name = ""
                # Sometimes multiple accessions will be mapped to a single sequence.
                # In this case, they will be separated by the \x01 char.
                # To get the accession name, just match until the first \x01.
                s = re.match('^>([^ ]*) ([^\x01]*).*', line)
                if s:
                    accession_id = s.group(1)
                    accession_name = s.group(2)
            else:
                seq_len += len(line)
                seq_bp_len += len(line.strip())
        if seq_len > 0 and len(accession_id) > 0:
            loc_dict[accession_id] = [seq_offset, header_len, seq_len]
        if seq_bp_len > 0 and len(accession_name) > 0:
            info_dict[accession_id] = [accession_name, seq_bp_len]
    loc_dict.close()
    info_dict.close()


if __name__ == '__main__':
    db_file = sys.argv[1]
    loc_db_file = sys.argv[2]
    info_db_file = sys.argv[3]

    generate_loc_db(db_file, loc_db_file, info_db_file)
