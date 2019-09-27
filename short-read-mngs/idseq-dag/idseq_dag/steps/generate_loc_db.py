''' Generate loc db  '''
import dbm
import idseq_dag.util.log as log
import re
import shelve

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.command import run_in_subprocess
from idseq_dag.util.dict import IdSeqDictForUpdate, IdSeqDictValue
BATCH_INSERT_SIZE = 300


class PipelineStepGenerateLocDB(PipelineStep):
    ''' Generate Loc DB for NT/NR '''

    def run(self):
        """
        Generate Loc DB for NT/NR
        """
        db_file = self.input_files_local[0][0]

        loc_db_file_sqlite = self.output_files_local()[0]
        info_db_file_sqlite = self.output_files_local()[1]
        self.generate_loc_db_for_sqlite(db_file, loc_db_file_sqlite, info_db_file_sqlite)

        loc_db_file = self.output_files_local()[2]
        # TODO: (gdingle): this needs to be added
        # info_db_file = self.output_files_local()[3]
        self.generate_loc_db_for_shelf(db_file, loc_db_file)

    @staticmethod
    def generate_loc_db_work(db_file, loc_db, info_db):
        loc_batch_list = []
        info_batch_list = []
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
                    log.write(f"{lines/1000000.0}M lines")
                if line[0] == '>':  # header line
                    if seq_len > 0 and len(accession_id) > 0:
                        loc_batch_list.append((accession_id, [seq_offset, header_len, seq_len]))
                        if len(loc_batch_list) >= BATCH_INSERT_SIZE:
                            loc_db.batch_inserts(loc_batch_list)
                            loc_batch_list = []
                    if seq_bp_len > 0 and len(accession_name) > 0:
                        info_batch_list.append((accession_id, [accession_name, seq_bp_len]))
                        if len(info_batch_list) >= BATCH_INSERT_SIZE:
                            info_db.batch_inserts(info_batch_list)
                            info_batch_list = []

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
                loc_batch_list.append((accession_id, [seq_offset, header_len, seq_len]))
            if seq_bp_len > 0 and len(accession_name) > 0:
                info_batch_list.append((accession_id, [accession_name, seq_bp_len]))
            loc_db.batch_inserts(loc_batch_list)
            info_db.batch_inserts(info_batch_list)

    @run_in_subprocess
    def generate_loc_db_for_sqlite(self, db_file, loc_db_file, info_db_file):
        with IdSeqDictForUpdate(loc_db_file, IdSeqDictValue.VALUE_TYPE_ARRAY) as loc_db, \
                IdSeqDictForUpdate(info_db_file, IdSeqDictValue.VALUE_TYPE_ARRAY) as info_db:
            PipelineStepGenerateLocDB.generate_loc_db_work(db_file, loc_db, info_db)

    @run_in_subprocess
    def generate_loc_db_for_shelf(self, db_file, loc_db_file):
        loc_dict = shelve.Shelf(dbm.ndbm.open(loc_db_file.replace(".db", ""), 'c'))
        with open(db_file) as dbf:
            seq_offset = 0
            seq_len = 0
            header_len = 0
            lines = 0
            accession_id = ""
            for line in dbf:
                lines += 1
                if lines % 100000 == 0:
                    log.write(f"{lines/1000000.0}M lines")
                if line[0] == '>':  # header line
                    if seq_len > 0 and len(accession_id) > 0:
                        loc_dict[accession_id] = (seq_offset, header_len, seq_len)
                    seq_offset = seq_offset + header_len + seq_len
                    header_len = len(line)
                    seq_len = 0
                    s = re.match('^>([^ ]*).*', line)
                    if s:
                        accession_id = s.group(1)
                else:
                    seq_len += len(line)
            if seq_len > 0 and len(accession_id) > 0:
                loc_dict[accession_id] = (seq_offset, header_len, seq_len)
        loc_dict.close()

    def count_reads(self):
        ''' Count reads '''
        pass
