import os
import threading
import time
import re
import traceback
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.s3 as s3
import idseq_dag.util.m8 as m8

from idseq_dag.util.dict import IdSeqDictValue, open_file_db_by_extension
from idseq_dag.util.trace_lock import TraceLock

MIN_ACCESSIONS_WHOLE_DB_DOWNLOAD = 5000
MAX_ACCESSION_SEQUENCE_LEN = 100000000
ALLOW_S3MI = False  # Allow s3mi only if running on an instance with enough RAM to fit NT and NR together...

class PipelineStepDownloadAccessions(PipelineStep):
    '''
        Download Accessions Based on the hit summary
    '''
    def run(self):
        (_align_m8, _deduped_m8, hit_summary, _orig_counts) = self.input_files_local[0]
        output_reference_fasta = self.output_files_local()[0]
        loc_db = s3.fetch_reference(
            self.additional_files["loc_db"],
            self.ref_dir_local,
            auto_unzip=True,  # This is default for references, but let's be explicit.
            allow_s3mi=ALLOW_S3MI)
        db_s3_path = self.additional_attributes["db"]
        db_type = self.additional_attributes["db_type"]
        (_read_dict, accession_dict, _selected_genera) = m8.summarize_hits(hit_summary)
        with open_file_db_by_extension(loc_db, IdSeqDictValue.VALUE_TYPE_ARRAY) as loc_dict:
            if len(accession_dict) < MIN_ACCESSIONS_WHOLE_DB_DOWNLOAD:
                self.download_ref_sequences_from_s3(accession_dict, output_reference_fasta, db_type,
                                                    loc_dict, db_s3_path)
            else:
                # download the whole alignment db
                #
                db_path = s3.fetch_reference(
                    db_s3_path,
                    self.ref_dir_local,
                    auto_unzip=True,   # This is default for references, but let's be explicit
                    allow_s3mi=ALLOW_S3MI)
                self.download_ref_sequences_from_file(accession_dict, loc_dict, db_path, output_reference_fasta)

    FIX_COMMA_REGEXP = re.compile(r'^(?P<accession_id>[^ ]+) (?P<wrong_pattern>, *)?(?P<description>.*)$')

    @staticmethod
    def _fix_ncbi_record(accession_data: str) -> str:
        '''
            We found a ncbi record that is oddly annotated (title starts with a comma,
            as you can see here: https://www.ncbi.nlm.nih.gov/protein/XP_002289390.1)
            That produce wrong results in blastx (blastx understand this accession
            as being "XP_002289390.1," instead of "XP_002289390.1")

            This method detects and removes this trailling comma from the title, and
            it is being invoked before writing it to file assembly/nr.refseq.fasta,
            which is a subset of the original nr index.

            `accession_data` is a string that looks like this:

                accession_data = (
                    '>XP_002289390.1 , partial [Thalassiosira pseudonana CCMP1335]\x01EED92927.1 Conserved Hypothetical Protein, partial [Thalassiosira pseudonana CCMP1335]\n'
                    'GGREKKKLLKSQKDGSAKDRHNPRAFSVANIVRTQRNVQRNLDRAQKKEYVPLSDRRAARVEEGPPSLVAVVGPPGVGKS\n'
                    'TLIRSLVKLYTNHNLTNPTGPITVCTSQTKRITFLECPNTPTAMLDVAKIADLVLLCVDAKFGFEMETFEFLNMMQTHGF\n'
                    'PKVMGIFTHLDQFRTQKNLRKTKKLLKHRFWTEIYDGAKMFYFSGCVNGKYLKHEVKQLSLLLSRIKYRPLVWRNTHPYV\n'
                    # (...more lines with sequence data...)
                )
                # note entry above represents a single string and line breaks and ^A character are part of it.

            this method returns same string with title correction:

            result = (
                '>XP_002289390.1 partial [Thalassiosira pseudonana CCMP1335]\x01EED92927.1 Conserved Hypothetical Protein, partial [Thalassiosira pseudonana CCMP1335]\n'
                'GGREKKKLLKSQKDGSAKDRHNPRAFSVANIVRTQRNVQRNLDRAQKKEYVPLSDRRAARVEEGPPSLVAVVGPPGVGKS\n'
                'TLIRSLVKLYTNHNLTNPTGPITVCTSQTKRITFLECPNTPTAMLDVAKIADLVLLCVDAKFGFEMETFEFLNMMQTHGF\n'
                'PKVMGIFTHLDQFRTQKNLRKTKKLLKHRFWTEIYDGAKMFYFSGCVNGKYLKHEVKQLSLLLSRIKYRPLVWRNTHPYV\n'
                # (...more lines with sequence data...)
            )

            Right now this method is only being used to remove the heading comma,
            since this is the only case that we found so far affecting the pipeline.
            It may be extended the future to handle more exceptions if it is needed.
        '''
        def _fix_headers(line):
            if len(line) > 0 and line[0] == ">":
                # support for multiheader line separted by CTRL_A (https://en.wikipedia.org/wiki/FASTA_format#Description_line)
                header_items = line.lstrip(">").split("\x01")
                header_items = (PipelineStepDownloadAccessions.FIX_COMMA_REGEXP.sub(r"\g<accession_id> \g<description>", header_item) for header_item in header_items)
                return ">" + ("\x01".join(header_items))
            return line

        lines = accession_data.split("\n")
        lines = (_fix_headers(line) for line in lines)
        return "\n".join(lines)

    def download_ref_sequences_from_file(self, accession_dict, loc_dict, db_path,
                                         output_reference_fasta):
        db_file = open(db_path, 'r')
        with open(output_reference_fasta, 'w') as orf:
            for accession, _taxinfo in accession_dict.items():
                accession_data = self.get_sequence_by_accession_from_file(accession, loc_dict, db_file)
                if accession_data:
                    accession_data = self._fix_ncbi_record(accession_data)
                    orf.write(accession_data)
        db_file.close()

    def download_ref_sequences_from_s3(self, accession_dict, output_reference_fasta, db_type,
                                       loc_dict, db_s3_path):
        ''' Download accessions specified in the selected_genera '''
        threads = []
        error_flags = {}
        semaphore = threading.Semaphore(64)
        mutex = TraceLock("download_ref_sequences_from_s3", threading.RLock())

        bucket, key = db_s3_path[5:].split("/", 1)
        accession_dir = os.path.join(self.output_dir_local, db_type, 'accessions')
        command.make_dirs(accession_dir)
        for accession, _taxinfo in accession_dict.items():
            entry = loc_dict.get(accession)
            if not entry:
                continue
            accession_out_file = os.path.join(accession_dir, accession)
            semaphore.acquire()
            t = threading.Thread(
                target=PipelineStepDownloadAccessions.fetch_sequence_for_thread,
                args=[
                    error_flags, accession, accession_out_file, list(entry),
                    bucket, key, semaphore, mutex
                ])
            t.start()
            threads.append(t)
        for t in threads:
            t.join()
        if error_flags:
            raise RuntimeError("Error in getting sequences by accession list.")
        # Combine all the downloaded accessions to a fasta file
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''find "${accession_dir}" -type f | xargs -n 32 -P 1 cat >> "${output_reference_fasta}";''',
                named_args={
                    'accession_dir': os.path.join(accession_dir, ""),  # adds a slash at the end if it doesn't have one
                    'output_reference_fasta': output_reference_fasta
                }
            )
        )

    @staticmethod
    def get_sequence_by_accession_from_file(accession_id, loc_dict, db_file):
        entry = loc_dict.get(accession_id)
        if entry:
            range_start = int(entry[0])
            seq_len = int(entry[1]) + int(entry[2])
            if seq_len <= MAX_ACCESSION_SEQUENCE_LEN:
                db_file.seek(range_start, 0)
                return db_file.read(seq_len)
        return None

    @staticmethod
    def fetch_sequence_for_thread(error_flags, accession, accession_out_file,
                                  entry, bucket, key,
                                  semaphore, mutex):
        ''' fetch sequence from S3 for the specific accession'''
        try:
            assert entry
            range_start, name_length, seq_len = [int(e) for e in entry]
            range_end = range_start + name_length + seq_len - 1
            if seq_len <= MAX_ACCESSION_SEQUENCE_LEN:
                num_retries = 3
                for attempt in range(num_retries):
                    try:
                        s3.fetch_byterange(range_start, range_end, bucket, key, accession_out_file)
                        break
                    except Exception as e:
                        if attempt + 1 < num_retries:  # Exponential backoff
                            time.sleep(1.0 * (4**attempt))
                        else:
                            msg = f"All retries failed for getting sequence by accession ID {accession}: {e}"
                        raise RuntimeError(msg)
        except:
            with mutex:
                if not error_flags:
                    traceback.print_exc()
                error_flags["error"] = 1
        finally:
            semaphore.release()
