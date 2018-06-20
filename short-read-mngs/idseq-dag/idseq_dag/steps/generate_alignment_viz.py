import json
import os
import re
import shelve
import subprocess
import threading
import time
import traceback
from collections import defaultdict

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.taxid_lineage import INVALID_CALL_BASE_ID
import idseq_dag.util.log as log
import idseq_dag.util.command as command
import idseq_dag.util.s3 as s3


class PipelineStepGenerateAlignmentViz(PipelineStep):
    """Pipeline step to generate JSON file for read alignment visualizations to
    be consumed by the web app.
    """

    def count_reads(self):
        pass

    REF_DISPLAY_RANGE = 100
    MAX_SEQ_DISPLAY_SIZE = 6000

    def run(self):
        # Setup
        nt_db = s3.fetch_from_s3(
            self.additional_files["nt_db"],
            self.ref_dir_local,
            allow_s3mi=True)
        nt_loc_db = s3.fetch_from_s3(
            self.additional_files["nt_loc_db"],
            self.ref_dir_local,
            allow_s3mi=True)
        db_type = "nt"  # Only NT supported for now
        # TODO: Design a way to map in/out files more robustly, e.g. by name/type
        annotated_m8 = self.input_files_local[0][1]
        annotated_fasta = self.input_files_local[1][0]
        output_json_dir = os.path.join(self.output_dir_local, "align_viz")

        # Go through annotated_fasta with a db_type (NT/NR match). Infer the
        # family/genus/species info
        read2seq = PipelineStepGenerateAlignmentViz.parse_reads(
            annotated_fasta, db_type)
        log.write(f"Read to Seq dictionary size: {len(read2seq)}")

        db_path = nt_loc_db.replace(".db", "")
        nt_loc_dict = shelve.open(db_path)
        groups, line_count = self.process_reads_from_m8_file(
            annotated_m8, read2seq)

        if nt_db.startswith("s3://"):
            PipelineStepGenerateAlignmentViz.get_sequences_by_accession_list_from_s3(
                groups, nt_loc_dict, nt_db)
        else:
            PipelineStepGenerateAlignmentViz.get_sequences_by_accession_list_from_file(
                groups, nt_loc_dict, nt_db)

        for accession_id, ad in groups.items():
            ad['coverage_summary'] = PipelineStepGenerateAlignmentViz.calculate_alignment_coverage(
                ad)

        result_dict, to_be_deleted = self.populate_reference_sequences(groups)

        # Delete temp files
        def safe_multi_delete(files):
            for f in files:
                try:
                    os.remove(f)
                except:
                    pass

        deleter_thread = threading.Thread(
            target=safe_multi_delete, args=[to_be_deleted])
        deleter_thread.start()

        self.dump_align_viz_json(output_json_dir, db_type, result_dict)

        deleter_thread.join()

        # Write summary file
        summary_msg = f"Read2Seq Size: {len(read2seq)}, M8 lines {line_count}, " \
                  f"{len(groups)} unique accession ids "
        summary_file_name = f"{output_json_dir}.summary"
        with open(summary_file_name, 'w') as summary_f:
            summary_f.write(summary_msg)

    def process_reads_from_m8_file(self, annotated_m8, read2seq):
        # Go through m8 file and infer the alignment info. Grab the fasta
        # sequence, lineage info.
        groups = {}
        line_count = 0
        with open(annotated_m8, 'r') as m8f:
            for line in m8f:
                line_count += 1
                if line_count % 100000 == 0:
                    log.write(f"{line_count} lines in the m8 file processed.")

                line_columns = line.rstrip().split("\t")
                read_id = line_columns[0]
                seq_info = read2seq.get(read_id)
                if seq_info:
                    accession_id = line_columns[1]
                    metrics = line_columns[2:]
                    # "ad" is short for "accession_dict" aka "accession_info"
                    ad = groups.get(accession_id, {'reads': []})
                    sequence, ad['family_id'], ad['genus_id'], ad[
                        'species_id'] = seq_info

                    ref_start = int(metrics[-4])
                    ref_end = int(metrics[-3])
                    if ref_start > ref_end:  # SWAP
                        ref_start, ref_end = ref_end, ref_start
                    ref_start -= 1

                    prev_start = ref_start - self.REF_DISPLAY_RANGE
                    if prev_start < 0:
                        prev_start = 0
                    post_end = ref_end + self.REF_DISPLAY_RANGE
                    markers = prev_start, ref_start, ref_end, post_end
                    ad['reads'].append([read_id, sequence, metrics, markers])
                    base_url = "https://www.ncbi.nlm.nih.gov/nuccore"
                    ad['ref_link'] = f"{base_url}/{accession_id}?report=fasta"
                    groups[accession_id] = ad

        log.write(f"{line_count} lines in the m8 file")
        log.write(f"{len(groups)} unique accession ids")
        return groups, line_count

    def populate_reference_sequences(self, groups):
        result_dict = {}
        to_be_deleted = []
        error_count = 0  # Cap max errors

        # "ad" is short for "accession_dict" aka "accession_info"
        for accession_id, ad in groups.items():
            try:
                tmp_file = f'accession-{accession_id}'
                if ad['ref_seq_len'] <= self.MAX_SEQ_DISPLAY_SIZE and 'ref_seq' not in ad:
                    if ad['ref_seq_len'] == 0:
                        ad['ref_seq'] = "REFERENCE SEQUENCE NOT FOUND"
                    else:
                        with open(tmp_file, "rb") as tf:
                            ad['ref_seq'] = tf.read()
                        to_be_deleted.append(tmp_file)

                if 'ref_seq' in ad:
                    ref_seq = ad['ref_seq']
                    for read in ad['reads']:
                        prev_start, ref_start, ref_end, post_end = read[3]
                        read[3] = [
                            ref_seq[prev_start:ref_start],
                            ref_seq[ref_start:ref_end],
                            ref_seq[ref_end:post_end]
                        ]
                else:
                    # The reference sequence is too long to read entirely in RAM,
                    # so we only read the mapped segments.
                    with open(tmp_file, "rb") as tf:
                        for read in ad['reads']:
                            prev_start, ref_start, ref_end, post_end = read[3]
                            tf.seek(prev_start, 0)
                            segment = tf.read(post_end - prev_start)
                            read[3] = [
                                segment[0:(ref_start - prev_start)],
                                segment[(ref_start - prev_start):(
                                    ref_end - prev_start)],
                                segment[(ref_end - prev_start):(
                                    post_end - prev_start)]
                            ]
                    to_be_deleted.append(tmp_file)
                if ad['ref_seq_len'] > self.MAX_SEQ_DISPLAY_SIZE:
                    ad['ref_seq'] = '...Reference Seq Too Long ...'
            except:
                ad['ref_seq'] = "ERROR ACCESSING REFERENCE SEQUENCE FOR ACCESSION " \
                                "ID {}".format(accession_id)
                if error_count == 0:
                    # Print stack trace for first error
                    traceback.print_exc()
                error_count += 1
            finally:
                family_id = ad.pop('family_id')
                genus_id = ad.pop('genus_id')
                species_id = ad.pop('species_id')
                family_dict = result_dict.get(family_id, {})
                genus_dict = family_dict.get(genus_id, {})
                species_dict = genus_dict.get(species_id, {})
                species_dict[accession_id] = ad
                genus_dict[species_id] = species_dict
                family_dict[genus_id] = genus_dict
                result_dict[family_id] = family_dict

        if error_count > 10:
            # Fail this many and the job is toast
            msg = "Sorry, could not access reference sequences for over " \
                  "{error_count} accession IDs.".format(error_count=error_count)
            raise RuntimeError(msg)

        return result_dict, to_be_deleted

    def dump_align_viz_json(self, output_json_dir, db_type, result_dict):
        def align_viz_name(tag, lin_id):
            return f"{output_json_dir}/{db_type}.{tag}.{int(lin_id)}.align_viz.json"

        # Generate JSON files for the align_viz folder
        command.execute(f"mkdir -p {output_json_dir}")
        for (family_id, family_dict) in result_dict.items():
            fn = align_viz_name("family", family_id)
            with open(fn, 'w') as out_f:
                json.dump(family_dict, out_f)
            self.additional_files_to_upload.append(fn)

            for (genus_id, genus_dict) in family_dict.items():
                fn = align_viz_name("genus", genus_id)
                with open(fn, 'w') as out_f:
                    json.dump(genus_dict, out_f)
                self.additional_files_to_upload.append(fn)

                for (species_id, species_dict) in genus_dict.items():
                    fn = align_viz_name("species", species_id)
                    with open(fn, 'w') as out_f:
                        json.dump(species_dict, out_f)
                    self.additional_files_to_upload.append(fn)

    @staticmethod
    def parse_reads(annotated_fasta, db_type):
        read2seq = {}
        search_string = f"species_{db_type}"
        adv_search_string = "family_%s:([-\d]+):.*genus_%s:([-\d]+):.*species_%s:(" \
                            "[-\d]+).*NT:[^:]*:(.*)" % (
                                db_type, db_type, db_type)

        with open(annotated_fasta, 'r') as af:
            read_id = ''
            for line in af:
                if line[0] == '>':
                    read_id = line
                else:
                    sequence = line
                    m = re.search("%s:([\d-]*)" % search_string, read_id)
                    if m:
                        species_id = int(m.group(1))
                        if species_id > 0 or species_id < INVALID_CALL_BASE_ID:
                            # Match found
                            ma = re.search(adv_search_string, read_id)
                            if ma:
                                read2seq[ma.group(4).rstrip()] = [
                                    sequence.rstrip(),
                                    ma.group(1),
                                    ma.group(2),
                                    ma.group(3)
                                ]
        return read2seq

    @staticmethod
    def get_sequences_by_accession_list_from_file(accession2seq, nt_loc_dict,
                                                  nt_file):
        with open(nt_file) as ntf:
            for accession_id, accession_info in accession2seq.items():
                ref_seq, seq_name = PipelineStepGenerateAlignmentViz.get_sequence_by_accession_id_ntf(
                    accession_id, nt_loc_dict, ntf)
                accession_info['ref_seq'] = ref_seq
                accession_info['ref_seq_len'] = len(ref_seq)
                accession_info['name'] = seq_name

    @staticmethod
    def get_sequences_by_accession_list_from_s3(accession_id_groups,
                                                nt_loc_dict, nt_s3_path):
        threads = []
        error_flags = {}
        semaphore = threading.Semaphore(64)
        mutex = threading.RLock()
        nt_bucket, nt_key = nt_s3_path[5:].split("/", 1)
        for accession_id, accession_info in accession_id_groups.items():
            semaphore.acquire()
            t = threading.Thread(
                target=PipelineStepGenerateAlignmentViz.
                get_sequence_for_thread,
                args=[
                    error_flags, accession_info, accession_id, nt_loc_dict,
                    nt_bucket, nt_key, semaphore, mutex
                ])
            t.start()
            threads.append(t)
        for t in threads:
            t.join()
        if error_flags:
            raise RuntimeError("Error in getting sequences by accession list.")

    @staticmethod
    def get_sequence_for_thread(error_flags,
                                accession_info,
                                accession_id,
                                nt_loc_dict,
                                nt_bucket,
                                nt_key,
                                semaphore,
                                mutex,
                                seq_count=[0]):  #pylint: disable=dangerous-default-value
        try:
            ref_seq_len, seq_name = PipelineStepGenerateAlignmentViz.get_sequence_by_accession_id_s3(
                accession_id, nt_loc_dict, nt_bucket, nt_key)
            with mutex:
                accession_info['ref_seq_len'] = ref_seq_len
                accession_info['name'] = seq_name
                seq_count[0] += 1
                if seq_count[0] % 100 == 0:
                    msg = f"{seq_count[0]} sequences fetched, most recently " \
                          f"{accession_id}"
                    log.write(msg)
        except:
            with mutex:
                if not error_flags:
                    traceback.print_exc()
                error_flags["error"] = 1
        finally:
            semaphore.release()

    @staticmethod
    def get_sequence_by_accession_id_ntf(accession_id, nt_loc_dict, ntf):
        ref_seq = ''
        seq_name = ''
        entry = nt_loc_dict.get(accession_id)
        if entry:
            range_start = entry[0]
            seq_len = entry[1] + entry[2]
            ntf.seek(range_start, 0)
            seq_name, ref_seq = ntf.read(seq_len).split("\n", 1)
            ref_seq = ref_seq.replace("\n", "")
            seq_name = seq_name.split(" ", 1)[1]
        return ref_seq, seq_name

    @staticmethod
    def get_sequence_by_accession_id_s3(accession_id, nt_loc_dict, nt_bucket,
                                        nt_key):
        try:
            from subprocess import DEVNULL
        except ImportError:
            DEVNULL = open(os.devnull, "r+b")

        seq_len = 0
        seq_name = ''
        entry = nt_loc_dict.get(accession_id)
        if not entry:
            return seq_len, seq_name

        range_start, name_length, seq_len = entry

        accession_file = f'accession-{accession_id}'
        num_retries = 3
        for attempt in range(num_retries):
            try:
                pipe_file = f'pipe-{attempt}-accession-{accession_id}'
                os.mkfifo(pipe_file)
                range_end = range_start + name_length + seq_len - 1
                get_range = f"aws s3api get-object " \
                            f"--range bytes={range_start}-{range_end} " \
                            f"--bucket {nt_bucket} " \
                            f"--key {nt_key} {pipe_file}"
                get_range_proc = subprocess.Popen(
                    get_range, shell=True, stdout=DEVNULL)

                cmd = f"cat {pipe_file} | " \
                      f"tee >(tail -n+2 | tr -d '\n' > {accession_file}) | " \
                      f"head -1"
                seq_name = subprocess.check_output(
                    cmd, executable='/bin/bash', shell=True).split(" ", 1)[1]
                exitcode = get_range_proc.wait()
                assert exitcode == 0, "Error in getting sequence by accession ID."
                seq_len = os.stat(accession_file).st_size
                break
            except:
                if attempt + 1 < num_retries:  # Exponential backoff
                    time.sleep(1.0 * (4**attempt))
                else:
                    msg = "All retries failed for getting sequence by accession ID."
                    raise RuntimeError(msg)
            finally:
                try:
                    os.remove(pipe_file)
                except:
                    pass
        return seq_len, seq_name

    @staticmethod
    def compress_coverage(coverage):
        keys = sorted(coverage.keys())
        if len(keys) <= 1:
            return coverage
        output = {}

        start = keys[0]
        current = start
        val = coverage[start]

        for k in keys[1:]:
            if (k - current) == 1 and coverage[k] == val:
                current = k
            else:
                output[f"{start}-{current}"] = val
                start = k
                current = k
                val = coverage[k]

        output[f"{start}-{current}"] = val
        return output

    @staticmethod
    def calculate_alignment_coverage(alignment_data):
        ref_len = alignment_data['ref_seq_len']
        # Setup. Can be implemented more cleanly.
        coverage = defaultdict(lambda: 0)
        output = {
            'ref_seq_len': ref_len,
            'total_read_length': 0,
            'total_aligned_length': 0,
            'total_mismatched_length': 0,
            'num_reads': 0
        }
        if ref_len == 0:
            return output

        reads = alignment_data['reads']
        for read in reads:
            seq = read[1]
            m8_metrics = read[2]
            ref_start = int(m8_metrics[-4])
            ref_end = int(m8_metrics[-3])
            if ref_start > ref_end:  # SWAP
                ref_start, ref_end = ref_end, ref_start
            ref_start -= 1

            output['total_read_length'] += len(seq)
            output['total_aligned_length'] += (ref_end - ref_start)
            output['total_mismatched_length'] += int(m8_metrics[2])
            output['num_reads'] += 1
            for bp in range(ref_start, ref_end):
                coverage[bp] += 1
        output['distinct_covered_length'] = len(coverage)
        output[
            'coverage'] = PipelineStepGenerateAlignmentViz.compress_coverage(
                coverage)
        return output
