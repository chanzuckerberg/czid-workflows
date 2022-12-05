import json
import os
import re
import traceback
from collections import defaultdict
from urllib.parse import urlparse

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from s3quilt import download_chunks

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.lineage import INVALID_CALL_BASE_ID
import idseq_dag.util.log as log
import idseq_dag.util.command as command

from idseq_dag.util.dict import open_file_db_by_extension

class PipelineStepGenerateAlignmentViz(PipelineStep):
    """Pipeline step to generate JSON file for read alignment visualizations to
    be consumed by the web app.
    """

    def count_reads(self):
        pass

    REF_DISPLAY_RANGE = 100
    MAX_SEQ_DISPLAY_SIZE = 6000

    def run(self):
        nt_s3_path = self.additional_attributes["nt_db"]
        nt_loc_db = self.additional_files["nt_loc_db"]
        db_type = "nt"  # Only NT supported for now
        # TODO: Design a way to map in/out files more robustly, e.g. by name/type
        annotated_m8 = self.input_files_local[0][0]
        annotated_nt_fasta = self.input_files_local[1][0]
        annotated_nr_fasta = self.input_files_local[1][2]
        output_json_dir = os.path.join(self.output_dir_local, "align_viz")
        output_longest_reads_dir = os.path.join(self.output_dir_local, "longest_reads")

        # Go through annotated_nt_fasta with a db_type (NT/NR match). Infer the
        # family/genus/species info
        read2seq = PipelineStepGenerateAlignmentViz.parse_reads(
            annotated_nt_fasta, db_type)
        log.write(f"Read to Seq dictionary size: {len(read2seq)}")

        groups, line_count = self.process_reads_from_m8_file(
            annotated_m8, read2seq)

        with open_file_db_by_extension(nt_loc_db, "QII", stringify=False) as nt_loc_dict:
            log.write("Getting sequences by accession list from file...")
            PipelineStepGenerateAlignmentViz.get_sequences_by_accession_list_from_file(
                groups, nt_loc_dict, nt_s3_path)

        for ad in groups.values():
            ad['coverage_summary'] = PipelineStepGenerateAlignmentViz.calculate_alignment_coverage(
                ad)

        result_dict = self.populate_reference_sequences(groups)

        self.dump_align_viz_json(output_json_dir, db_type, result_dict)

        read2seq = {}  # to free up memory
        command.make_dirs(output_longest_reads_dir)
        PipelineStepGenerateAlignmentViz.output_n_longest_reads("nt", annotated_nt_fasta, output_longest_reads_dir)
        PipelineStepGenerateAlignmentViz.output_n_longest_reads("nr", annotated_nr_fasta, output_longest_reads_dir)

        # Write summary file
        summary_msg = f"Read2Seq Size: {len(read2seq)}, M8 lines {line_count}, " \
            f"{len(groups)} unique accession ids "
        summary_file_name = f"{output_json_dir}.summary"
        with open(summary_file_name, 'w') as summary_f:
            summary_f.write(summary_msg)

    @staticmethod
    def output_n_longest_reads(db_type: str, annotated_fasta: str, output_longest_reads_dir: str, n=5):
        n_longest = {}
        for level in ["family", "genus", "species"]:
            n_longest[level] = defaultdict(list)

        read2seq = PipelineStepGenerateAlignmentViz.parse_reads(annotated_fasta, db_type)
        for read_id, seq_info in read2seq.items():
            ids = {}
            sequence, ids["family"], ids["genus"], ids["species"] = seq_info
            read = SeqRecord(Seq(sequence), id=read_id, description="")
            for level in ["family", "genus", "species"]:
                n_longest_reads = n_longest[level][ids[level]]
                duplicate = False
                for i, r in enumerate(n_longest_reads):
                    if read.seq == r.seq:
                        duplicate = True
                        break
                    if len(read.seq) > len(r.seq):
                        n_longest[level][ids[level]] = n_longest_reads[:i] + [read] + n_longest_reads[i:n - 1]
                        break
                else:
                    if len(n_longest[level][ids[level]]) < n and not duplicate:
                        n_longest[level][ids[level]].append(read)

        for level in n_longest:
            for taxid, sequences in n_longest[level].items():
                fn = f"{output_longest_reads_dir}/{db_type}.{level}.{taxid}.longest_5_reads.fasta"
                with open(fn, "w") as f:
                    writer = FastaWriter(f, wrap=None)
                    writer.write_file(sequences)

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
        error_count = 0  # Cap max errors

        # "ad" is short for "accession_dict" aka "accession_info"
        for accession_id, ad in groups.items():
            try:
                ref_seq = ad.get('ref_seq')
                for read in ad['reads']:
                    prev_start, ref_start, ref_end, post_end = read[3]
                    if ref_seq:
                        read[3] = [
                            ref_seq[prev_start:ref_start],
                            ref_seq[ref_start:ref_end],
                            ref_seq[ref_end:post_end]
                        ]
                    else:
                        read[3] = ['', '', '']

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

        return result_dict

    def dump_align_viz_json(self, output_json_dir, db_type, result_dict):
        def align_viz_name(tag, lin_id):
            return f"{output_json_dir}/{db_type}.{tag}.{int(lin_id)}.align_viz.json"

        # Generate JSON files for the align_viz folder
        command.make_dirs(output_json_dir)
        for (family_id, family_dict) in result_dict.items():
            fn = align_viz_name("family", family_id)
            with open(fn, 'w') as out_f:
                json.dump(family_dict, out_f)

            for (genus_id, genus_dict) in family_dict.items():
                fn = align_viz_name("genus", genus_id)
                with open(fn, 'w') as out_f:
                    json.dump(genus_dict, out_f)

                for (species_id, species_dict) in genus_dict.items():
                    fn = align_viz_name("species", species_id)
                    with open(fn, 'w') as out_f:
                        json.dump(species_dict, out_f)
        self.additional_output_folders_hidden.append(output_json_dir)

    @staticmethod
    def parse_reads(annotated_fasta, db_type):
        read2seq = {}
        search_string = f"species_{db_type}"
        adv_search_string = r"family_%s:([-\d]+):.*genus_%s:([-\d]+):.*species_%s:(" \
                            r"[-\d]+).*NT:[^:]*:(.*)" % (
                                db_type, db_type, db_type)

        with open(annotated_fasta, 'r') as af:
            read_id = ''
            for line in af:
                if line[0] == '>':
                    read_id = line
                else:
                    sequence = line
                    m = re.search(r"%s:([\d-]*)" % search_string, read_id)
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
                                                  nt_s3_path):
        parsed = urlparse(nt_s3_path)
        accession_ids = [a_id for a_id in accession2seq.keys() if a_id in nt_loc_dict]
        accession_ranges = [nt_loc_dict[a_id] for a_id in accession_ids]
        sequences = download_chunks(
            parsed.hostname,
            parsed.path[1:],
            (s + hl for s, hl, _ in accession_ranges),
            (sl for _, _, sl in accession_ranges),
        )

        """
        we want to create the equivelent pairs to:

          for accession_id, data in zip(accession_ids, data)

        but we want to remove data from the `sequences` list as we add it to `accession2seq`
        to avoid using too much memory double-storing the sequences. `pop` removes items
        from the end so calling it on a whole list is the equivilent to iterating through it
        in reverse. Since order doesn't matter it just needs to preserve the pairs between
        `accession_ids` and `sequences` all we need to do is reverse iterate through
        `accession_ids` to create the same result.
        """
        for accession_id in reversed(accession_ids):
            data = sequences.pop()
            ref_seq = data.replace("\n", "")
            accession2seq[accession_id]['ref_seq'] = ref_seq
            accession2seq[accession_id]['ref_seq_len'] = len(ref_seq)

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
