from os.path import basename, dirname, join

import idseq_dag.util.fasta as fasta
import idseq_dag.util.m8 as m8

from idseq_dag.engine.pipeline_step import PipelineCountingStep
from idseq_dag.util.czid_dedup_clusters import parse_clusters_file
from idseq_dag.util.count import READ_COUNTING_MODE, ReadCountingMode
from idseq_dag.util.parsing import BlastnOutput6NTRerankedReader

UNMAPPED_HEADER_PREFIX = '>NR::NT::'


def _old_read_name(new_read_name):
    # Inverse of new_read_name creation above.  Needed to cross-reference to original read_id
    # in order to identify all duplicate reads for this read_id.
    return new_read_name.split(":", 4)[-1]


def _annotate_fasta_with_accessions(
    merged_input_fasta,
    nt_m8,
    nr_m8,
    output_fasta,
    output_unmapped_fasta,
    clusters_dict=None,
    unique_output_fa=None
):
    def get_map(blastn_6_path):
        with open(blastn_6_path) as blastn_6_f:
            return {row["qseqid"]: row["sseqid"] for row in BlastnOutput6NTRerankedReader(blastn_6_f, filter_invalid=True)}

    nt_map = get_map(nt_m8)
    nr_map = get_map(nr_m8)

    unique_output_file = open(unique_output_fa, "w") if clusters_dict else None

    with open(merged_input_fasta, 'r', encoding='utf-8') as input_fasta_f, open(output_fasta, 'w') as output_fasta_f, open(output_unmapped_fasta, "w") as output_unmapped_fasta_f:
        sequence_name = input_fasta_f.readline()
        sequence_data = input_fasta_f.readline()

        # TODO: (tmorse) fasta parsing
        while sequence_name and sequence_data:
            read_id = sequence_name.rstrip().lstrip('>')
            # Need to annotate NR then NT in this order for alignment viz
            # Its inverse is old_read_name()
            nr_accession = nr_map.get(read_id, '')
            nt_accession = nt_map.get(read_id, '')
            new_read_name = ">NR:{nr_accession}:NT:{nt_accession}:{read_id}".format(
                nr_accession=nr_accession,
                nt_accession=nt_accession,
                read_id=read_id)
            output_file = output_fasta_f if nt_accession or nr_accession else output_unmapped_fasta_f
            output_file.write(new_read_name + '\n')
            output_file.write(sequence_data)

            if not (nt_accession or nr_accession) and clusters_dict and unique_output_file:
                output_clusters(
                    output_unmapped_fasta_f,
                    unique_output_file,
                    new_read_name,
                    sequence_data,
                    clusters_dict,
                )
            sequence_name = input_fasta_f.readline()
            sequence_data = input_fasta_f.readline()
    if unique_output_file:
        unique_output_file.close()


def output_clusters(output_unmapped_fasta_f, unique_output_file, read_header, read_sequence, clusters_dict):
    unique_output_file.write(read_header + "\n")
    unique_output_file.write(read_sequence)
    line = read_header
    header_suffix = ""
    if line[-2:-1] == "/":  # /1 or /2
        line, header_suffix = line[:-2], line[-2:]
        assert header_suffix in ('/1', '/2')
        assert len(read_header) == len(line) + len(header_suffix)

    key = line.split(UNMAPPED_HEADER_PREFIX)[1]
    if key not in clusters_dict:
        key = key + header_suffix
        header_suffix = ""

    other_keys = clusters_dict[key][1:]  # key should always be present
    for other_key in other_keys:
        other_header = UNMAPPED_HEADER_PREFIX + other_key + header_suffix
        output_unmapped_fasta_f.write(other_header + "\n")
        output_unmapped_fasta_f.write(read_sequence)  # write duplicate seq


def generate_annotated_fasta(
    pre_alignment_fa_path: str,
    nt_m8_path: str,
    nr_m8_path: str,
    annotated_fasta_path: str,
    unidentified_fasta_path: str,
    unique_unidentified_fasta: str = None,
    duplicate_clusters_path: str = None,
):
    # See app/lib/dags/postprocess.json.jbuilder in idseq-web
    if duplicate_clusters_path:
        assert READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL
        # NOTE: this will load the set of all original read headers, which
        # could be several GBs in the worst case.
        clusters_dict = parse_clusters_file(duplicate_clusters_path)
    else:
        clusters_dict = None

    _annotate_fasta_with_accessions(
        pre_alignment_fa_path,
        nt_m8_path,
        nr_m8_path,
        annotated_fasta_path,
        unidentified_fasta_path,
        clusters_dict=clusters_dict,
        unique_output_fa=unique_unidentified_fasta,
    )


class PipelineStepGenerateAnnotatedFasta(PipelineCountingStep):
    """
    Generate annotated fasta and unidentified fasta
    """

    _unique_unidentified_fasta = None

    def run(self):
        # See app/lib/dags/postprocess.json.jbuilder in idseq-web
        duplicate_clusters_path = None
        if len(self.input_files_local) == 5:
            duplicate_clusters_path = self.input_files_local[3][0]

        annotated_fasta_path = self.output_files_local()[0]
        unidentified_fasta_path = self.output_files_local()[1]
        unique_unidentified_fasta = join(
            dirname(unidentified_fasta_path),
            'unique_' + basename(unidentified_fasta_path)
        )
        self._unique_unidentified_fasta = unique_unidentified_fasta

        generate_annotated_fasta(
            pre_alignment_fa_path=self.input_files_local[0][-1],
            nt_m8_path=self.input_files_local[1][1],
            nr_m8_path=self.input_files_local[2][1],
            annotated_fasta_path=annotated_fasta_path,
            unidentified_fasta_path=unidentified_fasta_path,
            duplicate_clusters_path=duplicate_clusters_path,
            unique_unidentified_fasta=unique_unidentified_fasta
        )

        if duplicate_clusters_path:
            self.additional_output_files_visible.append(unique_unidentified_fasta)

    def count_reads(self):
        # The webapp expects this count to be called "unidentified_fasta"
        # To keep backwards compatibility, we use unique reads. See
        # generate_unidentified_fasta.
        super()._count_reads_work(
            cluster_key=_old_read_name,
            counter_name="unidentified_fasta",
            fasta_files=[self._unique_unidentified_fasta]
        )
