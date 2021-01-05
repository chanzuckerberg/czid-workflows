from os.path import basename, dirname, join

import idseq_dag.util.fasta as fasta
import idseq_dag.util.m8 as m8

from idseq_dag.engine.pipeline_step import PipelineCountingStep
from idseq_dag.util.idseq_dedup_clusters import parse_clusters_file
from idseq_dag.util.count import READ_COUNTING_MODE, ReadCountingMode
from idseq_dag.util.parsing import BlastnOutput6NTRerankedReader

UNMAPPED_HEADER_PREFIX = '>NR::NT::'


class PipelineStepGenerateAnnotatedFasta(PipelineCountingStep):
    """
    Generate annotated fasta and unidentified fasta
    """

    def _annotated_fasta(self):
        return self.output_files_local()[0]

    def _unidentified_fasta(self):
        return self.output_files_local()[1]

    def _unique_unidentified_fasta(self):
        return join(
            dirname(self._unidentified_fasta()),
            'unique_' + basename(self._unidentified_fasta())
        )

    def run(self):
        merged_fasta = self.input_files_local[0][-1]
        gsnap_m8 = self.input_files_local[1][1]
        rapsearch2_m8 = self.input_files_local[2][1]

        # See app/lib/dags/postprocess.json.jbuilder in idseq-web
        if len(self.input_files_local) == 5:
            assert READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL
            idseq_dedup_clusters = self.input_files_local[3][0]
            # NOTE: this will load the set of all original read headers, which
            # could be several GBs in the worst case.
            clusters_dict = parse_clusters_file(idseq_dedup_clusters)
        else:
            clusters_dict = None

        annotated_fasta = self._annotated_fasta()
        unidentified_fasta = self._unidentified_fasta()
        unique_unidentified_fasta = self._unique_unidentified_fasta()
        self.annotate_fasta_with_accessions(merged_fasta, gsnap_m8, rapsearch2_m8, annotated_fasta)
        self.generate_unidentified_fasta(
            annotated_fasta,
            unidentified_fasta,
            clusters_dict,
            unique_unidentified_fasta
        )
        if clusters_dict:
            self.additional_output_files_visible.append(unique_unidentified_fasta)

    def count_reads(self):
        # The webapp expects this count to be called "unidentified_fasta"
        # To keep backwards compatibility, we use unique reads. See
        # generate_unidentified_fasta.
        super()._count_reads_work(
            cluster_key=PipelineStepGenerateAnnotatedFasta.old_read_name,
            counter_name="unidentified_fasta",
            fasta_files=[self._unique_unidentified_fasta()]
        )

    @staticmethod
    def annotate_fasta_with_accessions(merged_input_fasta, nt_m8, nr_m8, output_fasta):
        def get_map(blastn_6_path):
            with open(blastn_6_path) as blastn_6_f:
                return {row["qseqid"]: row["sseqid"] for row in BlastnOutput6NTRerankedReader(blastn_6_f, filter_invalid=True)}

        nt_map = get_map(nt_m8)
        nr_map = get_map(nr_m8)

        with open(merged_input_fasta, 'r', encoding='utf-8') as input_fasta_f:
            with open(output_fasta, 'w') as output_fasta_f:
                sequence_name = input_fasta_f.readline()
                sequence_data = input_fasta_f.readline()

                # TODO: (tmorse) fasta parsing
                while sequence_name and sequence_data:
                    read_id = sequence_name.rstrip().lstrip('>')
                    # Need to annotate NR then NT in this order for alignment viz
                    # Its inverse is old_read_name()
                    new_read_name = "NR:{nr_accession}:NT:{nt_accession}:{read_id}".format(
                        nr_accession=nr_map.get(read_id, ''),
                        nt_accession=nt_map.get(read_id, ''),
                        read_id=read_id)
                    output_fasta_f.write(">%s\n" % new_read_name)
                    output_fasta_f.write(sequence_data)
                    sequence_name = input_fasta_f.readline()
                    sequence_data = input_fasta_f.readline()

    @staticmethod
    def old_read_name(new_read_name):
        # Inverse of new_read_name creation above.  Needed to cross-reference to original read_id
        # in order to identify all duplicate reads for this read_id.
        return new_read_name.split(":", 4)[-1]

    def generate_unidentified_fasta(
        self,
        input_fa,
        output_fa,
        clusters_dict=None,
        unique_output_fa=None
    ):
        """
        Generates files with all unmapped reads. If COUNT_ALL, which was added
        in v4, then include non-unique reads extracted upstream by idseq-dedup.

        unique_output_fa exists primarily for counting. See count_reads above.
        """
        unique_output_file = open(unique_output_fa, "w") if clusters_dict else None
        with open(output_fa, "w") as output_file:
            for read in fasta.iterator(input_fa):
                if not read.header.startswith(UNMAPPED_HEADER_PREFIX):
                    continue

                output_file.write(read.header + "\n")
                output_file.write(read.sequence + "\n")
                if unique_output_file:
                    unique_output_file.write(read.header + "\n")
                    unique_output_file.write(read.sequence + "\n")

                if clusters_dict:
                    # get inner part of header like
                    # '>NR::NT::NB501961:14:HM7TLBGX2:4:23511:18703:20079/2'
                    line = read.header
                    header_suffix = ""
                    if line[-2:-1] == "/":  # /1 or /2
                        line, header_suffix = line[:-2], line[-2:]
                        assert header_suffix in ('/1', '/2')
                        assert len(read.header) == len(line) + len(header_suffix)

                    key = line.split(UNMAPPED_HEADER_PREFIX)[1]
                    other_keys = clusters_dict[key][1:]  # key should always be present
                    for other_key in other_keys:
                        other_header = UNMAPPED_HEADER_PREFIX + other_key + header_suffix
                        output_file.write(other_header + "\n")
                        output_file.write(read.sequence + "\n")  # write duplicate seq
