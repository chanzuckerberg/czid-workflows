from idseq_dag.engine.pipeline_step import PipelineCountingStep
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.m8 as m8

class PipelineStepGenerateAnnotatedFasta(PipelineCountingStep):
    '''
    generate annotated fasta
    '''
    def _annotated_fasta(self):
        return self.output_files_local()[0]

    def _unidentified_fasta(self):
        return self.output_files_local()[1]

    def run(self):
        ''' annotate fasta '''
        merged_fasta = self.input_files_local[0][-1]
        gsnap_m8 = self.input_files_local[1][1]
        rapsearch2_m8 = self.input_files_local[2][1]
        annotated_fasta = self._annotated_fasta()
        unidentified_fasta = self._unidentified_fasta()
        self.annotate_fasta_with_accessions(merged_fasta, gsnap_m8, rapsearch2_m8, annotated_fasta)
        self.generate_unidentified_fasta(annotated_fasta, unidentified_fasta)

    def count_reads(self):
        # The webapp expects this count to be called "unidentified_fasta"
        super()._count_reads_work(
            cluster_key=PipelineStepGenerateAnnotatedFasta.old_read_name,
            counter_name="unidentified_fasta",
            fasta_files=[self._unidentified_fasta()]
        )

    @staticmethod
    def annotate_fasta_with_accessions(merged_input_fasta, nt_m8, nr_m8, output_fasta):
        def get_map(m8_file):
            return dict((read_id, accession_id)
                        for read_id, accession_id, _percent_id,
                        _alignment_length, _e_value, _bitscore, _line in m8.iterate_m8(
                            m8_file, 0, "annotate_fasta_with_accessions"))

        nt_map = get_map(nt_m8)
        nr_map = get_map(nr_m8)

        with open(merged_input_fasta, 'r', encoding='utf-8') as input_fasta_f:
            with open(output_fasta, 'w') as output_fasta_f:
                sequence_name = input_fasta_f.readline()
                sequence_data = input_fasta_f.readline()

                while sequence_name and sequence_data:
                    read_id = sequence_name.rstrip().lstrip('>')
                    # Need to annotate NR then NT in this order for alignment viz
                    new_read_name = "NR:{nr_accession}:NT:{nt_accession}:{read_id}".format(  # Its inverse is old_read_name()
                        nr_accession=nr_map.get(read_id, ''),
                        nt_accession=nt_map.get(read_id, ''),
                        read_id=read_id)
                    output_fasta_f.write(">%s\n" % new_read_name)
                    output_fasta_f.write(sequence_data)
                    sequence_name = input_fasta_f.readline()
                    sequence_data = input_fasta_f.readline()

    @staticmethod
    def generate_unidentified_fasta(input_fa, output_fa):
        # TODO  remove annotated fasta intermediate file and replace > with : below
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''grep -A 1 '>NR::NT::' "$1" | sed '/^--$/d' > "$2";''',
                args=[
                    input_fa,
                    output_fa
                ]
            )
        )

    @staticmethod
    def old_read_name(new_read_name):
        # Inverse of new_read_name creation above.  Needed to cross-reference to original read_id
        # in order to identify all duplicate reads for this read_id.
        return new_read_name.split(":", 4)[-1]
