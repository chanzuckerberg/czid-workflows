''' Generate GSNAP index given NT '''

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.log as log
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

class PipelineStepGenerateRapsearch2Index(PipelineStep):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.upload_results_with_checksum = True

    ''' Generate   RAPSearch index from NR '''
    def run(self):
        """
          Generate Rapsearch2 index. To be called from idseq-infra

        """
        nr_db = self.input_files_local[0][0]
        output_nr_index = self.output_files_local()[0]
        output_nr_info_file = output_nr_index + '.info'
        log.write(f"input: {nr_db} output: {output_nr_index}")
        command.execute(
            command_patterns.SingleCommand(
                cmd="prerapsearch",
                args=[
                    "-d",
                    nr_db,
                    "-n",
                    output_nr_index
                ]
            )
        )
        self.additional_files_to_upload.append(output_nr_info_file)

    def count_reads(self):
        ''' Count reads '''
        pass
