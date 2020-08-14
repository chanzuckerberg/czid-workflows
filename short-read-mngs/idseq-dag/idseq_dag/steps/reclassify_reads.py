import json
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count

class PipelineStepReclassifyReads(PipelineStep):
    '''
        Reclassify reads after alignment
    '''
    def run(self):
        '''
            Dummy implementation. just copy the files over.
            Real thing to be implemented later.
        '''
        input_files = self.input_files_local[0]
        output_files = self.output_files_local()
        for i in range(len(input_files)):
            command.copy_file(input_files[i], output_files[i])

        command.write_text_to_file('1234', output_files[4])

    def count_reads(self):
        pass
