from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command


class PipelineStepRunStar(PipelineStep):
    def run(self):
        input_file_0 = self.input_files_local[0][0]
        input_file_1 = self.input_files_local[0][1]
        output_file_0 = self.output_files_local()[0]
        output_file_1 = self.output_files_local()[1]
        command.execute("gzip -dc %s |head -20 > %s" % (input_file_0, output_file_0))
        command.execute("gzip -dc %s |head -20 > %s" % (input_file_1, output_file_1))
