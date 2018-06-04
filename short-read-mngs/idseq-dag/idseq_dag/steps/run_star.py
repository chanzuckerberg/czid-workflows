import subprocess
from idseq_dag.engine.pipeline_step import PipelineStep
class PipelineStepRunStar(PipelineStep):
    def run(self):
        input_file_0 = self.input_files_local[0][0]
        input_file_1 = self.input_files_local[0][1]
        output_file_0 = self.output_files_local()[0]
        output_file_1 = self.output_files_local()[1]
        subprocess.check_call("gzip -dc %s |head -20 > %s" % (input_file_0, output_file_0), shell=True)
        subprocess.check_call("gzip -dc %s |head -20 > %s" % (input_file_1, output_file_1), shell=True)


