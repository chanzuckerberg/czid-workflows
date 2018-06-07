from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command


class PipelineStepRunGsnapFilter(PipelineStep):
    def run(self):
        for f in self.output_files_local():
            command.execute("date > %s" % f)
