from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command


class PipelineStepRunCDHitDup(PipelineStep):
    def run(self):
        for f in self.output_files_local():
            command.execute("date > %s" % f)
