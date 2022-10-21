from idseq_dag.steps.run_star import PipelineStepRunStar

class PipelineStepRunStarUpstream(PipelineStepRunStar):
    """ Runs STAR as a step with a single input target:
    namely the output from the validation step, which includes
    both sanitized versions of the original sequence files and validation counts.
    """

    def run(self):
        self.sequence_input_files = self.input_files_local[0][1:3]
        self.validated_input_counts_file = self.input_files_local[0][0]
        super().run()
