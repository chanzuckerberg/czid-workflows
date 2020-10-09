from idseq_dag.steps.run_star import PipelineStepRunStar
from idseq_dag.util.count import load_duplicate_cluster_sizes, reads_in_group

class PipelineStepRunStarDownstream(PipelineStepRunStar):
    """ Runs STAR as a step with two input targets:
    (1) sequence files (output from any step)
        --> can be either (a) a single file, or (b) paired files plus a merged file which will be ignored
    (2) validation counts (output from validation step)
        --> first file is a count json; remaining files are sanitized sequence files which we don't need
            because we use the files from (1)
    """
    def __init__(self, *args, **kwargs):
        # Don't collect insert size metrics at this stage
        super().__init__(*args, disable_insert_size_metrics=True, **kwargs)

    def run(self):
        self.sequence_input_files = self.input_files_local[0][:2]
        self.validated_input_counts_file = self.input_files_local[1][0]
        super().run()

    # The code below is copied from PipelineCountingStep
    # Couldn't inherit PipelineCountingStep because PipelineStepRunStar shouldn't inherit from it
    def input_cluster_sizes_path(self):
        # The last last input to PipelineCountingStep is cluster_sizes.tsv
        tsv = self.input_files_local[-1][-1]
        assert tsv.endswith(".tsv"), str(self.input_files_local)
        return tsv

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = reads_in_group(
            file_group=self.output_files_local()[0:2],
            cluster_sizes=load_duplicate_cluster_sizes(self.input_cluster_sizes_path()),
            cluster_key=lambda x: x)
