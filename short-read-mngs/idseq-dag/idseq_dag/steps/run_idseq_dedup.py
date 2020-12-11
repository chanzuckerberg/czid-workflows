import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.fasta as fasta
import idseq_dag.util.log as log

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.exceptions import InsufficientReadsError


class PipelineStepRunIDSeqDedup(PipelineStep):  # Deliberately not PipelineCountingStep
    """Identifies duplicate reads.

    ```
    idseq-dedup
    -i {input_fasta}
    -o {output_fasta}
    -l 70
    ```

    Require exact match on first 70 nucleotides (`-l 70`) to deem fragments identical.
    Only 70, because sequencer errors increase toward the end of a read.
    For an illustration of this effect, look [here](https://insilicoseq.readthedocs.io/en/latest/iss/model.html).

    Per the idseq-dedup Documentation, available [here](https://github.com/chanzuckerberg/idseq-dedup),
    the idseq-dedup command above will output two or three non-empty files.
    The first output is named exactly as directed via the “-o” option, and contains all representative cluster read IDs.
    The second output with extension “.csv” relates each read ID to its representative cluster read ID.
    For paired end reads, a third output lists the cluster representatives for R2 reads.
    """

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0], 2):
            raise InsufficientReadsError("Insufficient reads before idseq-dedup")

    def run(self):
        input_fas = self.input_files_local[0]
        output_files = self.output_files_local()
        assert len(output_files) == len(input_fas) + 2, f"Context: {input_fas} -> {output_files}."
        output_fas = output_files[:len(input_fas)]
        duplicate_cluster_sizes_path = output_files[-1]
        assert duplicate_cluster_sizes_path.endswith(".tsv"), str(output_files)
        duplicate_clusters_path = output_files[-2]
        assert duplicate_clusters_path.endswith(".csv"), str(output_files)

        # See docstring above for explanation of these options.
        idseq_dedup_params = [
            '-i', input_fas[0], '-o', output_fas[0],
            '-l', '70',
            '-c', duplicate_clusters_path,
            '--cluster-size-output', duplicate_cluster_sizes_path,
        ]
        if len(input_fas) == 2:
            idseq_dedup_params += ['-i', input_fas[1], '-o', output_fas[1]]
        command.execute(
            command_patterns.SingleCommand(
                cmd='idseq-dedup',
                args=idseq_dedup_params
            )
        )

    def count_reads(self):
        self.should_count_reads = True
        # Here we intentionally count unique reads.
        self.counts_dict[self.name] = count.reads_in_group(
            self.output_files_local()[:-2])  # last two outputs are not fastas
