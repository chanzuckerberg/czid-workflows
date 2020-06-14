import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.fasta as fasta

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.exceptions import InsufficientReadsError
from idseq_dag.util.cdhit_clusters import parse_clusters_file
from idseq_dag.util.count import save_cdhit_cluster_sizes


class PipelineStepRunCDHitDup(PipelineStep):  # Deliberately not PipelineCountingStep
    """Identifies duplicate reads.

    ```
    cd-hit-dup
    -i {input_fasta}
    -o {output_fasta}
    -e 0.0
    -u 70
    ```

    Require exact match (-e 0.0) on first 70 nucleotides (-u 70) to deem fragments identical.
    Only 70, because sequencer errors increase toward the end of a read.  For an illustration
    of this effect, look [here](https://insilicoseq.readthedocs.io/en/latest/iss/model.html).

    Per the CDHit Documentation, available [here](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#cdhitdup),
    the cd-hit-dup command above will output two or three non-empty files.  The first output is named exactly as
    directed via the “-o” option, and contains all cluster (or duplicate) representatives.
    The second output with extension “.clstr” relates each duplicate read ID to its
    distinct representative.  For paired end reads, a third output named by the “-o2” option
    lists the cluster representatives for R2 reads.

    (If the “-f” option were to be specified, a second ”.clstr” file would show
    all chimeric reads;  but that output is empty when “-f” is absent.)

    This step also outputs a TSV file mapping each cluster representative
    to its cluster size, which is useful when converting counts of clusters
    to counts of original fragments.
    """

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0], 2):
            raise InsufficientReadsError("Insufficient reads before CDHitDup")

    def run(self):
        ''' Invoking cd-hit-dup '''
        input_fas = self.input_files_local[0]

        output_files = self.output_files_local()
        assert len(output_files) == len(input_fas) + 2, f"Context: {input_fas} -> {output_files}."
        output_fas = output_files[:len(input_fas)]
        cdhit_cluster_sizes_path = output_files[-1]
        assert cdhit_cluster_sizes_path.endswith(".tsv"), str(output_files)
        cdhit_clusters_path = output_files[-2]
        assert cdhit_clusters_path.endswith(".clstr"), str(output_files)

        # See docstring above for explanation of these options.
        cdhitdup_params = [
            '-i', input_fas[0], '-o', output_fas[0],
            '-e', '0.0', '-u', '70'
        ]
        if len(input_fas) == 2:
            cdhitdup_params += ['-i2', input_fas[1], '-o2', output_fas[1]]
        else:
            # This is a workaround for a cdhitdup bug with unpaired end reads.  Basically
            # it doesn't work with unpaired reads.  It emits the correct number of reads,
            # but not the cluster representatives --- instead emits multiple reads from some
            # clusters, and no reads at all for other clusters.  To work around this issue
            # we simply run cdhitdup with a second input argument that is identical to the
            # first, and place the second output in a dummy.  Note that cdhit actually will
            # segfault if you give it /dev/null for the second output file.
            # TODO: Run "cmp dedup1.fa ignore_me_please.fa" each time as a sanity check.
            cdhitdup_params += ['-i2', input_fas[0], '-o2', "ignore_me_please.fa"]
        command.execute(
            command_patterns.SingleCommand(
                cmd='cd-hit-dup',
                args=cdhitdup_params
            )
        )

        # Emit cluster sizes.  One line per cluster.  Format "<cluster_size> <cluster_read_id>".
        # This info is loaded in multiple subsequent steps using m8.load_cdhit_cluster_sizes,
        # and used to convert unique read counts to original read counts, and also to compute
        # per-taxon DCRs emitted alongside taxon_counts.
        clusters_dict = parse_clusters_file(cdhit_clusters_path, output_fas[0])
        save_cdhit_cluster_sizes(cdhit_cluster_sizes_path, clusters_dict)

    def count_reads(self):
        self.should_count_reads = True
        # Here we intentionally count unique reads.
        self.counts_dict[self.name] = count.reads_in_group(
            self.output_files_local()[:-2])  # last two outputs are not fastas
