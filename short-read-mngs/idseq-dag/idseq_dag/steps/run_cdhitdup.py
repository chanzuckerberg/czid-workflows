import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.fasta as fasta

from idseq_dag.engine.pipeline_step import InputFileErrors, PipelineStep
from idseq_dag.util.count import save_cdhit_cluster_sizes, DAG_SURGERY_HACKS_FOR_READ_COUNTING


class PipelineStepRunCDHitDup(PipelineStep):  # Deliberately not PipelineCountingStep
    """ Removes duplicate reads.

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
            self.input_file_error = InputFileErrors.INSUFFICIENT_READS

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
        PipelineStepRunCDHitDup._emit_cluster_sizes(cdhit_cluster_sizes_path, cdhit_clusters_path, output_fas[0])

        # TODO: When the matching idseq-web request is deployed, remove this line, because those would
        # then become bona-fide outputs of the step and thus would not need to be added here.
        if DAG_SURGERY_HACKS_FOR_READ_COUNTING:
            self.additional_output_files_visible.extend([cdhit_clusters_path, cdhit_cluster_sizes_path])

    @staticmethod
    def _emit_cluster_sizes(cdhit_cluster_sizes_path, cdhit_clusters_path, deduped_fasta_path):
        # Emit cluster sizes.  One line per cluster.  Format "<cluster_size> <cluster_read_id>".
        # This info is loaded in multiple subsequent steps using m8.load_cdhit_cluster_sizes,
        # and used to convert unique read counts to original read counts, and also to compute
        # per-taxon DCRs emitted alongside taxon_counts.

        # First identify the cluster representative reads emitted by cd-hit-dup.  Originally we
        # used the ".clstr" output of cd-hit-dup for this, but turns out that for unpaired reads
        # the actual deduped output of cdhit contains different representatives.
        cluster_sizes_dict = {}
        for read in fasta.iterator(deduped_fasta_path):
            # A read header looks someting like
            #
            #    >M05295:357:000000000-CRPNR:1:1101:22051:10534 OPTIONAL RANDOM STUFF"
            #     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            #
            # where the first character on the line is '>' and the read ID (underlined above)
            # extends from '>' to the first whitespace character, not including '>' itself.
            #
            # The fasta iterator already asserts that read.header[0] is '>'.
            #
            read_id = read.header.split(None, 1)[0][1:]
            cluster_sizes_dict[read_id] = None  # not yet known

        def record_cluster_size(cluster_size, emitted_reads_from_cluster, line_number, _read_id):
            assert emitted_reads_from_cluster, f"If this assertion fails, CD-HIT-DUP has forgotten to emit a read for this cluster.  In that case, just use the current read_id as cluster_representative.  Everything will work fine, aside from reduced sensitivity. {line_number}"
            assert len(emitted_reads_from_cluster) == 1, f"If this assertion fails, CD-HIT-DUP has emitted multiple reads from the same cluster.  Feel free to comment out this assertion if that happens a lot in practice.  Everything will run fine, but read counts contributed by that cluster will be exaggerated.  If you want to fix that, make the cluster sizes a float --- divide the actual cluster size by the number of reads emitted for the cluster, i.e. by len(emitted_reads_from_cluster). Probably an even better way of fixing it would be to emit your own fasta based on the .clstr file if that's reliable, or use a tool other than cdhit that doesn't have this bug.  {line_number}: {emitted_reads_from_cluster}"
            cluster_representative = emitted_reads_from_cluster.pop()
            assert cluster_representative in cluster_sizes_dict, "If this fails it's our bug here."
            cluster_sizes_dict[cluster_representative] = cluster_size

        # Example input lines that form a cluster:
        #
        #    "0       140nt, >M05295:357:000000000-CRPNR:1:2119:16143:8253... *"
        #    "1       140nt, >M05295:357:000000000-CRPNR:1:1101:22051:10534... at 1:140:1:140/+/100.00%"
        #    "2       140nt, >M05295:357:000000000-CRPNR:1:1102:15401:7483... at 1:140:1:140/+/100.00%"
        #    ...
        #    "2334    140nt, >M05295:357:000000000-CRPNR:1:1102:13405:3483... at 1:140:1:140/+/100.00%"
        #
        # Corresponding output line for that cluster:
        #
        #    "2335    140nt, >M05295:357:000000000-CRPNR:1:2119:16143:8253"
        #
        # Please note that "..." above does not indicate truncation. CD-HIT-DUP appends "..." to the read
        # IDs even if the read IDs have not been truncated.
        #
        # Per CD-HIT-DUP docs, when a "-d" argument is not specified, each read ID is obtained by
        # splitting the corresponding FASTA line on whitespace and taking the first item.  That is also
        # how we do it throughout this pipeline.  The code below assumes (and relies upon) the read IDs
        # being free from whitespace.
        #
        # Furthremore, we expect read IDs to be unique in a sequencing run.  Violating that assumption,
        # if it does not break cdhit itself, might produce slightly bogus numbers, because of how it
        # affects subsampling and per-read DCR correction.  Preserving this uniqueness is why we do not
        # specify a "-d" flag to cdhit, and allow it thus to use the entire read id.
        #
        # Further note that, although CD-HIT-DUP documentation states the read marked with '*' is the
        # one chosen as representative, we have found, particularly with unpaired fasta, the emitted
        # deduplicated read is not the one marked with '*'.  Hence we handle that case correctly below,
        # and put a lot of assertions to make sure problems with cdhit output are detected.
        with open(cdhit_clusters_path, "r") as clusters_file:
            emitted_reads_from_cluster = set()  # set of reads in both dedup1.fa and current cluster; cardinality 1!
            cluster_size = 0
            read_id = None
            line_number = 0
            for line in clusters_file:
                line_number += 1
                if line.startswith(">"):
                    continue
                parts = line.strip().split()
                serial = int(parts[0])
                assert parts[2][0] == ">", line
                assert parts[2].endswith("..."), line
                if serial == 0 and cluster_size > 0:
                    # We've just encountered the first read of a new cluster.  Emit all data held for the old cluster.
                    record_cluster_size(cluster_size, emitted_reads_from_cluster, line_number, read_id)
                    emitted_reads_from_cluster = set()
                    cluster_size = 0
                assert cluster_size == serial, f"{line_number}: {cluster_size}, {serial}, {line}"
                read_id = parts[2][1:-3]
                cluster_size += 1
                if read_id in cluster_sizes_dict:
                    emitted_reads_from_cluster.add(read_id)
            if cluster_size > 0:
                record_cluster_size(cluster_size, emitted_reads_from_cluster, line_number, read_id)

        save_cdhit_cluster_sizes(cdhit_cluster_sizes_path, cluster_sizes_dict)

    def count_reads(self):
        self.should_count_reads = True
        # Here we intentionally count unique reads.
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[:-2])  # last two outputs are not fastas
