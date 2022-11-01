import csv

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.log as log

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.exceptions import InsufficientReadsError, InvalidInputFileError
from idseq_dag.util.czid_dedup_clusters import parse_clusters_file
from idseq_dag.util.count import save_duplicate_cluster_sizes


class PipelineStepRunCZIDDedup(PipelineStep):  # Deliberately not PipelineCountingStep
    """Identifies duplicate reads.

    ```
    czid-dedup
    -i {input_fasta}
    -o {output_fasta}
    -l 70
    ```

    Require exact match on first 70 nucleotides (`-l 70`) to deem fragments identical.
    Only 70, because sequencer errors increase toward the end of a read.
    For an illustration of this effect, look [here](https://insilicoseq.readthedocs.io/en/latest/iss/model.html).

    Per the czid-dedup Documentation, available [here](https://github.com/chanzuckerberg/czid-dedup),
    the czid-dedup command above will output two or three non-empty files.
    The first output is named exactly as directed via the “-o” option, and contains all representative cluster read IDs.
    The second output with extension “.csv” relates each read ID to its representative cluster read ID.
    For paired end reads, a third output lists the cluster representatives for R2 reads.
    """

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0], 2):
            raise InsufficientReadsError("There was an insufficient number of reads in the sample after the host and quality filtering steps.")

    def run(self):
        input_fas = self.input_files_local[0]
        output_files = self.output_files_local()
        assert len(output_files) == len(input_fas) + 2, f"Context: {input_fas} -> {output_files}."
        output_fas = output_files[:len(input_fas)]
        duplicate_cluster_sizes_path = output_files[-1]
        assert duplicate_cluster_sizes_path.endswith(".tsv"), str(output_files)
        duplicate_clusters_path = output_files[-2]
        assert duplicate_clusters_path.endswith(".csv"), str(output_files)

        duplicate_clusters_raw_path = 'clusters_raw.csv'
        # See docstring above for explanation of these options.
        czid_dedup_params = [
            '-i', input_fas[0], '-o', output_fas[0],
            '-l', '70',
            '-c', duplicate_clusters_raw_path,
        ]
        if len(input_fas) == 2:
            czid_dedup_params += ['-i', input_fas[1], '-o', output_fas[1]]
        command.execute(
            command_patterns.SingleCommand(
                cmd='czid-dedup',
                args=czid_dedup_params
            )
        )

        read_ids = set()
        # Add leading single quote to every cell in the CSV that starts with a special character
        #   to guard against potential csv injection
        with open(duplicate_clusters_raw_path) as r, open(duplicate_clusters_path, 'w') as w:
            for row in csv.reader(r):
                if row[1] in read_ids:
                    raise InvalidInputFileError("The input file(s) has duplicate read IDs.")
                read_ids.add(row[1])
                csv.writer(w).writerow(c if c[0].isalnum() else f"'{c}" for c in row)

        # Emit cluster sizes.  One line per cluster.  Format "<cluster_size> <cluster_read_id>".
        # This info is loaded in multiple subsequent steps using m8.load_duplicate_cluster_sizes,
        # and used to convert unique read counts to original read counts, and also to compute
        # per-taxon DCRs emitted alongside taxon_counts.
        log.write("parsing duplicate clusters file")
        clusters_dict = parse_clusters_file(duplicate_clusters_path)
        log.write("saving duplicate cluster sizes")
        save_duplicate_cluster_sizes(duplicate_cluster_sizes_path, clusters_dict)
        log.write("saved duplicate cluster sizes")

    def count_reads(self):
        self.should_count_reads = True
        # Here we intentionally count unique reads.
        self.counts_dict[self.name] = count.reads_in_group(
            self.output_files_local()[:-2])  # last two outputs are not fastas
