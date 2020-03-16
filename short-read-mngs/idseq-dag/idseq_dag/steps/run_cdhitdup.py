import os

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.log as log

from idseq_dag.engine.pipeline_step import InputFileErrors, PipelineStep


class PipelineStepRunCDHitDup(PipelineStep):
    """ Removes duplicate reads.

    ```
    cd-hit-dup
    -i {input_fasta}
    -o {output_fasta}
    -e 0.05
    -u 70
    ```

    Per the CDHit Documentation, available [here](https://github.com/weizhongli/cdhit/wiki),
    this command uses a 0.05 threshold for the number of mismatches - indicating that if
    two reads are > 95% similar they will be considered duplicates. It uses only the
    first/last 70 bases of each read to do the analysis on sequence similarity.

    cd-hit-dup will output three or four files. Two of them are the same as the
    output files of CD-HIT: one (named exactly the same as the file name
    specified by the “-o” option) is the cluster (or duplicate) representatives,
    the other is the clustering file (xxx.clstr) relating each duplicate to its
    representative. For paired end reads, another file by the “-o2” option is
    the cluster representatives for R2 reads. The last file (xxx2.clstr)
    contains the chimeric clusters. In this file, the description for each
    chimeric cluster contains cluster ids of its parent clusters from the
    clustering file xxx.clstr.
    """

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0], 2):
            self.input_file_error = InputFileErrors.INSUFFICIENT_READS

    def run(self):
        ''' Invoking cd-hit-dup '''
        input_fas = self.input_files_local[0]
        output_fas = self.output_files_local()
        cdhitdup_params = [
            '-i', input_fas[0], '-o', output_fas[0],
            '-e', '0.05', '-u', '70'
        ]
        if len(input_fas) == 2:
            cdhitdup_params += ['-i2', input_fas[1], '-o2', output_fas[1]]
        command.execute(
            command_patterns.SingleCommand(
                cmd='cd-hit-dup',
                args=cdhitdup_params
            )
        )
        self._add_clstr_files()

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(
            self.output_files_local()[0:2])

    def _add_clstr_files(self):
        output_fas = self.output_files_local()[0]
        clstr_file = output_fas + '.clstr'  # clusters
        clstr_file2 = output_fas + '2.clstr'  # chimeric clusters
        if os.path.isfile(clstr_file) and os.path.isfile(clstr_file2):
            self.additional_output_files_visible += [clstr_file, clstr_file2]
        else:
            log.write(
                f"WARNING: Files not found: {clstr_file} and {clstr_file2}")
