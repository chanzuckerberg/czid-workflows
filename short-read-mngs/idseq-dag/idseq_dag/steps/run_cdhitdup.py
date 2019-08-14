from idseq_dag.engine.pipeline_step import PipelineStep, InputFileErrors
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count

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

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])
