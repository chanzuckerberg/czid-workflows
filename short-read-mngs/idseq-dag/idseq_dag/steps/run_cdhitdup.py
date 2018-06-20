from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count

class PipelineStepRunCDHitDup(PipelineStep):
    '''
    CD-HIT-DUP is used to identify duplicates from single or paired reads.
    Two FASTA inputs means paired reads.
    See: http://weizhongli-lab.org/cd-hit/
    '''
    def run(self):
        ''' Invoking cd-hit-dup '''
        input_fas = self.input_files_local[0]
        output_fas = self.output_files_local()
        cdhitdup_params = [
            'cd-hit-dup', '-i', input_fas[0], '-o', output_fas[0],
            '-e', '0.05', '-u', '70'
        ]
        if len(input_fas) == 2:
            cdhitdup_params += ['-i2', input_fas[1], '-o2', output_fas[1]]
        command.execute(" ".join(cdhitdup_params))

    def count_reads(self):
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

