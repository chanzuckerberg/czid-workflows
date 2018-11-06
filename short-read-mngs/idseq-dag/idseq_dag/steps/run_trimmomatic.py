''' Run Trimmomatic '''
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count
import idseq_dag.util.s3 as s3

class PipelineStepRunTrimmomatic(PipelineStep):
    ''' Trimmomatic PipelineStep implementation '''
    def run(self):
        """
        Trim any residual Illumina adapters.
        Discard any reads that become too short.

        See: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
        """
        input_files = self.input_files_local[0][0:2]
        output_files = self.output_files_local()
        is_paired = (len(input_files) == 2)
        adapter_fasta = s3.fetch_from_s3(
            self.additional_files["adapter_fasta"],
            self.ref_dir_local)

        if is_paired:
            paired_arg = "PE"
            output_args = [output_files[0],                  # R1, paired, to be kept
                           f"{output_files[0]}__unpaired",   # R1, no longer paired, to be discarded
                           output_files[1],                  # R2, paired, to be kept
                           f"{output_files[1]}__unpaired"]   # R2, no longer paired, to be discarded
        else:
            paired_arg = "SE"
            output_args = output_files

        cmd = " ".join([
            "java -jar /usr/local/bin/trimmomatic-0.38.jar",
            paired_arg,
            "-phred33",
            *input_files,
            *output_args,
            f"ILLUMINACLIP:{adapter_fasta}:2:30:10",
            # Remove Illumina adapters provided in the fasta file. Initially, look for seed matches
            # allowing maximally *2* mismatches. These seeds will be extended and clipped if in the case of paired end
            # reads a score of *30* is reached, or in the case of single ended reads a
            # score of *10*.
            "MINLEN:75"
            # Discard reads which are less than *75* bases long after these steps.
        ])
        command.execute(cmd)

    def count_reads(self):
        ''' Count reads '''
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])
