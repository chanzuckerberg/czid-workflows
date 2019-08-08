''' Run Trimmomatic '''
from idseq_dag.engine.pipeline_step import PipelineStep, InputFileErrors
import idseq_dag.util.command as command
import idseq_dag.util.count as count
import idseq_dag.util.s3 as s3
import idseq_dag.util.fasta as fasta


class PipelineStepRunTrimmomatic(PipelineStep):
    """ Removes adapter sequences.

    ```
    java -jar /usr/local/bin/trimmomatic-0.38.jar 
    PE|SE 
    -phred33 
    [input_files] 
    [output_files] 
    ILLUMINACLIP:{adapter_fasta}:2:30:10:8:true
    MINLEN:35
    ```

    Where, the adapter_fasta for single-end reads is: __illumina_TruSeq3-SE.fasta__

    ```
    >TruSeq3_IndexedAdapter
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    >TruSeq3_UniversalAdapter
    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
    ```

    And paired-end reads: __illumina_TruSeq3-PE-2_NexteraPE-PE.fasta__

    ```
    >PrefixPE/1
    TACACTCTTTCCCTACACGACGCTCTTCCGATCT
    >PrefixPE/2
    GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
    >PE1
    TACACTCTTTCCCTACACGACGCTCTTCCGATCT
    >PE1_rc
    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
    >PE2
    GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
    >PE2_rc
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC>PrefixNX/1
    AGATGTGTATAAGAGACAG
    >PrefixNX/2
    AGATGTGTATAAGAGACAG
    >Trans1
    TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
    >Trans1_rc
    CTGTCTCTTATACACATCTGACGCTGCCGACGA
    >Trans2
    GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
    >Trans2_rc
    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
    ```

    __note:__ the output reads at this step include full-length reads that have no adapter, 
    plus reads where the adapter has been identified and lopped off, leaving a shorter read 
    (but not shorter than 35nt). In cases where the insert size is small, resulting in adapter 
    read-through, R2 will be a direct reverse complement of R1; the "true" parameter enables 
    these reads to be saved for downstream analysis.
    """
    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0][0:2], 1):
            self.input_file_error = InputFileErrors.INSUFFICIENT_READS

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

        if fasta.input_file_type(input_files[0]) != 'fastq':
            # Not fastq
            for in_file, out_file in zip(input_files, output_files):
                command.execute(f"cp {in_file} {out_file}")
            return

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
            f"ILLUMINACLIP:{adapter_fasta}:2:30:10:8:true",
            # Remove Illumina adapters provided in the fasta file. Initially, look for seed matches
            # allowing maximally *2* mismatches. These seeds will be extended and clipped if in the case of paired end
            # reads a score of *30* is reached, or in the case of single ended reads a
            # score of *10*.
            # additional parameters: minAdapterLength = 8, keepBothReads = true; these are set to require pairs to be
            #    kept even when an adapter read-through occurs and R2 is a direct reverse complement of R1.
            "MINLEN:35"
            # Discard reads which are less than *75* bases long after these steps.
        ])
        command.execute(cmd)

    def count_reads(self):
        ''' Count reads '''
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])
