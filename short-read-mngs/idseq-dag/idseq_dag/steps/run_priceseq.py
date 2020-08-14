''' Run PriceSeq Filter '''
import os
from subprocess import run

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.exceptions import InsufficientReadsError
import idseq_dag.util.log as log
import idseq_dag.util.count as count
import idseq_dag.util.fasta as fasta

class PipelineStepRunPriceSeq(PipelineStep):
    """ Removes low-quality reads to reduce their impact on downstream alignment and analysis.

    The default input is a .fastq file, which enables the use of -rqf parameter for filtering on base quality.

    ```
    PriceSeqFilter
    -a 12
    -rnf 90
    -log c
    -fp {input files}
    -op {output files}
    -rqf 85 0.98
    ```

    Per the PriceSeq Documentation, available [here](http://derisilab.ucsf.edu/software/price/PriceDocumentation130506/independentQualityFilter.html),
    this command uses 12 threads to

    1. filter out pairs of reads if either has an unacceptably high number of uncalled nucleotides (N's),
    requiring that 90% of nucleotides are called per read.
    2. filter out sequences with an unacceptably high number of low-quality nucleotides, as defined by
    the provided quality score, requiring that 85% of nucleotides have a probability of being correct of
    greater than 0.98.
    """
    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0][0:2], 1):
            raise InsufficientReadsError("Insufficient reads")

    def run(self):
        """PriceSeqFilter is used to filter input data based on quality. Two FASTQ
        inputs means paired reads.

        See: http://derisilab.ucsf.edu/software/price/
        """
        input_files = self.input_files_local[0][0:2]
        output_files = self.output_files_local()
        is_paired = (len(input_files) == 2)

        # PriceSeqFilter determines input type based on extension. It will
        # throw an exception if output extension doesn't match input
        # extension.
        file_type = fasta.input_file_type(input_files[0])

        price_out = [
            f"{f}_priceseqfilter_output.{file_type}" for f in input_files
        ]
        self.run_priceseqfilter(input_files, price_out, is_paired, file_type)

        # Some FASTQ files come with very low quality scores, e.g. some PacBio files
        # are entirely 0s for quality score, so in the case that everything gets
        # filtered out by PriceSeqFilter, try running again w/o quality scores.
        # We accomplish this by just converting to FASTA.
        if file_type == 'fastq':
            convert_and_rerun = False
            for file in price_out:
                statinfo = os.stat(file)
                if statinfo.st_size == 0:
                    # if any output files are size zero, then we need to rerun
                    convert_and_rerun = True
                    break

            if convert_and_rerun:
                log.write("0 reads left after PriceSeqFilter, converting input to FASTA and re-running")
                step = "FASTQ to FASTA conversion"
                log.write(f"Starting {step}...")
                convert_out = [
                    f"{f}_converted.fasta" for f in input_files
                ]
                fasta.fq2fa(input_files[0], convert_out[0])
                if is_paired:
                    fasta.fq2fa(input_files[1], convert_out[1])
                log.write(f"Finished {step}.")
                file_type = 'fasta'
                price_out = [
                    f"{f}_priceseqfilter_output.{file_type}" for f in input_files
                ]
                self.run_priceseqfilter(convert_out, price_out, is_paired, file_type)

        # After PriceSeqFilter, all files should be in FASTA format
        if file_type != 'fasta':
            step = "FASTQ to FASTA conversion"
            log.write(f"Starting {step}...")
            fasta.fq2fa(price_out[0], output_files[0])
            if is_paired:
                fasta.fq2fa(price_out[1], output_files[1])
            log.write(f"Finished {step}.")
        else:
            step = "Multi-line FASTA to single-line FASTA conversion"
            log.write(f"Starting {step}...")
            fasta.multilinefa2singlelinefa(price_out[0], output_files[0])
            if is_paired:
                fasta.multilinefa2singlelinefa(price_out[1], output_files[1])
            log.write(f"Finished {step}.")

    def run_priceseqfilter(self, in_files, out_files, is_paired, file_type="fastq"):
        command = ["PriceSeqFilter", "-a", "12", "-rnf", "90", "-log", "c"]
        params = {"fasta": [], "fastq": ["-rqf", "85", "0.98"]}
        os.symlink(in_files[0], "priceseq_in1." + file_type)
        if is_paired:
            assert len(in_files) == 2 and len(out_files) == 2
            os.symlink(in_files[1], "priceseq_in2." + file_type)
            command.extend(["-fp", "priceseq_in1." + file_type, "priceseq_in2." + file_type])
            command.extend(["-op", "priceseq_out1." + file_type, "priceseq_out2." + file_type])
            command.extend(params[file_type])
        else:
            assert len(in_files) == 1 and len(out_files) == 1
            command.extend(["-f", "priceseq_in1." + file_type, "-o", "priceseq_out1." + file_type])
        run(command)
        os.rename("priceseq_out1." + file_type, out_files[0])
        if is_paired:
            os.rename("priceseq_out2." + file_type, out_files[1])

    def count_reads(self):
        ''' Count reads '''
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])
