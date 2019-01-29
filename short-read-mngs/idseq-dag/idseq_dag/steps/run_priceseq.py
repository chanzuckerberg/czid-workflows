''' Run PriceSeq Filter '''
import os

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count
import idseq_dag.util.fasta as fasta

class PipelineStepRunPriceSeq(PipelineStep):
    ''' PriceSeq PipelineStep implementation '''
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
                log.write(f"0 reads left after PriceSeqFilter, \
                            converting input to FASTA and re-running")
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

    def run_priceseqfilter(self, in_files, out_files, is_paired, file_type):
        params = ["PriceSeqFilter", '-a', '12', '-rnf', '90', '-log', 'c']
        if is_paired:
            params.extend([
                '-fp', in_files[0], in_files[1], '-op', out_files[0],
                out_files[1]
            ])
        else:
            params.extend(['-f', in_files[0], '-o', out_files[0]])
        if file_type != "fasta":  # Default fastq. Explicitly specify fasta.
            params.extend(['-rqf', '85', '0.98'])
        cmd = " ".join(params)
        command.execute(cmd)


    def count_reads(self):
        ''' Count reads '''
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

