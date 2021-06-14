import os
import multiprocessing
from idseq_dag.engine.pipeline_step import PipelineCountingStep
from idseq_dag.exceptions import InsufficientReadsError
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.convert as convert
import idseq_dag.util.log as log
from idseq_dag.util.s3 import fetch_reference


class PipelineStepRunBowtie2(PipelineCountingStep):
    """ Removes remaining host reads.

    While STAR provides an initial, rapid removal of host sequences,
    it is possible that some host sequences will make it past this step.
    Therefore, to improve filter sensitivity, an additional host removal step is
    performed using Bowtie2 to clean up any remaining sequences.

    ```
    bowtie2
    -q
    -x {genome_basename}
    -f
    --very-sensitive-local
    -S {output_sam_file}
    --seed random_seed
    -p {number_of_cpus}
    -1 [input R1]
    -2 [input R2]
    ```

    If the input is single-end the `-U [input R, if not paired]` argument is used instead of `-1` and `-2`.
    Bowtie2 documentation can be found [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner)
    """

    def input_fas(self):
        return self.input_files_local[0][0:2]

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_fas(), 1):
            raise InsufficientReadsError("Insufficient reads before bowtie2")

    def run(self):
        input_fas = self.input_fas()
        output_fas = self.output_files_local()
        genome_dir = fetch_reference(
            self.additional_files["bowtie2_genome"],
            self.ref_dir_local,
            allow_s3mi=True,
            auto_untar=True)
        output_sam_file = os.path.join(
            self.output_dir_local,
            self.additional_attributes["output_sam_file"])
        self.additional_output_files_hidden.append(output_sam_file)
        # The file structure looks like
        # "bowtie2_genome/GRCh38.primary_assembly.genome.3.bt2"
        genome_basename = command.glob(f"{genome_dir}/*.bt2*", max_results=1)[0]
        # remove two extensions: ex: hg38_phiX_rRNA_mito_ERCC.3.bt2 -> hg38_phiX_rRNA_mito_ERCC
        genome_basename = os.path.splitext(os.path.splitext(genome_basename)[0])[0]

        bowtie2_params = [
            '-q', '-x', genome_basename, '-f',
            '--very-sensitive-local', '-S', output_sam_file
        ]

        # --seed cannot be used with -p multithreading
        # We have observed the lack of multithreading resulting in
        # severe performance degradation in some cases. So for the
        # time being multithreading is being chosen over determinism.
        # To seed bowtie2 do something similar to:
        # bowtie2_params.extend(['--seed', '4'])
        bowtie2_params.extend(['-p', str(multiprocessing.cpu_count())])

        if len(input_fas) == 2:
            bowtie2_params.extend(['-1', input_fas[0], '-2', input_fas[1]])
        else:
            bowtie2_params.extend(['-U', input_fas[0]])

        # Example:
        # bowtie2 -q -x /mnt/idseq/ref/bowtie2_genome/hg38_phiX_rRNA_mito_ERCC -f \
        #         --very-sensitive-local -S /mnt/idseq/results/589/bowtie2_human.sam \
        #         -p 32 \
        #         -1 /mnt/idseq/results/589/unmapped_human_1.fa -2 /mnt/idseq/results/589/unmapped_human_2.fa
        command.execute(
            command_patterns.SingleCommand(
                cmd='bowtie2',
                args=bowtie2_params
            )
        )
        log.write("Finished Bowtie alignment.")

        if len(input_fas) == 2:
            convert.generate_unmapped_pairs_from_sam(output_sam_file,
                                                     output_fas)
        else:
            convert.generate_unmapped_singles_from_sam(output_sam_file,
                                                       output_fas[0])
