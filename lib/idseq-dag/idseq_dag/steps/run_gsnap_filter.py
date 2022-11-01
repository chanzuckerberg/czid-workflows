import os
import subprocess
from idseq_dag.engine.pipeline_step import PipelineCountingStep
from idseq_dag.exceptions import InsufficientReadsError
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.convert as convert
import idseq_dag.util.log as log
import idseq_dag.util.count as count
from idseq_dag.util.s3 import fetch_reference


class PipelineStepRunGsnapFilter(PipelineCountingStep):
    """ Regardless of specified “host” organism, it is essential to remove all potentially-human
    sequences for privacy reasons. Thus, a final GSNAP alignment is performed against the human
    genome for samples from all host types.

    ```
    gsnapl
    -A sam
    --batch=0
    --use-shared-memory=0
    --gmap-mode=all
    --npaths=1
    --ordered
    -t 32
    --max-mismatches=40
    -D {gsnap_base_dir}
    -d {gsnap_index_name}
    -o {output_sam_file}
    {input_fas}
    ```

    Two input FASTAs means paired reads.
    The GSNAP documentation can be found [here](http://research-pub.gene.com/gmap/src/README).
    """
    # Two input FASTAs means paired reads.

    def input_fas(self):
        return self.input_files_local[0][0:2]

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_fas(), 1):
            raise InsufficientReadsError("There was an insufficient number of reads in the sample after the host and quality filtering steps.")

    def validate_output_files(self):
        if not count.files_have_min_reads(self.output_files_local(), 1):
            raise InsufficientReadsError("There was an insufficient number of reads in the sample after the host and quality filtering steps.")

    def run(self):
        input_fas = self.input_fas()
        output_fas = self.output_files_local()
        output_sam_file = os.path.join(self.output_dir_local,
                                       self.additional_attributes["output_sam_file"])
        self.additional_output_files_hidden.append(output_sam_file)

        genome_dir = fetch_reference(self.additional_files["gsnap_genome"],
                                     self.ref_dir_local,
                                     allow_s3mi=True,
                                     auto_untar=True)
        gsnap_base_dir = os.path.dirname(genome_dir)
        gsnap_index_name = os.path.basename(genome_dir)
        # Hack to determine gsnap vs gsnapl
        error_message = subprocess.run(
            ['gsnapl', '-D', gsnap_base_dir, '-d', gsnap_index_name],
            input='>'.encode('utf-8'),
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE
        ).stderr
        gsnap_exe = "gsnap" if 'please run gsnap instead' in error_message.decode('utf-8') else "gsnapl"
        # Run Gsnap
        gsnap_params = [
            '-A', 'sam', '--batch=0', '--use-shared-memory=0',
            '--gmap-mode=all', '--npaths=1', '--ordered', '-t', 32,
            '--max-mismatches=40', '-D', gsnap_base_dir, '-d', gsnap_index_name,
            '-o',
            output_sam_file
        ] + input_fas
        command.execute(
            command_patterns.SingleCommand(
                cmd=gsnap_exe,
                args=gsnap_params
            )
        )
        log.write("Finished GSNAP alignment.")

        # Extract out unmapped files from sam
        if len(input_fas) == 2:
            convert.generate_unmapped_pairs_from_sam(
                output_sam_file, output_fas)
        else:
            convert.generate_unmapped_singles_from_sam(
                output_sam_file, output_fas[0])
        self.validate_output_files()
