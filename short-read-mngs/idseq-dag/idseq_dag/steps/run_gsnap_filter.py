import os
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.convert as convert
import idseq_dag.util.log as log
import idseq_dag.util.count as count
from idseq_dag.util.s3 import fetch_from_s3


class PipelineStepRunGsnapFilter(PipelineStep):
    """GSNAP is an aligner we use to filter out reads mapping to the host
    genome. We add this additional step after STAR/Bowtie-based filtering to
    increase our sensitivity to host filtering. Currently only runs on the
    human host with a chimp reference genome to increase sensitivity to human
    matches. Two input FASTAs means paired reads.

    http://research-pub.gene.com/gmap/
    """
    def run(self):
        input_fas = self.input_files_local[0][0:2]
        output_fas = self.output_files_local()
        output_sam_file = os.path.join(self.output_dir_local,
                                       self.additional_attributes["output_sam_file"])
        self.additional_files_to_upload.append(output_sam_file)

        genome_dir = fetch_from_s3(self.additional_files["gsnap_genome"],
                                   self.ref_dir_local,
                                   allow_s3mi=True,
                                   auto_untar=True)
        gsnap_base_dir = os.path.dirname(genome_dir)
        gsnap_index_name = os.path.basename(genome_dir)
        # Run Gsnap
        gsnap_params = [
            'gsnapl', '-A sam', '--batch=0', '--use-shared-memory=0',
            '--gmap-mode=all', '--npaths=1', '--ordered', '-t 32',
            '--max-mismatches=40', '-D', gsnap_base_dir, '-d', gsnap_index_name,
            '-o',
            output_sam_file
        ] + input_fas
        command.execute(" ".join(gsnap_params))
        log.write("Finished GSNAP alignment.")

        # Extract out unmapped files from sam
        if len(input_fas) == 2:
            convert.generate_unmapped_pairs_from_sam(
                output_sam_file, output_fas)
        else:
            convert.generate_unmapped_singles_from_sam(
                output_sam_file, output_fas[0])

    def count_reads(self):
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])
