''' Generate GSNAP index given NT '''

import os
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.log as log
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

class PipelineStepGenerateGsnapIndex(PipelineStep):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.upload_results_with_checksum = True

    ''' Generate  gsnap index '''
    def run(self):
        """
          Generate GSNAP index. To be called from idseq-infra
        """
        nt_db = self.input_files_local[0][0]
        output_nt_index_tar = self.output_files_local()[0]
        output_nt_index_parent_dir = os.path.dirname(output_nt_index_tar)
        output_tar_base = os.path.basename(output_nt_index_tar)
        output_nt_index_dir_base = output_tar_base[:-4]
        k = self.additional_attributes.get("k", 16)  # kmer k
        log.write(f"input: {nt_db} output: {output_nt_index_tar}")
        command.execute(
            command_patterns.SingleCommand(
                cmd="gmap_build",
                args=[
                    "-D",
                    output_nt_index_parent_dir,
                    "-d",
                    output_nt_index_dir_base,
                    "-k",
                    k,
                    nt_db
                ]
            )
        )

        output_nt_index_dir = os.path.join(output_nt_index_parent_dir, output_nt_index_dir_base)
        self.additional_output_folders_hidden.append(output_nt_index_dir)

        command.execute(
            command_patterns.SingleCommand(
                cd=output_nt_index_parent_dir,
                cmd="tar",
                args=[
                    "cvf",
                    output_tar_base,
                    output_nt_index_dir_base
                ]
            )
        )

    def count_reads(self):
        ''' Count reads '''
        pass
