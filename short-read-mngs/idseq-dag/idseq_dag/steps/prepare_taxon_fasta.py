''' Generate Phylogenetic tree '''
import os

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.s3 as s3
import idseq_dag.util.count as count

class PipelineStepPrepareTaxonFasta(PipelineStep):
    ''' 
    Fetch fasta file(s) of reads mapping to a taxid (NT/NR) in a sample.
    To be used as input for phylo_trees.
    '''
    def run(self):
        output_file = self.output_files_local()[0]
        byterange_dict = self.additional_attributes["taxon_byteranges"]

        # Retrieve IDseq taxon fasta files
        partial_fasta_files = []
        for hit_type, byterange in byterange_dict.items():
            first_byte = byterange[0]
            last_byte = byterange[1]
            s3_file = byterange[2]
            local_basename = f"{hit_type}_{os.path.basename(output_file)}.fasta"
            bucket, key = s3.split_identifiers(s3_file)
            local_file = os.path.join(self.output_dir_local, local_basename)
            s3.fetch_byterange(first_byte, last_byte, bucket, key, local_file)
            partial_fasta_files.append(local_file)
        self.fasta_union(partial_fasta_files, output_file)
        for fasta in partial_fasta_files + [output_file]:
            print(f"{count.reads(fasta)} reads in {fasta}")

        # Trim Illumina adapters
        # TODO: consider moving this to the beginning of the main pipeline
        self.trim_adapters_in_place(output_file)

    def count_reads(self):
        pass

    @staticmethod
    def trim_adapters_in_place(local_file):
        local_file_trimmed = os.path.join(os.path.dirname(local_file), "trimmed_" + os.path.basename(local_file))
        command.execute(
            command_patterns.SingleCommand(
                cmd='cutadapt',
                args=[
                    "-a",
                    "AGATCGGAAGAGCACACGTCT",
                    "-o",
                    local_file_trimmed,
                    local_file
                ]
            )
        )
        command.move_file(local_file_trimmed, local_file)

    @staticmethod
    def fasta_union(partial_fasta_files, full_fasta_file):
        ''' Takes a list of fasta file paths and writes the union of the fasta records to full_fasta_file. '''
        if len(partial_fasta_files) == 1:
            command.execute(
                command_patterns.SingleCommand(
                    cmd='ln',
                    args=[
                        "-s",
                        partial_fasta_files[0],
                        full_fasta_file
                    ]
                )
            )
            return
        # For now, just load the files into memory. They are relatively small and
        # the same approach is used in the web app to serve taxon fasta files.
        # TODO: consider changing it to a solution that requires less memory.
        full_fasta = set()
        for fasta_file in partial_fasta_files:
            with open(fasta_file, 'r') as f:
                fasta = set(f.read().split(">"))
            full_fasta.update(fasta)
        full_fasta.discard('')
        with open(full_fasta_file, 'w') as f:
            f.write('>' + '>'.join(full_fasta))
