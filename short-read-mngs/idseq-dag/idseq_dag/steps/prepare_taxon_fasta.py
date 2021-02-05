''' Generate Phylogenetic tree '''
import os
import json

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
        os.mkdir("taxon_fastas")
        json_input_path = self.input_files_local[0][0]
        with open(json_input_path) as f:
            json_input = json.load(f)
        for sample in json_input:
            output_file = os.path.join("taxon_fastas", str(sample["run_id"]) + ".fasta")
            # Retrieve IDseq taxon fasta files
            partial_fasta_files = []
            for hit_type in ["nr", "nt"]:
                byterange = sample[f"taxon_byterange_{hit_type}"]
                first_byte = byterange["first_byte"]
                last_byte = byterange["last_byte"]
                s3_file = byterange["refined_taxid_annotated_sorted_fasta_s3_path"]
                local_file = f"{hit_type}_{os.path.basename(output_file)}"
                bucket, key = s3.split_identifiers(s3_file)
                s3.fetch_byterange(first_byte, last_byte, bucket, key, local_file)
                partial_fasta_files.append(local_file)
            self.fasta_union(partial_fasta_files, output_file)
            for fasta in partial_fasta_files + [output_file]:
                print(f"{count.reads(fasta)} unique reads in {fasta}")

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
