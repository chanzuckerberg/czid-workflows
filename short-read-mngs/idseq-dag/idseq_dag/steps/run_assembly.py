''' Run Assembly (spades) '''
import json
import os
import traceback
from collections import defaultdict
from Bio import SeqIO
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

from idseq_dag.util.m8 import MIN_CONTIG_SIZE
from idseq_dag.util.count import get_read_cluster_size, load_duplicate_cluster_sizes, READ_COUNTING_MODE, ReadCountingMode


class PipelineStepRunAssembly(PipelineStep):
    """ To obtain longer contigs for improved sensitivity in mapping, short reads must be
    de novo assembled using SPADES.
    The SPADES output loses the information about which contig each individual read belongs to.
    Therefore, we use  bowtie2 to map the original reads onto their assembled contigs.

    1. The short reads are assembled into contigs using SPADES.

    ```
    spades.py
    -1 {input_fasta}
    -2 {input_fasta2}
    -o {assembled_dir}
    -m {memory}
    -t 32
    —only-assembler
    ```

    SPADES documentation can be found [here](http://cab.spbu.ru/software/spades/)

    2. The single-read identity of reads merged into each contig are lost by SPADES.
    To recover this information and identify which contig each read belongs to,
    the contigs are then used to build a Bowtie2 database:

    ```
    bowtie2-build {assembled_contig} {bowtie_index_path}
    ```

    3. The original reads are mapped back to their assembled contigs:

    ```
    bowtie2
    -x {bowtie_index_path}
    -f
    -U {fasta_file}
    --very-sensitive
    -p 32 > {output_bowtie_sam}
    ```
    """

    def run(self):
        """
           Run Assembly
        """
        input_fasta = self.input_files_local[0][-1]
        bowtie_fasta = self.input_files_local[0][-1]
        input_fasta2 = None
        if len(self.input_files_local[0]) >= 2:
            input_fasta = self.input_files_local[0][0]
            input_fasta2 = self.input_files_local[0][1]
        duplicate_cluster_sizes_path, = self.input_files_local[1]
        assert duplicate_cluster_sizes_path.endswith(".tsv"), self.input_files_local[1]

        assembled_contig, assembled_contig_all, assembled_scaffold, bowtie_sam, contig_stats = self.output_files_local()
        read2contig = {}
        memory = self.additional_attributes.get('memory', 100)
        min_contig_length = int(self.additional_attributes.get('min_contig_length', 0))
        self.assemble(input_fasta, input_fasta2, bowtie_fasta, duplicate_cluster_sizes_path,
                      assembled_contig, assembled_contig_all, assembled_scaffold,
                      bowtie_sam, contig_stats, read2contig, int(memory), min_contig_length)

    @staticmethod
    def assemble(input_fasta,
                 input_fasta2,
                 bowtie_fasta,  # fasta file for running bowtie against contigs
                 duplicate_cluster_sizes_path,
                 assembled_contig,
                 assembled_contig_all,
                 assembled_scaffold,
                 bowtie_sam,
                 contig_stats,
                 read2contig,
                 memory=100,
                 min_contig_length=0):
        basedir = os.path.dirname(assembled_contig)
        assembled_dir = os.path.join(basedir, 'spades')
        command.make_dirs(assembled_dir)
        assembled_contig_tmp = os.path.join(assembled_dir, 'contigs.fasta')
        assembled_scaffold_tmp = os.path.join(assembled_dir, 'scaffolds.fasta')

        try:
            if input_fasta2:
                command.execute(
                    command_patterns.SingleCommand(
                        cmd="spades.py",
                        args=[
                            "-1",
                            input_fasta,
                            "-2",
                            input_fasta2,
                            "-o",
                            assembled_dir,
                            "-m",
                            memory,
                            "-t",
                            32,
                            "--only-assembler"
                        ]
                    )
                )
            else:
                command.execute(
                    command_patterns.SingleCommand(
                        cmd="spades.py",
                        args=[
                            "-s",
                            input_fasta,
                            "-o",
                            assembled_dir,
                            "-m",
                            memory,
                            "-t",
                            32,
                            "--only-assembler"
                        ]
                    )
                )
            command.move_file(assembled_contig_tmp, assembled_contig_all)
            command.move_file(assembled_scaffold_tmp, assembled_scaffold)

            if min_contig_length:
                # apply contig length filter
                SeqIO.write(
                    (r for r in SeqIO.parse(assembled_contig_all, "fasta")
                        if len(r.seq) >= min_contig_length),
                    assembled_contig,
                    "fasta",
                )
            else:
                command.copy_file(assembled_contig_all, assembled_contig)

            PipelineStepRunAssembly.generate_read_to_contig_mapping(assembled_contig, bowtie_fasta,
                                                                    read2contig, duplicate_cluster_sizes_path, bowtie_sam, contig_stats)
        except:
            # Assembly failed. create dummy output files
            command.write_text_to_file(';ASSEMBLY FAILED', assembled_contig)
            command.write_text_to_file(';ASSEMBLY FAILED', assembled_contig_all)
            command.write_text_to_file(';ASSEMBLY FAILED', assembled_scaffold)
            command.write_text_to_file('@NO INFO', bowtie_sam)
            command.write_text_to_file('{}', contig_stats)
            traceback.print_exc()
        command.remove_rf(assembled_dir)

    @staticmethod
    def generate_read_to_contig_mapping(assembled_contig,
                                        fasta_file,
                                        read2contig,
                                        duplicate_cluster_sizes_path,
                                        output_bowtie_sam,
                                        output_contig_stats):
        ''' read -> contig mapping through bowtie2 alignment '''
        base_output_dir = os.path.dirname(fasta_file)
        # build bowtie index based on assembled_contig
        bowtie_index_path = os.path.join(base_output_dir, 'bowtie-contig')
        command.make_dirs(bowtie_index_path)
        command.execute(
            command_patterns.SingleCommand(
                cmd='bowtie2-build',
                args=[
                    assembled_contig,
                    bowtie_index_path
                ]
            )
        )
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''bowtie2 -x "${bowtie_index_path}" -f -U "${fasta_file}" --very-sensitive -p 32 > "${output_bowtie_sam}";''',
                named_args={
                    'bowtie_index_path': bowtie_index_path,
                    'fasta_file': fasta_file,
                    'output_bowtie_sam': output_bowtie_sam
                }
            )
        )
        contig_stats = PipelineStepRunAssembly.generate_info_from_sam(output_bowtie_sam, read2contig, duplicate_cluster_sizes_path)
        with open(output_contig_stats, 'w') as ocf:
            json.dump(contig_stats, ocf)

    @staticmethod
    def generate_info_from_sam(bowtie_sam_file, read2contig, duplicate_cluster_sizes_path):
        contig_stats = defaultdict(int)
        contig_unique_counts = defaultdict(int)
        duplicate_cluster_sizes = load_duplicate_cluster_sizes(duplicate_cluster_sizes_path)
        with open(bowtie_sam_file, "r", encoding='utf-8') as samf:
            for line in samf:
                if line[0] == '@':
                    continue
                fields = line.split("\t")
                read = fields[0]
                contig = fields[2]
                contig_stats[contig] += get_read_cluster_size(duplicate_cluster_sizes, read)  # these are non-unique read counts now
                contig_unique_counts[contig] += 1
                if contig != '*':
                    read2contig[read] = contig
        for contig, unique_count in contig_unique_counts.items():  # TODO can't we just filter those out after spades, IN ONE PLACE
            if unique_count < MIN_CONTIG_SIZE:
                del contig_stats[contig]
            elif READ_COUNTING_MODE == ReadCountingMode.COUNT_UNIQUE:
                contig_stats[contig] = unique_count
        return contig_stats

    def count_reads(self):
        ''' count reads '''
        pass
