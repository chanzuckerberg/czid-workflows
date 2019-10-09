''' Generate Host Genome given a fasta'''
import os
import multiprocessing
from psutil import virtual_memory

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3

class PipelineStepGenerateHostGenome(PipelineStep):
    ''' Generate Host Genome PipelineStep implementation '''
    def run(self):
        """
        Generate host genome indexes for STAR and bowtie2
        """
        # Set up
        input_fasta_path = self.input_files_local[0][0]
        ercc_fasta_path = s3.fetch_from_s3(self.additional_files["ercc_fasta"],
                                           self.output_dir_local,
                                           allow_s3mi=True,
                                           auto_unzip=True)
        if input_fasta_path[-3:] == '.gz':
            # unzip the file
            dest_path = input_fasta_path[:-3]
            command.execute(
                command_patterns.ShellScriptCommand(
                    script=r'''gzip -dc "${input_fasta_path}" > "${dest_path}";''',
                    named_args={
                        'input_fasta_path': input_fasta_path,
                        'dest_path': dest_path
                    }
                )
            )

            input_fasta_path = dest_path

        input_gtf_path = None
        if self.additional_files.get("input_gtf"):
            input_gtf_path = s3.fetch_from_s3(self.additional_files["input_gtf"],
                                              self.output_dir_local,
                                              allow_s3mi=True)

        ercc_gtf_path = s3.fetch_from_s3(self.additional_files["ercc_gtf"],
                                         self.output_dir_local,
                                         allow_s3mi=True,
                                         auto_unzip=True)

        host_name = self.additional_attributes["host_name"]
        max_star_part_size = self.additional_attributes.get("max_star_part_size")
        input_fasta_with_ercc = f"{input_fasta_path}.with_ercc"
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''cat "${ercc_fasta_path}" "${input_fasta_path}" > "${input_fasta_with_ercc}";''',
                named_args={
                    'ercc_fasta_path': ercc_fasta_path,
                    'input_fasta_path': input_fasta_path,
                    'input_fasta_with_ercc': input_fasta_with_ercc
                }
            )
        )

        input_gtf_with_ercc = ercc_gtf_path
        if input_gtf_path:
            input_gtf_with_ercc = f"{input_gtf_path}.with_ercc"
            command.execute(
                command_patterns.ShellScriptCommand(
                    script=r'''cat "${ercc_gtf_path}" "${input_gtf_path}" > "${input_gtf_with_ercc}";''',
                    named_args={
                        'ercc_gtf_path': ercc_gtf_path,
                        'input_gtf_path': input_gtf_path,
                        'input_gtf_with_ercc': input_gtf_with_ercc
                    }
                )
            )

        output_fasta_file, output_gtf_file, output_star_index, output_bowtie2_index = self.output_files_local()

        command.copy_file(input_fasta_with_ercc, output_fasta_file)
        command.copy_file(input_gtf_with_ercc, output_gtf_file)

        # make STAR index
        self.make_star_index(
            input_fasta_with_ercc, input_gtf_with_ercc, output_star_index, max_star_part_size)

        # make bowtie2 index
        self.make_bowtie2_index(host_name, input_fasta_with_ercc, output_bowtie2_index)

    @staticmethod
    def split_fasta(fasta_file, max_fasta_part_size):
        fasta_file_list = []
        part_idx = 0
        current_size = 0
        current_output_file_name = "%s.%d" % (fasta_file, part_idx)
        current_output_file = open(current_output_file_name, 'wb')
        fasta_file_list.append(current_output_file_name)
        with open(fasta_file, 'rb') as input_f:
            current_read = input_f.readline()
            for line in input_f:
                # Check if we have to switch different output fasta file
                if current_size > max_fasta_part_size:
                    current_output_file.close()
                    part_idx += 1
                    current_size = 0
                    current_output_file_name = "%s.%d" % (fasta_file, part_idx)
                    current_output_file = open(current_output_file_name, 'wb')
                    fasta_file_list.append(current_output_file_name)

                if line[0] == '>':  # got a new read
                    current_output_file.write(current_read)
                    current_size += len(current_read)
                    current_read = line
                else:
                    current_read += line
            current_output_file.write(current_read)
            current_output_file.close()
        return fasta_file_list

    @staticmethod
    def make_star_index(fasta_file, gtf_file, output_star_genome_path, max_star_part_size):
        star_genome_dir_name = output_star_genome_path[:-4]

        # star genome organization
        # STAR_genome/part-${i}, parts.txt
        fasta_file_list = []
        if max_star_part_size and os.path.getsize(fasta_file) > max_star_part_size:
            fasta_file_list = PipelineStepGenerateHostGenome.split_fasta(
                fasta_file, max_star_part_size)
        else:
            fasta_file_list.append(fasta_file)

        for i in range(len(fasta_file_list)):
            log.write("start making STAR index part %d" % i)
            gtf_command_part = []
            if i == 0 and gtf_file:
                gtf_command_part = ["--sjdbGTFfile", gtf_file]

            star_genome_part_dir = f"{star_genome_dir_name}/part-{i}"

            command.make_dirs(star_genome_part_dir)
            star_command_params = [
                '--runThreadN',
                str(multiprocessing.cpu_count()), '--runMode', 'genomeGenerate',
                *gtf_command_part, '--genomeDir', star_genome_part_dir,
                '--genomeFastaFiles', fasta_file_list[i],
                '--limitGenomeGenerateRAM', virtual_memory().available
            ]
            command.execute(
                command_patterns.SingleCommand(
                    cmd='STAR',
                    args=star_command_params
                )
            )
            log.write(f"finished making STAR index part {i}")
        # record # parts into parts.txt
        command.write_text_to_file(len(fasta_file_list), os.path.join(star_genome_dir_name, "parts.txt"))
        star_genome = os.path.basename(star_genome_dir_name)
        star_work_dir = os.path.dirname(star_genome_dir_name)
        command.execute(
            command_patterns.SingleCommand(
                cmd="tar",
                args=[
                    "cvf",
                    output_star_genome_path,
                    "-C",
                    star_work_dir,
                    star_genome
                ]
            )
        )

    @staticmethod
    def make_bowtie2_index(host_name, fasta_file, output_bowtie2_index):
        bowtie2_genome_dir_name = output_bowtie2_index[:-4]
        command.make_dirs(bowtie2_genome_dir_name)
        command.execute(
            command_patterns.SingleCommand(
                cd=bowtie2_genome_dir_name,
                cmd='bowtie2-build',
                args=[
                    fasta_file,
                    host_name
                ]
            )
        )
        log.write("finished making bowtie2 index")
        # archive
        bowtie_genome = os.path.basename(bowtie2_genome_dir_name)
        bowtie_work_dir = os.path.dirname(bowtie2_genome_dir_name)
        command.execute(
            command_patterns.SingleCommand(
                cmd="tar",
                args=[
                    "cvf",
                    output_bowtie2_index,
                    "-C",
                    bowtie_work_dir,
                    bowtie_genome
                ]
            )
        )


    def count_reads(self):
        ''' Count reads '''
        pass
