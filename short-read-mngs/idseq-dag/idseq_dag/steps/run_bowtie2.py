import os
import multiprocessing
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.convert as convert
import idseq_dag.util.log as log
import idseq_dag.util.count as count
from idseq_dag.util.s3 import fetch_from_s3

class PipelineStepRunBowtie2(PipelineStep):
    def run(self):
        input_fas = self.input_files_local[0][0:2]
        output_fas = self.output_files_local()
        genome_dir = fetch_from_s3(self.additional_files["bowtie2_genome"],
                                   self.ref_dir_local,
                                   allow_s3mi=True,
                                   auto_untar=True)
        output_sam_file = os.path.join(self.output_dir_local,
                                       self.additional_attributes["output_sam_file"])
        self.additional_files_to_upload.append(output_sam_file)
        # The file structure looks like
        # "bowtie2_genome/GRCh38.primary_assembly.genome.3.bt2"
        # The code below will handle up to "bowtie2_genome/GRCh38.primary_assembly.
        # genome.99.bt2" but not 100.
        cmd = "ls {genome_dir}/*.bt2*".format(genome_dir=genome_dir)
        local_genome_dir_ls = command.execute_with_output(cmd).decode('utf-8')
        genome_basename = local_genome_dir_ls.split("\n")[0][:-6]
        if genome_basename[-1] == '.':
            genome_basename = genome_basename[:-1]
        bowtie2_params = [
            'bowtie2', '-q', '-p',
            str(multiprocessing.cpu_count()), '-x', genome_basename, '-f',
            '--very-sensitive-local', '-S', output_sam_file
        ]

        if len(input_fas) == 2:
            bowtie2_params.extend(['-1', input_fas[0], '-2', input_fas[1]])
        else:
            bowtie2_params.extend(['-U', input_fas[0]])
        command.execute(" ".join(bowtie2_params))
        log.write("Finished Bowtie alignment.")

        if len(input_fas) == 2:
            convert.generate_unmapped_pairs_from_sam(
                output_sam_file, output_fas)
        else:
            convert.generate_unmapped_singles_from_sam(
                output_sam_file, output_fas[0])

    def count_reads(self):
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])
