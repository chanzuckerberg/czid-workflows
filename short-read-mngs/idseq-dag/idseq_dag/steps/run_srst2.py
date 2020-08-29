import os
import pandas as pd
import shutil
from functools import reduce

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.s3 import fetch_reference
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log

MATCHED_READS_FILE = "matched_reads.tsv"

class PipelineStepRunSRST2(PipelineStep):
    '''
    Short Read Sequence Typing for Bacterial Pathogens

    This program is designed to take Illumina sequence data, a MLST database and/or a database
    of gene sequences (e.g. resistance genes, virulence genes, etc) and report the presence of STs and/or reference genes.
    See: https://github.com/katholt/srst2
    '''

    def run(self):
        ''' Invoking srst2 '''
        OUTPUT_LOG = 'output.log'
        OUTPUT_GENES = 'output__genes__ARGannot_r2__results.txt'
        OUTPUT_FULL_GENES = 'output__fullgenes__ARGannot_r2__results.txt'
        is_paired = (len(self.input_files_local[0]) == 2)
        is_fasta = (self.additional_attributes['file_ext'] == 'fasta')
        is_zipped = (self.input_files_local[0][0][-3:] == '.gz')
        self.execute_srst2(is_paired, is_fasta, is_zipped)
        log = os.path.join(self.output_dir_local, OUTPUT_LOG)
        log_dest = self.output_files_local()[0]
        results = os.path.join(self.output_dir_local, OUTPUT_GENES)
        results_dest = self.output_files_local()[1]
        shutil.move(log, log_dest)
        shutil.move(results, results_dest)
        if not os.path.exists(os.path.join(self.output_dir_local, OUTPUT_FULL_GENES)):
            for f in self.output_files_local()[2:6]:
                PipelineStepRunSRST2.fill_file_path(f)
        else:
            # Post processing of amr data
            self.normalize_bam_file()
            self.generate_mapped_reads_tsv()
            total_reads = self.get_total_reads(is_zipped, is_fasta)
            results_full = os.path.join(self.output_dir_local, OUTPUT_FULL_GENES)
            results_full_dest = self.output_files_local()[2]
            shutil.move(results_full, results_full_dest)
            self.process_amr_results(results_full_dest, total_reads)

    # Inherited method
    def count_reads(self):
        pass

    def execute_srst2(self, is_paired, is_fasta, is_zipped):
        """Executes srst2 with appropriate parameters based on whether input files are zipped,
           paired reads and on file type."""
        srst2_params = []
        srst2_params.extend(self.get_common_params())
        if is_fasta:
            file_ext = '.fasta.gz' if is_zipped else '.fasta'
            srst2_params.extend(['--read_type', 'f'])
        else:
            file_ext = '.fastq.gz' if is_zipped else '.fastq'
        if is_paired:
            srst2_params.extend(['--input_pe'])
        else:
            srst2_params.extend(['--input_se'])
        for i, rd in enumerate(self.input_files_local[0]):
            link_name = f"_R{i+1}_001{file_ext}"
            command.execute(
                command_patterns.SingleCommand(
                    cmd='ln',
                    args=[
                        '-sf',
                        rd,
                        link_name
                    ]
                )
            )
            srst2_params.append(link_name)
        if is_paired:
            srst2_params.extend(['--forward', '_R1_001', '--reverse', '_R2_001'])
        command.execute(
            command_patterns.SingleCommand(
                cmd='srst2',
                args=srst2_params
            )
        )

    def get_common_params(self):
        """Helper that gets srst2 parameters common to both paired and single rds."""
        # TODO:  Why is this not fetch_reference?   So it can be cached.
        db_file_path = fetch_reference(self.additional_files["resist_gene_db"], self.ref_dir_local, allow_s3mi=False)  # too small for s3mi
        min_cov = str(self.additional_attributes['min_cov'])
        # srst2 expects this to be a string, in dag could be passed in as a number
        n_threads = str(self.additional_attributes['n_threads'])
        return ['--min_coverage', min_cov, '--threads', n_threads,
                '--output', os.path.join(self.output_dir_local, 'output'), '--log', '--gene_db', db_file_path]

    def normalize_bam_file(self):
        """Ensure files needed are actually present"""
        if os.path.exists(self.output_files_local()[5]):
            return
        # For unpaired fastq inputs, srst2 gives a different name to the sorted bam file that it outputs
        # We rename the bam file to what we expect (as specified in the dag)
        unpaired_bam_path = f'{self.output_dir_local}/output___R1_001.ARGannot_r2.sorted.bam'
        if os.path.exists(unpaired_bam_path):
            command.execute(
                command_patterns.SingleCommand(
                    cmd='mv',
                    args=[
                        unpaired_bam_path,
                        self.output_files_local()[5]
                    ]
                )
            )

    def generate_mapped_reads_tsv(self):
        """Use bedtools to generate a table of mapped reads for each genome in the ARG ANNOT database.
            If a new resistance gene db is used, the .bed file will need to be updated manually."""
        bed_file_path = fetch_reference(self.additional_files["resist_genome_bed"], self.ref_dir_local, allow_s3mi=False)
        sample_bam_file_path = self.output_files_local()[5]

        tmp_sort_dir = os.path.join(self.output_dir_local, "tmp_sort")
        command.make_dirs(tmp_sort_dir)

        # Convert the sorted.bam output from SRST2 to the bed format, then sort the bed file.
        # This allows us to use the "sorted" mode of bedtools coverage, which is memory-efficient.
        # Otherwise, large sorted.bam files will cause our machines to run out of RAM.
        #
        # Note that despite being called "sorted.bam", the bam is not sorted the way we need it to be.
        #
        # env LC_ALL=C ensures that the sort command uses the same sort order on all machines.
        #
        # The -T flag with tmp_sort_dir ensures that we make tmp files inside /mnt, which is where our huge AWS volumes are mounted.
        # By default, the sort command creates temp files in /tmp, which has very little disk space.
        command.execute(
            command_patterns.ShellScriptCommand(
                script='''
                    set -o pipefail;
                    bedtools bamtobed -i "$1" |
                    env LC_ALL=C sort -T "$2" -k1,1 -k2,2n |
                    bedtools coverage -sorted -a "$3" -b stdin > "$4";''',
                args=[
                    sample_bam_file_path,
                    tmp_sort_dir,
                    bed_file_path,
                    os.path.join(self.output_dir_local, MATCHED_READS_FILE)
                ]
            )
        )

        command.remove_rf(tmp_sort_dir)

    def get_total_reads(self, is_zipped, is_fasta):
        """Gets the total number of reads in the sample by counting them directly from the
            fastq or fasta files."""
        # TODO: factor out into utility function, see nonhost_fastq
        input_filenames = self.input_files_local[0]
        if is_zipped:
            unzipped_filenames = []
            for filename in input_filenames:
                if not os.path.exists(filename[:len(filename) - 3]):
                    command.execute(
                        command_patterns.SingleCommand(
                            cmd='gunzip',
                            args=[
                                '-k',
                                filename
                            ]
                        )
                    )
                unzipped_filenames.append(filename[:len(filename) - 3])
            input_filenames = unzipped_filenames
        if is_fasta:  # Number of lines per read can vary, so we use grep
            grep_output = command.execute_with_output(
                command_patterns.SingleCommand(
                    cmd='grep',
                    args=[
                        '-c',
                        '^>',  # fastas start reads with "^>".
                        *input_filenames
                    ]
                )
            )
            output_lines = [line for line in grep_output.split("\n") if line != '']
            if ":" in output_lines[0]:
                # for paired fastas - when run on just one file, grep outputs only
                # a number. But when this command is run on two files, grep outputs
                # a string formatted as filename:count for each file, with count being
                # what we want to add up.
                read_counts = map(lambda line: int(line.split(":")[1]), output_lines)
                return reduce(lambda x, y: x + y, list(read_counts))
            else:
                return int(output_lines[0])
        else:  # fastqs have 4 lines for every read, so we count lines and divide by 4
            wc_params = ['wc', '-l']
            wc_params.extend(input_filenames)
            wc_output = command.execute_with_output(" ".join(wc_params))
            # take the set of characters from the last line, which is the total number of lines
            # for paired reads or the only line for unpaired reads
            wc_lines = [line for line in wc_output.split("\n") if line != '']
            wc_target_line = [line for line in wc_lines[-1].split(" ") if line != '']
            total_line_count = int(wc_target_line[0])
            return total_line_count / 4

    @staticmethod
    def _append_dpm_to_results(amr_results, total_reads):
        """Calculates the depth per million for each gene in the result and appends it to the
            results dataframe."""
        amr_results["dpm"] = amr_results.apply(lambda row: row["depth"] * 1000000 / total_reads, axis=1)
        return amr_results

    @staticmethod
    def _append_rpm_to_results(proc_amr_results, matched_reads_path, total_reads):
        """Reads in a table of matched reads generated by generate_mapped_reads_tsv() and appends
            the reads per million (rpm) and total reads to the dataframe for processed amr results"""
        matched_reads = pd.read_csv(matched_reads_path, delimiter="\t", names=["allele", "reads"], usecols=[0, 3])
        matched_reads["allele"] = matched_reads.apply(lambda row: "_".join(row["allele"].split("__")[2:]), axis=1)
        rpm_list, total_reads_list = PipelineStepRunSRST2._calculate_rpms(matched_reads, proc_amr_results, total_reads)
        proc_amr_results["total_reads"] = total_reads_list
        proc_amr_results["rpm"] = rpm_list
        return proc_amr_results

    @staticmethod
    def _calculate_rpms(rpm_df, amr_df, total_reads):
        """Matches each gene in the amr_results dataframe to it's associated number of reads in
            the matched reads table and calculates the reads per million. Both total reads and
            reads per million are appended in order in lists and then bulk appended back in
            _append_rpm_to_results() to the proc_amr_results dataframe."""
        rpm_list = []
        total_reads_list = []
        for row in amr_df.itertuples():
            if len(rpm_df[rpm_df["allele"] == row.allele]["reads"].values) == 0:
                # This should never happen, but it has happened before
                # because there were typos in allele names in ARGannot_r2.fasta that caused mismatches with argannot_genome.bed.
                # Log an error. The following line will crash the pipeline step, which is intended.
                # We prefer failing the pipeline step to showing incorrect or missing data while failing silently to the user.
                log.write(f"AmrAlleleMismatchError: {row.allele} (from ARGannot_r2.fasta) could not be found in argannot_genome.bed")

            reads_for_allele = rpm_df[rpm_df["allele"] == row.allele]["reads"].values[0]
            total_reads_list.append(reads_for_allele)
            rpm_for_allele = reads_for_allele * 1000000 / total_reads
            rpm_list.append(rpm_for_allele)

        return [rpm_list, total_reads_list]

    @staticmethod
    def fill_file_path(file_path):
        """Helper function to open an "empty" file at a given file location.
           Note that aws s3 cannot upload a 0 byte file from local to s3;
           see  https://github.com/aws/aws-cli/issues/2403.
           Possible problems seem to be from python versions mismatching on aws cli and linux
           system. This doesn't seem to be the case on the staging machine, though.
           So the official recommendation to install aws-cli from pip does not seem to apply since
           the Python versions match up.
           I tried many suggestions on the link above + others -- for now, the following
           seems a reasonable workaround.
           We use os file functions to avoid overhead of using Python's File object
           functions to open an empty file  and using a cp command to move to destination.
           TODO: See if there are aws CLI installation errors that can be fixed."""
        fd = os.open(file_path, os.O_RDWR | os.O_CREAT)
        os.write(fd, b"\n")
        os.close(fd)

    @staticmethod
    def _get_pre_proc_amr_results(amr_raw_path):
        """ Reads in raw amr results file outputted by srst2, and does initial processing of marking gene family."""
        amr_results = pd.read_csv(amr_raw_path, delimiter="\t")
        # Parse out gene family as substring after '_', e.g. Aph_AGly's gene family would be AGly
        amr_results['gene_family'] = amr_results.apply(lambda row: row.gene.split('_', 1)[1], axis=1)
        return amr_results

    @staticmethod
    def _summarize_amr_gene_families(amr_results):
        """Returns a gene family-level summary of total_genes_hit, total_coverage, and total_depth."""
        amr_summary = amr_results.groupby(['gene_family']).agg({'gene_family': ['size'], 'coverage': ['sum'], 'depth': ['sum']})
        amr_summary.columns = [' '.join(col) for col in amr_summary.columns]
        amr_summary = amr_summary.rename(columns={'gene_family size': 'total_gene_hits', 'coverage sum': 'total_coverage', 'depth sum': 'total_depth'}).reset_index()
        return amr_summary

    def process_amr_results(self, amr_results_path, total_reads):
        """ Writes processed amr result table with total_genes_hit, total_coverage,
            and total_depth column values filled in for all genes to output files; and likewise with
            the gene-family level summary of amr results table. """
        amr_results = PipelineStepRunSRST2._get_pre_proc_amr_results(amr_results_path)
        amr_summary = PipelineStepRunSRST2._summarize_amr_gene_families(amr_results)
        amr_summary.to_csv(
            self.output_files_local()[4],
            mode='w',
            index=False,
            encoding='utf-8')
        sorted_amr = amr_results.sort_values(by=['gene_family'])
        proc_amr = pd.merge_ordered(sorted_amr, amr_summary, fill_method='ffill', left_by=['gene_family'])
        proc_amr_with_rpm = PipelineStepRunSRST2._append_rpm_to_results(proc_amr, os.path.join(self.output_dir_local, MATCHED_READS_FILE), total_reads)
        proc_amr_with_rpm_and_dpm = PipelineStepRunSRST2._append_dpm_to_results(proc_amr_with_rpm, total_reads)
        proc_amr_with_rpm_and_dpm.to_csv(
            self.output_files_local()[3],
            mode='w',
            index=False,
            encoding='utf-8')
