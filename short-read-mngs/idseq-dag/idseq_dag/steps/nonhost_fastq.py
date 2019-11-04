import os
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

class PipelineStepNonhostFastq(PipelineStep):
    # Either one or two input read files can be supplied.
    # Works for both FASTA and FASTQ, although non-host FASTQ is more useful.
    def run(self):
        scratch_dir = os.path.join(self.output_dir_local, "scratch_nonhost_fastq")
        command.make_dirs(scratch_dir)
        self.nonhost_headers = [
            os.path.join(scratch_dir, "nonhost_headers_r1.txt"),
            os.path.join(scratch_dir, "nonhost_headers_r2.txt")
        ]

        # Assumed to be [R1.fastq, R2.fastq] if there are two read files.
        fastqs = self.input_files_local[0]

        nonhost_fasta = self.input_files_local[1][0]
        output_fastqs = self.output_files_local()
        fastqs = self.unzip_files(fastqs)

        self.generate_nonhost_headers(nonhost_fasta)

        for i in range(len(fastqs)):
            self.generate_nonhost_fastq(self.nonhost_headers[i], fastqs[i], output_fastqs[i])

        # Clean up scratch files.
        for nonhost_headers in self.nonhost_headers:
            os.remove(nonhost_headers)

    @staticmethod
    # Unzip files with gunzip if necessary.
    def unzip_files(fastqs):
        new_fastqs = []

        for fastq in fastqs:
            if fastq[-3:] == '.gz':
                command.execute(
                    command_patterns.SingleCommand(
                        cmd="gunzip",
                        args=[
                            "-f",
                            "-k",
                            fastq
                        ]
                    )
                )

                new_fastqs.append(fastq[:-3])
            else:
                new_fastqs.append(fastq)

        return new_fastqs

    @staticmethod
    # Extract the original FASTQ header from the fasta file.
    #
    # Example fasta line:
    # >family_nr:4070:family_nt:1903414:genus_nr:4107:genus_nt:586:species_nr:4081:species_nt:587:NR:ABI34274.1:NT:CP029736.1:A00111:123:HCMCTDMXX:1:1111:5575:4382/1
    #
    # We just want the part after NT:XX:
    # We also split based on /1 and /2.
    def extract_header_from_line(line):
        line = line.strip()
        if line[-2:] == "/1":
            read_index = 0
            line = line[:-2]
        elif line[-2:] == "/2":
            read_index = 1
            line = line[:-2]
        else:
            # If there is no suffix /1 or /2, then only one read file was supplied.
            read_index = 0

        fragments = line.split(":")
        nt_index = fragments.index("NT")
        header = ":".join(fragments[nt_index + 2:])

        return {
            "header": header,
            "read_index": read_index
        }

    def generate_nonhost_headers(self, nonhost_fasta_file):
        nonhost_headers = [[], []]
        with open(nonhost_fasta_file, "r") as input_file:
            num = 0
            for line in input_file:
                num += 1
                # Assumes that the header line in the nonhost_fasta starts with ">"
                if line[0] == ">":
                    header = PipelineStepNonhostFastq.extract_header_from_line(line)
                    nonhost_headers[header["read_index"]].append(header["header"] + "\n")

        with open(self.nonhost_headers[0], "w") as output_file:
            output_file.writelines(nonhost_headers[0])

        with open(self.nonhost_headers[1], "w") as output_file:
            output_file.writelines(nonhost_headers[1])

    @staticmethod
    # Use seqtk, which is orders of magnitude faster than Python for this particular step.
    def generate_nonhost_fastq(nonhost_headers, fastq, output_file):
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''seqtk subseq "$1" "$2" > "$3";''',
                args=[
                    fastq,
                    nonhost_headers,
                    output_file
                ]
            )
        )

    def count_reads(self):
        pass
