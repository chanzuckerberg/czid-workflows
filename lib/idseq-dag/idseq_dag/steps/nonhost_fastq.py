import os

from typing import Dict, Optional, Sequence, Set, List, Tuple

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.czid_dedup_clusters import parse_clusters_file
from idseq_dag.util.count import READ_COUNTING_MODE, ReadCountingMode


class PipelineStepNonhostFastq(PipelineStep):
    # Either one or two input read files can be supplied.
    # Works for both FASTA and FASTQ, although non-host FASTQ is more useful.
    def run(self) -> None:
        clusters_dict = None
        if READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL:
            # NOTE: this will load the set of all original read headers, which
            # could be several GBs in the worst case.
            clusters_dict = parse_clusters_file(self.input_files_local[2][0])

        self.run_nonhost_fastq_generation(clusters_dict)

    def run_nonhost_fastq_generation(
        self,
        clusters_dict: Dict[str, List] = None,
    ) -> None:
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

        self.generate_nonhost_headers(nonhost_fasta, clusters_dict)

        for i in range(len(fastqs)):
            self.generate_nonhost_fastq(self.nonhost_headers[i], fastqs[i], output_fastqs[i])

        # Clean up scratch files.
        for nonhost_headers in self.nonhost_headers:
            os.remove(nonhost_headers)

    @staticmethod
    # Unzip files with gunzip if necessary.
    def unzip_files(fastqs: Sequence[str]) -> Sequence[str]:
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
    def extract_header_from_line(line: str) -> Tuple[int, str, Set[int]]:
        """
        Extract the original FASTQ header from the fasta file.

        Example fasta line:
        >family_nr:4070:family_nt:1903414:genus_nr:4107:genus_nt:586:species_nr:4081:species_nt:587:NR:ABI34274.1:NT:CP029736.1:A00111:123:HCMCTDMXX:1:1111:5575:4382/1

        We want to extract the taxids added upstream then we want only the part
        after NT:XX as the read ID. We also split R1 and R2 reads based on /1
        and /2.
        """
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

        return read_index, header

    def generate_nonhost_headers(
        self,
        nonhost_fasta_file: str,
        clusters_dict: Dict[str, List] = None,
    ):
        with open(nonhost_fasta_file, "r") as input_file, \
                open(self.nonhost_headers[0], "w") as output_file_0, \
                open(self.nonhost_headers[1], "w") as output_file_1:
            for line in input_file:
                # Assumes that the header line in the nonhost_fasta starts with ">"
                if line[0] != ">":
                    continue
                read_index, header = PipelineStepNonhostFastq.extract_header_from_line(line)
                if clusters_dict:
                    if header not in clusters_dict:
                        header += "/2" if read_index else "/1"
                    other_headers = clusters_dict[header][1:]
                else:
                    other_headers = []
                other_headers = clusters_dict[header][1:] if clusters_dict else []

                output_file = output_file_0 if read_index == 0 else output_file_1
                output_file.write(header + "\n")
                for other_header in other_headers:
                    output_file.write(other_header + "\n")

    @staticmethod
    # Use seqtk, which is orders of magnitude faster than Python for this particular step.
    def generate_nonhost_fastq(
        nonhost_headers: str,
        fastq: str,
        output_file: str
    ) -> None:
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
