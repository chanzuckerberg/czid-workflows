import json
import os

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.exceptions import InvalidFileFormatError, InsufficientReadsError
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.validate_constants as vc
import idseq_dag.util.s3 as s3

class PipelineStepRunValidateInput(PipelineStep):
    """ Validates that the input files are .fastq format and truncates to 75 million fragments
    (specifically, 75 million reads for single-end libraries or 150 million reads for paired-end libraries).
    The validation process counts the number of sequences that fall in specified length buckets, which will
    inform the parameters used downstream for initial host removal by STAR.
    """
    def run(self):
        # Setup
        input_files = self.input_files_local[0][0:2]
        num_inputs = len(input_files)
        assert num_inputs in [1, 2], 'Invalid number of input files'
        output_files = self.output_files_local()[1:3]
        summary_file = self.output_files_local()[0]
        max_fragments = self.additional_attributes["truncate_fragments_to"]

        file_ext = self.additional_attributes.get("file_ext")
        assert file_ext in ['fastq', 'fasta'], 'Invalid file extension'

        is_fastq = file_ext == 'fastq'

        try:
            for i in range(num_inputs):
                input_file = input_files[i]
                splited_input_file_name, splited_input_file_ext = os.path.splitext(input_file)

                num_lines = self.calc_max_num_lines(is_fastq, max_fragments)

                # unzip if .gz file
                if splited_input_file_ext == '.gz':
                    input_files[i] = splited_input_file_name
                    try:
                        # test if a valid gzip file
                        command.execute(
                            command_patterns.SingleCommand(
                                cmd="gzip",
                                args=[
                                    "-t",
                                    input_file
                                ]
                            )
                        )
                        # then decompress it
                        command.execute(
                            command_patterns.ShellScriptCommand(
                                script=r'''gzip -dc "${input_file}" | cut -c -"$[max_line_length+1]" | head -n "${num_lines}" | awk -f "${awk_script_file}" -v max_line_length="${max_line_length}" > "${output_file}";''',
                                named_args={
                                    "input_file": input_file,
                                    "awk_script_file": command.get_resource_filename("scripts/fastq-fasta-line-validation.awk"),
                                    "max_line_length": vc.MAX_LINE_LENGTH,
                                    "num_lines": num_lines,
                                    "output_file": splited_input_file_name
                                }
                            )
                        )
                    except:
                        raise InvalidFileFormatError("Invalid fastq/fasta/gzip file")
                else:
                    # Validate and truncate the input file to keep behavior consistent with gz input files
                    try:
                        tmp_file = splited_input_file_name + ".tmp"
                        command.execute(
                            command_patterns.ShellScriptCommand(
                                script=r'''cat "${input_file}" | cut -c -"$[max_line_length+1]" | head -n "${num_lines}" | awk -f "${awk_script_file}" -v max_line_length="${max_line_length}" > "${output_file}";''',
                                named_args={
                                    "input_file": input_file,
                                    "awk_script_file": command.get_resource_filename("scripts/fastq-fasta-line-validation.awk"),
                                    "max_line_length": vc.MAX_LINE_LENGTH,
                                    "num_lines": num_lines,
                                    "output_file": tmp_file
                                }
                            )
                        )
                        input_files[i] = tmp_file
                    except:
                        raise InvalidFileFormatError("Invalid fastq/fasta file")

            # keep a dictionary of the distribution of read lengths in the files
            self.summary_dict = {vc.BUCKET_TOO_SHORT: 0,
                                 vc.BUCKET_NORMAL: 0,
                                 vc.BUCKET_LONG: 0,
                                 vc.BUCKET_TOO_LONG: 0}
            # add a total reads output as source of truth (if filtering changes)
            self.total_output_reads = 0

            quick_check_passed = \
                self.quick_check_file(input_files[0], is_fastq) and \
                (num_inputs == 1 or self.quick_check_file(input_files[1], is_fastq))

            all_fragments = []

            for infile, outfile in zip(input_files, output_files):
                if quick_check_passed:
                    num_fragments = self.truncate_file(infile, outfile, is_fastq, max_fragments)
                else:
                    num_fragments = self._full_check_and_truncate_file(infile, outfile, is_fastq, max_fragments, num_inputs)
                all_fragments.append(num_fragments)

            if len(all_fragments) == 2 and abs(all_fragments[1] - all_fragments[0]) > 1000:
                raise InvalidFileFormatError("Paired input files need to contain the same number of reads")

            with open(summary_file, 'w') as summary_f:
                json.dump(self.summary_dict, summary_f)

        except Exception as e:
            with open(summary_file, 'w') as summary_f:
                json.dump({'Validation error': str(e)}, summary_f)
            s3_path = self.s3_path(summary_file)
            s3.upload_with_retries(summary_file, s3_path)
            raise e

        return

    # quick_check_file returns:
    #   True if the first 100 fragments all have the same length of reads and
    #   are well-formated single-line FASTA / FASTQ entries.
    #
    #   False if the entries are not formatted simply or the read lengths are not
    #   all identical or if there is another possibly recoverable error
    #
    #   Throws an exception in the case of an unrecoverable abnormality
    def quick_check_file(self, file, is_fastq, max_fragments_to_check=100):
        num_fragments = 0
        fragment_length = 0

        with open(file, 'r', encoding='utf-8') as input_f:
            while True:
                num_fragments += 1
                if num_fragments > max_fragments_to_check:
                    break

                identifier_l = input_f.readline()
                if len(identifier_l) == 0:  # EOF
                    if num_fragments == 1:
                        raise InsufficientReadsError("The input file contains 0 reads")
                    break

                read_l = input_f.readline()
                if len(read_l) == 0:  # unexpected EOF
                    raise InvalidFileFormatError("Invalid input file")

                if is_fastq:
                    identifier2_l = input_f.readline()
                    if len(identifier2_l) == 0:
                        raise InvalidFileFormatError("Invalid FASTQ file")

                    quality_l = input_f.readline()
                    if len(quality_l) == 0:
                        raise InvalidFileFormatError("Invalid FASTQ file")

                if is_fastq:
                    if identifier_l[0] != '@' or identifier2_l[0] != '+':
                        # may be FASTQ file with multi-line reads, requires full check
                        return False
                else:
                    if identifier_l[0] != '>':
                        # may be FASTA file with multi-line reads, requires full check
                        return False

                if fragment_length == 0:
                    fragment_length = len(read_l)
                    if fragment_length < vc.READ_LEN_CUTOFF_LOW or fragment_length > vc.READ_LEN_CUTOFF_MID:
                        # non-standard fragment lengths require more detailed examination
                        return False

                if fragment_length != len(read_l) or (is_fastq and fragment_length != len(quality_l)):
                    # file does not meet "quick check" requirements since fragments/quality
                    # scores are not all of same length
                    return False

        return True

    def calc_max_num_lines(self, is_fastq, max_fragments):
        if is_fastq:
            num_lines = max_fragments * 4
        else:
            num_lines = max_fragments * 2
        return num_lines

    def truncate_file(self, infile, outfile, is_fastq, max_fragments):
        num_lines = self.calc_max_num_lines(is_fastq, max_fragments)
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''head -n "${num_lines}" "${infile}" > "${outfile}";''',
                named_args={
                    'num_lines': num_lines,
                    'infile': infile,
                    'outfile': outfile
                }
            )
        )
        num_fragments = count.reads(outfile)
        self.summary_dict[vc.BUCKET_NORMAL] += num_fragments
        return num_fragments

    # _full_check_and_truncate_file does an exhaustive check of the input file, up to
    # max_fragments, and reformats the output to conform to what the rest of the
    # computational pipeline expects (single-line reads of max-length 10,000). After
    # viewing max_fragments reads or encountering EOF, the function returns.
    #
    # Throws an exception in the case of an unrecoverable abnormality
    def _full_check_and_truncate_file(self, infile, outfile, is_fastq, max_fragments, num_inputs):
        num_fragments = 0

        with open(infile, 'r', encoding='utf-8') as input_f, open(outfile, 'w') as output_f:
            next_line = input_f.readline()
            while True:
                num_fragments += 1
                if num_fragments > max_fragments:
                    break

                identifier_l = next_line
                if len(identifier_l) == 0:  # EOF
                    break

                read_l = input_f.readline()
                if len(read_l) == 0:
                    raise InvalidFileFormatError("Invalid input file")

                read_l = read_l.rstrip()
                next_line = input_f.readline()
                while len(next_line) > 0 and next_line[0] not in ['>', '@', '+']:
                    read_l += next_line.rstrip()
                    next_line = input_f.readline()

                if is_fastq:
                    identifier2_l = next_line
                    if len(identifier2_l) == 0:
                        raise InvalidFileFormatError("Invalid FASTQ file")

                    quality_l = input_f.readline()
                    if len(quality_l) == 0:
                        raise InvalidFileFormatError("Invalid FASTQ file")

                    quality_l = quality_l.rstrip()
                    next_line = input_f.readline()
                    while len(next_line) > 0 and next_line[0] not in ['>', '@', '+']:
                        quality_l += next_line.rstrip()
                        next_line = input_f.readline()

                if is_fastq:
                    if identifier_l[0] != '@':
                        raise InvalidFileFormatError("Invalid FASTQ file")
                    if identifier2_l[0] != '+':
                        raise InvalidFileFormatError("Invalid FASTQ file")
                else:
                    if identifier_l[0] != '>':
                        raise InvalidFileFormatError("Invalid FASTA file")

                # At this point, identifier_l and identifier2_l end in a newline and
                # read_l and quality_l do not end in a newline
                read_len = len(read_l)

                # Force read and quality lengths to be identical, either by padding quality
                # with the last quality score or truncating quality score
                if is_fastq:
                    if read_len > len(quality_l):
                        quality_l += (quality_l[-1] * (read_len - len(quality_l)))
                    elif read_len < len(quality_l):
                        quality_l = quality_l[0:read_len]

                if read_len < vc.READ_LEN_CUTOFF_LOW:
                    self.summary_dict[vc.BUCKET_TOO_SHORT] += 1
                elif read_len < vc.READ_LEN_CUTOFF_MID:
                    self.summary_dict[vc.BUCKET_NORMAL] += 1
                elif read_len < vc.READ_LEN_CUTOFF_HIGH:
                    self.summary_dict[vc.BUCKET_LONG] += 1
                else:
                    self.summary_dict[vc.BUCKET_TOO_LONG] += 1
                    read_l = read_l[0:vc.READ_LEN_CUTOFF_HIGH]
                    if is_fastq:
                        quality_l = quality_l[0:vc.READ_LEN_CUTOFF_HIGH]

                self.total_output_reads += 1
                output_f.write(identifier_l + read_l + "\n")
                if is_fastq:
                    output_f.write(identifier2_l + quality_l + "\n")

        return num_fragments

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = self.total_output_reads
