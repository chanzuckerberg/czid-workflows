import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

def reads_in_group(file_group, max_fragments=None):
    '''
    Count reads in a group of matching files, up to a maximum number of fragments.
    The term "fragment" refers to the physical DNA fragments the reads derive from.
    If the input is single, then 1 fragment == 1 read.
    If the input is paired, then 1 fragment == 2 reads.
    '''
    num_files = len(file_group)
    reads_in_first_file = reads(file_group[0], max_fragments)
    return num_files * reads_in_first_file

def reads(local_file_path, max_reads=None):
    '''
    Count reads in a local file based on file format inferred from extension,
    up to a maximum of max_reads.
    '''
    if local_file_path.endswith(".gz"):
        cmd = r'''zcat "${local_file_path}"'''
        file_format = local_file_path.split(".")[-2]
    else:
        cmd = r'''cat "${local_file_path}"'''
        file_format = local_file_path.split(".")[-1]

    named_args = {
        'local_file_path': local_file_path
    }

    if max_reads:
        max_lines = reads2lines(max_reads, file_format)
        assert max_lines is not None, "Could not convert max_reads to max_lines"
        cmd += r''' | head -n "${max_lines}"'''
        named_args.update({
            'max_lines': max_lines
        })

    cmd += " |  wc -l"

    cmd_output = command.execute_with_output(
        command_patterns.ShellScriptCommand(
            script=cmd,
            named_args=named_args
        )
    )
    line_count = int(cmd_output.strip().split(' ')[0])
    return lines2reads(line_count, file_format)

def lines2reads(line_count, file_format):
    '''
    Convert line count to read count based on file format.
    Supports fastq and SINGLE-LINE fasta formats.
    TODO: add support for m8 files here once the relevant steps are added in the pipeline engine.
    '''
    read_count = None
    if file_format in ["fq", "fastq"]:
        # Each read consists of exactly 4 lines, by definition.
        assert line_count % 4 == 0, "File does not follow fastq format"
        read_count = line_count // 4
    elif file_format in ["fa", "fasta"]:
        # Each read consists of exactly 2 lines (identifier and sequence).
        # ASSUMES the format is SINGLE-LINE fasta files, where each read's sequence
        # is on a single line rather than being spread over multiple lines.
        # This is fine for now, because we only count fasta files we generated
        # ourselves, following the single-line format.
        assert line_count % 2 == 0, "File does not follow single-line fasta format"
        read_count = line_count // 2
    else: 
        read_count = line_count
    return read_count

def reads2lines(read_count, file_format):
    '''
    Convert read count to line count based on file format.
    Currently supports fastq or SINGLE-LINE fasta.
    '''
    line_count = None
    if file_format in ["fq", "fastq"]:
        line_count = 4 * read_count
    elif file_format in ["fa", "fasta"]:
        # Assumes the format is single-line fasta rather than multiline fasta
        line_count = 2 * read_count
    return line_count

def files_have_min_reads(input_files, min_reads):
    """ Checks whether fa/fasta/fq/fastq have the minimum number of reads.
    Pipeline steps can use this method for input validation.
    """
    for input_file in input_files:
        num_reads = reads(input_file, max_reads=min_reads)
        if num_reads < min_reads:
            return False
    return True

