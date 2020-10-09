import gzip
import multiprocessing
from enum import Enum
from subprocess import run, PIPE

import idseq_dag.util.fasta as fasta
from idseq_dag import __version__
from idseq_dag.exceptions import InvalidFileFormatError

class ReadCountingMode(Enum):
    COUNT_UNIQUE = "COUNT UNIQUE READS"
    COUNT_ALL = "COUNT ALL READS"


READ_COUNTING_MODE = ReadCountingMode.COUNT_ALL

GZIP_MAGIC_HEADER = b'\037\213'

def count_reads(filename):
    '''
    Count reads in a given FASTA or FASTQ file.
    '''
    with open(filename, "rb") as gz_fh:
        is_gzipped = True if gz_fh.read(2).startswith(GZIP_MAGIC_HEADER) else False
    with gzip.open(filename) if is_gzipped else open(filename, mode="rb") as fmt_fh:
        chunk = fmt_fh.read(1)
        if len(chunk) == 0:
            return 0
        first_char = chunk.decode()[0]
    with open(filename, "rb") as fh:
        if first_char == ">":
            cmd = "grep -c '^>'"
            if is_gzipped:
                cmd = "gunzip | " + cmd
            return int(run(cmd, stdin=fh, stdout=PIPE, check=True, shell=True).stdout)
        elif first_char == "@":
            cmd = "wc -l"
            if is_gzipped:
                cmd = "gunzip | " + cmd
            num_lines = int(run(cmd, stdin=fh, stdout=PIPE, check=True, shell=True).stdout)
            if num_lines % 4 != 0:
                raise InvalidFileFormatError("File does not follow fastq format")
            return num_lines // 4
        raise InvalidFileFormatError("Unable to recognize file format")


reads = count_reads


def files_have_min_reads(input_files, min_reads):
    """ Checks whether fa/fasta/fq/fastq have the minimum number of reads.
    Pipeline steps can use this method for input validation.
    """
    for input_file in input_files:
        num_reads = count_reads(input_file)
        if num_reads < min_reads:
            return False
    return True


def _count_reads_expanding_duplicates(local_file_path, cluster_sizes, cluster_key):
    # See documentation for reads_in_group use case with cluster_sizes, below.
    unique_count, nonunique_count = 0, 0
    for read in fasta.iterator(local_file_path):
        # A read header looks someting like
        #
        #    >M05295:357:000000000-CRPNR:1:1101:22051:10534 OPTIONAL RANDOM STUFF"
        #     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #
        # where the first character on the line is '>' and the read ID (underlined above)
        # extends from '>' to the first whitespace character, not including '>' itself.
        #
        # The fasta iterator already asserts that read.header[0] is '>'.
        #
        # As we proceed down along the pipeline, read IDs get annotated with taxonomic information,
        # changing the above into something like
        #
        #   >NT:ABC2433.1:NR:ABC5656.2:M05295:357:000000000-CRPNR:1:1101:22051:10534 OPTIONAL RANDOM STUFF"
        #    ^^^^^^^^^^^^^^^^^^^^^^^^^^
        #
        # The underlined annotation has to be stripped out by the cluster_key function,
        # so that we can use the original read ID to look up the cluster size.
        #
        read_id = read.header.split(None, 1)[0][1:]
        unique_count += 1
        nonunique_count += get_read_cluster_size(cluster_sizes, cluster_key(read_id))
    return unique_count, nonunique_count


def get_read_cluster_size(duplicate_cluster_sizes, read_id):
    suffix = None
    cluster_size = duplicate_cluster_sizes.get(read_id)
    if cluster_size == None:
        prefix, suffix = read_id[:-2], read_id[-2:]
        cluster_size = duplicate_cluster_sizes.get(prefix)
    assert cluster_size != None and suffix in (
        None, "/1", "/2"), f"Read ID not found in duplicate_cluster_sizes dict: {read_id}"
    return cluster_size


def _load_duplicate_cluster_sizes_work(filename):
    duplicate_cluster_sizes = {}
    with open(filename, "r") as f:
        for line in f:
            cluster_size_str, read_id = line.split(None, 1)
            duplicate_cluster_sizes[read_id.strip()] = int(cluster_size_str)
    return duplicate_cluster_sizes


# Loading cluster sizes can be expensive prior to subsampling (for some exceptionally large
# samples with over 100 million reads).  To ameliorate this cost, we make sure it is only
# paid once per stage (not once per step).
_DUPLICATE_CLUSTER_SIZES_CACHE = {}
_DUPLICATE_CLUSTER_SIZES_LOCK = multiprocessing.RLock()


def load_duplicate_cluster_sizes(filename):
    with _DUPLICATE_CLUSTER_SIZES_LOCK:
        if filename not in _DUPLICATE_CLUSTER_SIZES_CACHE:
            _DUPLICATE_CLUSTER_SIZES_CACHE[filename] = _load_duplicate_cluster_sizes_work(filename)
        return _DUPLICATE_CLUSTER_SIZES_CACHE[filename]


def save_duplicate_cluster_sizes(filename, duplicate_clusters):
    with _DUPLICATE_CLUSTER_SIZES_LOCK:
        _DUPLICATE_CLUSTER_SIZES_CACHE[filename] = {}
    with open(filename, "w") as tsv:
        for read_id, clusters in duplicate_clusters.items():
            cluster_size = clusters[0]
            tsv.write(f"{cluster_size}\t{read_id}\n")
            _DUPLICATE_CLUSTER_SIZES_CACHE[filename][read_id] = cluster_size


def reads_in_group(file_group, max_fragments=None, cluster_sizes=None, cluster_key=None):
    '''
    OVERVIEW

    Count reads in a group of matching files, up to a maximum number of fragments,
    optionally expanding each fragment to its duplicate cluster size.  Inputs in
    FASTA and FASTQ format are supported, subject to restrictions below.

    DEFINITIONS

    The term "fragment" refers to the physical DNA fragments the reads derive from.
    If the input is single, then 1 fragment == 1 read.
    If the input is paired, then 1 fragment == 2 reads.

    The term "cluster" refers to the output of the pipeline step idseq-dedup, which
    identifies and groups together duplicate fragments into clusters.  Duplicates
    typically result from PCR or other amplificaiton.  It is cost effective for the
    pipeline to operate on a single representative fragment from each cluster, then
    infer back the original fragment count using the cluster sizes information
    emitted by the idseq-dedup step.  The caller should load that information via
    m8.load_cluster_sizes() before passing it to this function.

    RESTRICTIONS

    When cluster_sizes is specified, the input must be in FASTA format, and the max_reads
    argument must be unspecified (None).  Thus, the caller should choose between two
    distinct use cases:

       (*) Expand clusters when counting FASTA fragments.

       (*) Count FASTA or FASTQ fragments without expanding clusters,
           optionally truncating to max_fragments.  The input may
           even be gz compressed.

    PERFORMANCE

    When cluster_sizes is not specified, the implementation is very fast, via wc.

    When cluster_sizes is specified, the python implementation can process only about
    3-4 million reads per minute. Fortunately, only steps run_lzw and run_bowtie2 can
    have more than a million fragments, and even that is extremely unlikely (deeply
    sequenced microbiome samples are rare in idseq at this time).  All other steps either
    operate on at most 1 million fragments or do not specify cluster_sizes. If this
    performance ever becomes a concern, the relevant private function can be reimplemented
    easily in GO.

    INTERFACE

    It may be better to expose the two use cases as separate functions, so that users do
    not encounter surprising restrictions or performance differences.  However, that would
    fracture the documentation, obfuscate points of use, and cement a separation that is
    quite arbitrary and forced only by current implementation choices.  Instead, keeping it
    all in one would allow a future reimplementation in a more performant language that
    could lift all restrictions and make all code paths performant.  The choice, then, is
    this better future, and just fail an assert where the present falls short.
    '''
    assert None in (max_fragments, cluster_sizes), "Truncating to max_fragments is not supported at the same time as expanding cluster_sizes.  Consider setting max_fragments=None."
    assert (cluster_sizes == None) == (cluster_key == None), "Please specify cluster_key when using cluster_sizes."
    first_file = file_group[0]
    # This is so fast, just do it always as a sanity check.
    unique_fast = count_reads(first_file)
    if max_fragments is not None:
        unique_fast = min(unique_fast, max_fragments)
    if cluster_sizes:
        # Run this even if ReadCountingMode.COUNT_UNIQUE to get it well tested before release.  Dark launch.
        unique, nonunique = _count_reads_expanding_duplicates(first_file, cluster_sizes, cluster_key)
        assert unique_fast == unique, f"Different read counts from wc ({unique_fast}) and fasta.iterator ({unique}) for file {first_file}."
        assert unique <= nonunique, f"Unique count ({unique}) should not exceed nonunique count ({nonunique}) for file {first_file}."
    reads_in_first_file = unique_fast
    if cluster_sizes and READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL:
        reads_in_first_file = nonunique
    num_files = len(file_group)
    return num_files * reads_in_first_file
