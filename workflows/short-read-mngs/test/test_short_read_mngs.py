import json
import atexit
import re
import os.path
from Bio import SeqIO


def is_valid_fasta(filename):
    try:
        with open(filename) as f:
            for read in SeqIO.parse(f, "fasta"):
                pass
        return True
    except Exception:
        return False


def test_bench3_viral(short_read_mngs_bench3_viral_outputs):
    # short_read_mngs_bench3_viral_outputs is a fixture defined in ./conftest.py providing the JSON outputs of
    # short-read-mngs executed on the bench3 fastq's & viral reference databases. Using the (session-scoped) fixture
    # causes that test workflow run to execute.
    outp = short_read_mngs_bench3_viral_outputs
    atexit.register(lambda: print(f"short_read_mngs_bench3_viral run dir: {outp['dir']}"))

    # check correctness of the workflow outputs
    with open(
        outp["outputs"][
            "czid_short_read_mngs.postprocess.refined_taxon_count_out_assembly_refined_taxon_counts_with_dcr_json"
        ]
    ) as infile:
        taxon_counts = json.load(infile)["pipeline_output"]["taxon_counts_attributes"]

    taxa = set(entry["tax_id"] for entry in taxon_counts)
    assert abs(len(taxa) - 184) < 16

    for filename in outp["outputs"]:
        if filename.endswith(".fasta"):
            assert is_valid_fasta(filename), f"{filename} is not a valid fasta file"

    longest_reads = outp["outputs"]["czid_short_read_mngs.experimental.longest_reads"]
    basenames = [os.path.basename(fn) for fn in longest_reads]
    assert basenames, basenames
    assert all(re.match(r"n[rt]\.[a-z]+\.-?[0-9]+\.longest_5_reads.fasta", fn) for fn in basenames), basenames
    prefixes = set(fn[:2] for fn in basenames)
    assert "nt" in prefixes, f"'nt' not found in {prefixes}"
    assert "nr" in prefixes, f"'nr' not found in {prefixes}"

    for fn in longest_reads:
        with open(fn) as f:
            lines = list(f)
            assert 2 <= len(lines) <= 10, len(lines)
            prev = None
            for i, read in enumerate(lines):
                if i % 2 == 0:
                    assert read[0] == ">", read
                    continue
                assert prev is None or len(read) <= prev, (len(read), prev)
                prev = len(read)
                assert all(c in "ACTGUN" for c in read.strip()), read
    # TODO: further correctness tests
