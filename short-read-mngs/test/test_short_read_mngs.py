import json
import atexit


def test_bench3_viral(short_read_mngs_bench3_viral_outputs):
    # short_read_mngs_bench3_viral_outputs is a fixture defined in ./conftest.py providing the JSON outputs of
    # short-read-mngs executed on the bench3 fastq's & viral reference databases. Using the (session-scoped) fixture
    # causes that test workflow run to execute.
    outp = short_read_mngs_bench3_viral_outputs
    atexit.register(lambda: print(f"short_read_mngs_bench3_viral run dir: {outp['dir']}"))

    # check correctness of the workflow outputs
    with open(
        outp["outputs"][
            "idseq_short_read_mngs.postprocess.refined_taxon_count_out_assembly_refined_taxon_counts_with_dcr_json"
        ]
    ) as infile:
        taxon_counts = json.load(infile)["pipeline_output"]["taxon_counts_attributes"]

    taxa = set(entry["tax_id"] for entry in taxon_counts)
    assert len(taxa) == 149

    # TODO: further correctness tests
