import os
import json
from Bio import SeqIO


def test_RunAssembly_defaults(util, short_read_mngs_bench3_viral_outputs):
    """
    On default settings, the assembly_out_contigs_fasta and assembly_out_contigs_all_fasta files
    should be identical (because all contigs for this test dataset pass the default
    min_contig_length filter)
    """
    assembly_contigs_fasta = short_read_mngs_bench3_viral_outputs["outputs"][
        "idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_fasta"
    ]
    assembly_contigs_all_fasta = short_read_mngs_bench3_viral_outputs["outputs"][
        "idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_all_fasta"
    ]
    assembly_contigs = list(SeqIO.parse(assembly_contigs_fasta, "fasta"))
    assembly_contigs_all = list(SeqIO.parse(assembly_contigs_all_fasta, "fasta"))
    assert [(r.id, r.seq) for r in assembly_contigs] == [
        (r.id, r.seq) for r in assembly_contigs_all
    ]


def test_RunAssembly_filtered(util, short_read_mngs_bench3_viral_outputs):
    """
    Re-run the assembly task with a higher min_contig_length and observe that contigs are removed.
    """
    # load the task's inputs from the end-to-end workflow test
    task_name = "RunAssembly"
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"],
            "call-postprocess",
            f"call-{task_name}",
        )
    )

    # override min_contig_length
    inputs["min_contig_length"] = 150

    # rerun
    outp = util.miniwdl_run(
        util.repo_dir() / "short-read-mngs/postprocess.wdl",
        "--task",
        task_name,
        "-i",
        json.dumps(inputs),
    )

    # verify expected output subset
    assembly_contigs_fasta = os.path.join(
        outp["dir"], outp["outputs"][f"{task_name}.assembly_contigs_fasta"]
    )
    assembly_contigs_all_fasta = os.path.join(
        outp["dir"], outp["outputs"][f"{task_name}.assembly_contigs_all_fasta"]
    )
    contigs = list(SeqIO.parse(assembly_contigs_fasta, "fasta"))
    contigs_all = list(SeqIO.parse(assembly_contigs_all_fasta, "fasta"))
    assert contigs and all([len(r.seq) >= 150 for r in contigs])
    assert set(r.id for r in contigs_all) - set(r.id for r in contigs)
    assert all(r.seq == next(s.seq for s in contigs_all if s.id == r.id) for r in contigs)
