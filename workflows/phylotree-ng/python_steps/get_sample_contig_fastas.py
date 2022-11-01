import argparse
import json
from typing import TypedDict, Iterable, Set

from Bio import SeqIO


class Sample(TypedDict):
    sample_name: str
    workflow_run_id: int
    contig_fasta: str
    combined_contig_summary: str


def get_contig_ids(reference_taxid: str, contig_summary: str):
    with open(contig_summary) as f:
        for row in json.load(f):
            if row["taxid"] == reference_taxid:
                for contig_id in row["contig_counts"].keys():
                    yield contig_id


def get_records(contig_ids: Set[str], contig_fasta: str):
    for record in SeqIO.parse(contig_fasta, "fasta"):
        if record.id in contig_ids:
            yield record


def main(reference_taxid: str, samples: Iterable[Sample]):
    for sample in samples:
        contig_ids = set(get_contig_ids(reference_taxid, sample["combined_contig_summary"]))
        SeqIO.write(
            get_records(contig_ids, sample["contig_fasta"]),
            # It is critical that these be named with the workflow_run_id so the workflow_run_id
            # ends up in the resulting phylo tree
            f"{sample['workflow_run_id']}.fasta",
            "fasta",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-taxid")
    parser.add_argument("--samples")
    args = parser.parse_args()

    with open(args.samples) as f:
        samples: Iterable[Sample] = json.load(f)

    main(args.reference_taxid, samples)
