import argparse
import json
from typing import TypedDict, Iterable

from Bio import SeqIO


class Sample(TypedDict):
    sample_name: str
    workflow_run_id: int
    contig_fasta: str
    combined_contig_summary: str


def add_sample_names_to_records(variants: str, samples: Iterable[Sample]):
    sample_name_by_workflow_run_id = {str(s["workflow_run_id"]): s["sample_name"] for s in samples}
    for record in SeqIO.parse(variants, "fasta"):
        record.id = sample_name_by_workflow_run_id.get(record.id, record.id)
        yield record


def main(variants: str, output_variants: str, samples: Iterable[Sample]):
    SeqIO.write(add_sample_names_to_records(variants, samples), output_variants, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants")
    parser.add_argument("--samples")
    parser.add_argument("--output-variants")
    args = parser.parse_args()

    with open(args.samples) as f:
        samples: Iterable[Sample] = json.load(f)

    main(args.variants, args.output_variants, samples)
