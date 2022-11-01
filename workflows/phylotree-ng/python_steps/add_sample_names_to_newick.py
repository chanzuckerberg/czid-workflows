import argparse
import json
from typing import TypedDict, Iterable

from Bio.Phylo import NewickIO


class Sample(TypedDict):
    sample_name: str
    workflow_run_id: int
    contig_fasta: str
    combined_contig_summary: str


def main(newick: str, output_newick: str, samples: Iterable[Sample]):
    sample_name_by_workflow_run_id = {str(s["workflow_run_id"]): s["sample_name"] for s in samples}
    with open(newick) as i, open(output_newick, "w") as o:
        tree = next(NewickIO.parse(i))
        for node in tree.find_clades(order="level"):
            node.name = sample_name_by_workflow_run_id.get(node.name, node.name)
        NewickIO.write([tree], o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--newick")
    parser.add_argument("--samples")
    parser.add_argument("--output-newick")
    args = parser.parse_args()

    with open(args.samples) as f:
        samples: Iterable[Sample] = json.load(f)

    main(args.newick, args.output_newick, samples)
