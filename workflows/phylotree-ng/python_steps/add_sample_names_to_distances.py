import argparse
import json
from csv import DictReader, DictWriter
from typing import TypedDict, Iterable


class Sample(TypedDict):
    sample_name: str
    workflow_run_id: int
    contig_fasta: str
    combined_contig_summary: str


def main(distances: str, output_distances: str, samples: Iterable[Sample]):
    sample_name_by_workflow_run_id = {str(s["workflow_run_id"]): s["sample_name"] for s in samples}
    with open(distances) as r, open(output_distances, "w") as w:
        reader = DictReader(r, delimiter="\t")
        writer = None
        for i, row in enumerate(reader):
            if i == 0:
                assert reader.fieldnames
                writer = DictWriter(w, fieldnames=reader.fieldnames, delimiter="\t")
                writer.writeheader()
            assert writer
            row["Sample 1"] = sample_name_by_workflow_run_id.get(row["Sample 1"], row["Sample 1"])
            row["Sample 2"] = sample_name_by_workflow_run_id.get(row["Sample 2"], row["Sample 2"])
            writer.writerow(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--distances")
    parser.add_argument("--samples")
    parser.add_argument("--output-distances")
    args = parser.parse_args()

    with open(args.samples) as f:
        samples: Iterable[Sample] = json.load(f)

    main(args.distances, args.output_distances, samples)
