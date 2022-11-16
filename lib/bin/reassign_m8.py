import argparse

import pandas as pd


def main(
    m8_filepath: str,
    reads_to_contigs_filepath: str,
    output_filepath: str,
):
    m8 = pd.read_csv(m8_filepath, sep="\t", names=[
        "read_or_contig_id",
        "accession",
        "pid",
        "alignment_length",
        "mismatches",
        "gap_openings",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "evalue",
        "bit_score",
    ])

    reads_to_contigs = pd.read_csv(reads_to_contigs_filepath, sep="\t", names=[
        "read_id",
        "contig_id",
        "alignment",
    ], usecols=[
        "read_id",
        "contig_id",
    ])
    reads_to_contigs = reads_to_contigs[reads_to_contigs.contig_id != "*"]
    reassigned = pd.merge(m8, reads_to_contigs, left_on="read_or_contig_id", right_on="contig_id", how="left")
    reassigned["read_or_contig_id"] = reassigned.apply(lambda row: row["read_id"] if not pd.isnull(row["read_id"]) else row["read_or_contig_id"], axis=1)
    reassigned = reassigned.drop(columns=["read_id", "contig_id"])
    reassigned.to_csv(output_filepath, header=False, sep="\t", index=False)


parser = argparse.ArgumentParser()
parser.add_argument("--m8-filepath")
parser.add_argument("--reads-to-contigs-filepath")
parser.add_argument("--output-filepath")

if __name__ == "__main__":
    args = parser.parse_args()
    main(
        m8_filepath=args.m8_filepath,
        reads_to_contigs_filepath=args.reads_to_contigs_filepath,
        output_filepath=args.output_filepath,
    )
