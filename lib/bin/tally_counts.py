import argparse

import pandas as pd
from Bio import SeqIO


def main(
    reads_fastq_filepath: str,
    m8_filepath: str,
    hitsummary_filepath: str,
    reads_to_contigs_filepath: str,
    output_filepath: str,
):
    alignment_length = pd.read_csv(m8_filepath, sep="\t", index_col="read_id", names=[
        "read_id",
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
    ], usecols=["read_id", "alignment_length"])

    hitsummary = pd.read_csv(hitsummary_filepath, sep="\t", index_col="read_id", names=[
        "read_id",
        "x",
        "final_taxid",
        "accession",
        "species_taxid",
        "genus_taxid",
        "family_taxid",
    ], usecols=[
        "read_id",
        "species_taxid",
        "genus_taxid",
    ])
    hitsummary_n_rows = hitsummary.shape[0]

    # Add aln_len to the hitsummary df
    hitsummary = hitsummary.join(alignment_length[["alignment_length"]], how="inner", on="read_id")

    # Since we inner joined, if the resulting DataFrame is smaller it means there were some read_ids
    #   in hitsummary_df that were not in m8_df
    assert hitsummary.shape[0] == hitsummary_n_rows, "missing read_ids in m8 file"

    species_alignment_lengths = hitsummary[["species_taxid", "alignment_length"]].groupby(["species_taxid"]).sum()
    species_alignment_lengths.index.name = "taxid"
    species_alignment_lengths.columns = ["total_alignment_length"]
    genus_alignment_lengths = hitsummary[["genus_taxid", "alignment_length"]].groupby(["genus_taxid"]).sum()
    genus_alignment_lengths.index.name = "taxid"
    genus_alignment_lengths.columns = ["total_alignment_length"]

    reads_lengths = pd.DataFrame(
        {"read_id": read.id, "read_length": len(read.seq)} for read in SeqIO.parse(reads_fastq_filepath, "fastq")
    ).set_index("read_id")

    reads_to_contigs = pd.read_csv(reads_to_contigs_filepath, sep="\t", names=[
        "read_id",
        "contig_id",
        "alignment",
    ])
    reads_to_contigs = reads_to_contigs[reads_to_contigs.contig_id != "*"]
    reads_to_contigs["alignment_length"] = reads_to_contigs["alignment"].str.len()
    # we want to only keep the longest alignment for each read so reads are not double counted
    reads_to_contigs = reads_to_contigs.sort_values("alignment_length", ascending=False).drop_duplicates(["read_id"])
    reads_to_contigs = reads_to_contigs[["read_id", "contig_id"]].set_index("read_id")

    reads_to_contigs = reads_to_contigs.join(reads_lengths, how="inner", on="read_id")
    contig_sequence_lengths = reads_to_contigs[[
        "contig_id",
        "read_length",
    ]].groupby(["contig_id"]).sum()
    contig_sequence_lengths.columns = ["total_sequence_length"]

    taxid_sequence_lengths = pd.merge(
        hitsummary,
        contig_sequence_lengths,
        how="left",
        left_on="read_id",
        right_on="contig_id",
    )

    species_result_final = taxid_sequence_lengths[[
        "species_taxid",
        "total_sequence_length",
    ]].groupby(["species_taxid"]).sum().join(species_alignment_lengths, on="species_taxid")
    species_result_final.index.name = "taxid"
    species_result_final.insert(0, "level", "species")

    genus_result_final = taxid_sequence_lengths[[
        "genus_taxid",
        "total_sequence_length",
    ]].groupby(["genus_taxid"]).sum().join(genus_alignment_lengths, on="genus_taxid")
    genus_result_final.index.name = "taxid"
    genus_result_final.insert(0, "level", "genus")

    pd.concat([
        species_result_final,
        genus_result_final,
    ], axis=0).sort_values(by="total_alignment_length", ascending=False).to_csv(output_filepath)


parser = argparse.ArgumentParser()
parser.add_argument("--reads-fastq-filepath")
parser.add_argument("--m8-filepath")
parser.add_argument("--hitsummary-filepath")
parser.add_argument("--reads-to-contigs-filepath")
parser.add_argument("--output-filepath")

if __name__ == "__main__":
    args = parser.parse_args()
    main(
        reads_fastq_filepath=args.reads_fastq_filepath,
        m8_filepath=args.m8_filepath,
        hitsummary_filepath=args.hitsummary_filepath,
        reads_to_contigs_filepath=args.reads_to_contigs_filepath,
        output_filepath=args.output_filepath,
    )
