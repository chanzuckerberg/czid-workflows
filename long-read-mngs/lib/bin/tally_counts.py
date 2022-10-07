import argparse

import pandas as pd


def main(
    m8_filepath: str,
    hitsummary_filepath: str,
    reads_to_contigs_filepath: str,
    species_output_filepath: str,
    genus_output_filepath: str,
):
    m8_df = pd.read_csv(m8_filepath, sep="\t", index_col="read_id", names=[
        "read_id",
        "accession",
        "pid",
        "aln_len",
        "mismatches",
        "gap_openings",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "evalue",
        "bit_score",
    ])
    hitsummary_df = pd.read_csv(hitsummary_filepath, sep="\t", index_col="read_id", names=[
        "read_id",
        "x",
        "final_taxid",
        "accession",
        "species_taxid",
        "genus_taxid",
        "family_taxid",
    ])

    hitsummary_aln_df = hitsummary_df.join(m8_df[["aln_len"]], how="inner", on="read_id")

    assert hitsummary_df.shape[0] == hitsummary_aln_df.shape[0], "missing read_ids in m8 file"

    sp_data_collection = hitsummary_aln_df[["final_taxid", "aln_len"]]
    # sp_counts = sp_data_collection.groupby(["final_taxid"]).sum()

    gen_data_collection = hitsummary_aln_df[["genus_taxid", "aln_len"]]
    # gen_counts = gen_data_collection.groupby(["genus_taxid"]).sum()

    r2c_df = pd.read_csv(reads_to_contigs_filepath, sep='\t', names=['seqid', 'contig', 'sequence'])
    r2c_df['seq_len'] = [len(i) for i in r2c_df['sequence']]
    r2c_df = r2c_df[r2c_df.contig != '*']
    contig_tally_df = r2c_df[['contig', 'seq_len']]
    result = contig_tally_df.groupby(['contig']).sum()

    df = pd.merge(hitsummary_aln_df, result, how='outer', left_on="read_id", right_on="contig")
    sp_df = df[['final_taxid', 'seq_len']]
    gen_df = df[['genus_taxid', 'seq_len']]

    # these should contain the sorted species / genus counts by length (bp)
    sp_df.groupby(["final_taxid"]).sum().to_csv(species_output_filepath)
    gen_df.groupby(["genus_taxid"]).sum().to_csv(genus_output_filepath)


parser = argparse.ArgumentParser()
parser.add_argument("m8-filepath")
parser.add_argument("hitsummary-filepath")
parser.add_argument("reads-to-contigs-filepath")
parser.add_argument("species-output-filepath")
parser.add_argument("genus-output-filepath")

if __name__ == "__main__":
    args = parser.parse_args()
    main(**args)
