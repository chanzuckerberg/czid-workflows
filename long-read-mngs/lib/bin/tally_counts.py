import argparse

import pandas as pd


def main(
    m8_filepath: str,
    hitsummary_filepath: str,
    reads_to_contigs_filepath: str,
    output_filepath: str,
):
    alignment_length_df = pd.read_csv(m8_filepath, sep="\t", index_col="read_id", names=[
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

    hitsummary_df = pd.read_csv(hitsummary_filepath, sep="\t", index_col="read_id", names=[
        "read_id",
        "x",
        "final_taxid",
        "accession",
        "species_taxid",
        "genus_taxid",
        "family_taxid",
    ], usecols=[
        "read_id",
        "final_taxid",
        "genus_taxid",
    ])
    hitsummary_df_n_rows = hitsummary_df.shape[0]

    # Add aln_len to the hitsummary df
    hitsummary_df = hitsummary_df.join(alignment_length_df[["alignment_length"]], how="inner", on="read_id")
    del alignment_length_df

    # Since we inner joined, if the resulting DataFrame is smaller it means there were some read_ids
    #   in hitsummary_df that were not in m8_df
    assert hitsummary_df.shape[0] == hitsummary_df_n_rows, "missing read_ids in m8 file"

    sp_data_collection = hitsummary_df[["final_taxid", "alignment_length"]]
    sp_counts = sp_data_collection.groupby(["final_taxid"]).sum()

    gen_data_collection = hitsummary_df[["genus_taxid", "alignment_length"]]
    gen_counts = gen_data_collection.groupby(["genus_taxid"]).sum()

    reads_to_contigs_df = pd.read_csv(reads_to_contigs_filepath, sep='\t', names=['seqid', 'contig', 'sequence'])
    reads_to_contigs_df['seq_len'] = [len(i) for i in reads_to_contigs_df['sequence']]
    reads_to_contigs_df = reads_to_contigs_df[reads_to_contigs_df.contig != '*']
    contig_tally_df = reads_to_contigs_df[['contig', 'seq_len']]
    result = contig_tally_df.groupby(['contig']).sum()

    df = pd.merge(hitsummary_df, result, how='outer', left_on="read_id", right_on="contig")
    sp_df = df[['final_taxid', 'seq_len']]
    gen_df = df[['genus_taxid', 'seq_len']]
    print(sp_df.index.dtype)

    sp_result_final = sp_df.groupby(["final_taxid"]).sum().join(sp_counts, on="final_taxid")
    sp_result_final.index.name = "taxid"
    sp_result_final.insert(0, "level", "species")

    gen_result_final = gen_df.groupby(["genus_taxid"]).sum().join(gen_counts, on="genus_taxid")
    gen_result_final.index.name = "taxid"
    gen_result_final.insert(0, "level", "genus")

    pd.concat([sp_result_final, gen_result_final], axis=0).to_csv(output_filepath)


parser = argparse.ArgumentParser()
parser.add_argument("--m8-filepath")
parser.add_argument("--hitsummary-filepath")
parser.add_argument("--reads-to-contigs-filepath")
parser.add_argument("--output-filepath")

if __name__ == "__main__":
    args = parser.parse_args()
    main(
        m8_filepath=args.m8_filepath,
        hitsummary_filepath=args.hitsummary_filepath,
        reads_to_contigs_filepath=args.reads_to_contigs_filepath,
        output_filepath=args.output_filepath,
    )
