import re

import click
import pandas as pd
import pysam

# NOTE: Contigs are indexed by ARO accession, not gene_id. This is because contigs and reads
# come from different sources, so not every contig will have a corresponding gene_id,
# even if you join the two datasets on the ARO accession. Thus indexing on gene_id is not
# possible.

COLUMN_ARO_CONTIG = "ARO_contig_amr"
COLUMN_CONTIG_NAME = "Contig_contig_amr"
OUTPUT_BAM = "contig_amr_report.sorted.bam"


# pysam does not play nice with file objects, so we have to require file paths as inputs
@click.command()
@click.option("contigs_file", "--contigs", type=click.Path(exists=True), required=True)
@click.option("final_summary", "--final-summary", type=click.Path(exists=True), required=True)
def cli(contigs_file: click.Path, final_summary: click.Path):
    if not spades_assembly_succeeded(contigs_file):
        create_empty_bam_bai()
        return
    create_fasta_index(contigs_file)
    contigs_fasta = load_contigs_fasta(contigs_file)
    df = final_summary_dataframe(final_summary)
    write_bam_bai_pair(contigs_fasta, df)
    return


def spades_assembly_succeeded(contigs_file: click.Path) -> bool:
    with open(contigs_file) as f:
        first_line = f.readline()
        if first_line == ";ASSEMBLY FAILED\n":
            return False
    return True


# Create an empty BAM/BAI file if the SPADES assembly failed, then exit
def create_empty_bam_bai():
    output_bam = pysam.AlignmentFile(OUTPUT_BAM, "wb", reference_names=["NoGenes"], reference_lengths=[100])
    output_bam.close()
    pysam.index(OUTPUT_BAM)


# Create .fai index to enable querying fasta file
def create_fasta_index(contigs_file: click.Path):
    # Equivalent to the samtools faidx command
    pysam.faidx(contigs_file)


def load_contigs_fasta(contigs_file: click.Path) -> pysam.Fastafile:
    contigs_fasta = pysam.Fastafile(contigs_file)
    return contigs_fasta


def final_summary_dataframe(final_summary: click.Path) -> pd.DataFrame:
    # Load columns of interest from CSV and drop rows with at least one NaN
    df = pd.read_csv(final_summary, sep="\t", usecols=[COLUMN_ARO_CONTIG, COLUMN_CONTIG_NAME])
    df = df.dropna(subset=COLUMN_CONTIG_NAME)

    # Format accessions to follow the pattern 'ARO:3000000'; Remove extraneous _* at the end of contig names
    df[COLUMN_ARO_CONTIG] = df[COLUMN_ARO_CONTIG].apply(lambda x: f'ARO:{int(x)}')
    df[COLUMN_CONTIG_NAME] = df[COLUMN_CONTIG_NAME].apply(lambda x: re.sub(r'_\d*$', '', x))
    return df


def write_bam_bai_pair(contigs_fasta: pysam.Fastafile, df: pd.DataFrame):
    # Create BAM with mock reference lengths for the header. If contigss are found at all, have a mock accession
    # to make sure we can create the SAM file with no errors (web app will look for that file)
    contig_aros = df[COLUMN_ARO_CONTIG].dropna().unique().tolist() or ["NoGenes"]
    output_bam = pysam.AlignmentFile(
        OUTPUT_BAM,
        "wb",
        reference_names=contig_aros,
        reference_lengths=[100] * len(contig_aros)
    )

    # Go through each line of the TSV and create a SAM record (https://wckdouglas.github.io/2021/12/pytest-with-pysam)
    for index, row in df.iterrows():
        contig_aro = row[COLUMN_ARO_CONTIG]
        contig_name = row[COLUMN_CONTIG_NAME]
        contig_sequence = contigs_fasta.fetch(contig_name)

        # Create new alignment
        alignment = pysam.AlignedSegment(output_bam.header)
        alignment.reference_name = contig_aro
        alignment.query_name = contig_name
        alignment.query_sequence = contig_sequence
        alignment.reference_start = 1
        output_bam.write(alignment)

    output_bam.close()
    pysam.index(OUTPUT_BAM)


if __name__ == '__main__':
    cli()
