import argparse

from Bio import SeqIO

# Function to determine the appropriate file based on sequence length
def get_file_name(length):
    min_length = int(length/1000)*1000
    max_length = min_length + 1000

    # pad with leading zeros so that we can concatenate the files together in order with cat
    min_length_padded = str(min_length).zfill(12)
    max_length_padded = str(max_length).zfill(12)

    return f"sequences_{min_length_padded}-{max_length_padded}.fa"


def break_apart_fasta(fn):
    # Iterate through sequences and append to appropriate files
    with open(fn, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            seq_length = len(record.seq)
            file_name = get_file_name(seq_length)
            with open(file_name, "a") as output_handle:
                SeqIO.write(record, output_handle, "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Break apart fasta file into smaller files based on sequence length')
    parser.add_argument('--fasta-file', required=True)
    args = parser.parse_args()

    break_apart_fasta(args.fasta_file)
