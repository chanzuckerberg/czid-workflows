import csv
import re
import argparse

front_pattern = re.compile(r"^\d+[SH]")
end_pattern = re.compile(r"\d+[SH]$")


def split_cigar(cigar_str):
    """splits the count from type in a cigar string"""
    return int(cigar_str[:-1]), cigar_str[-1]


def get_clipping(cigar, alength, pattern):
    """Get the amount of clipping from the front or end of the cigar string"""
    clipping = 0

    match = re.search(pattern, cigar)
    if match:
        clip_length, clip_type = split_cigar(match.group())
        clipping += clip_length
        if clip_type == "H":
            alength += clip_length

    return clipping, alength


def main(reads_to_contig_sam, output_file, max_percent):
    with open(reads_to_contig_sam, "r") as csv_file, open(
        output_file, mode="w", newline=""
    ) as output_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")
        output = csv.writer(
            output_file,
            delimiter="\t",
            lineterminator="\n",
            quotechar="\t",
            quoting=csv.QUOTE_NONE,
        )

        # Loop through each row in the CSV file
        for row in csv_reader:
            if len(row) > 1 and len(row[0]) > 1 and row[0][0] == "@":
                output.writerow(row)
            elif row[5] != "*":
                cigar = row[5]
                alength = len(row[9])
                front_clipping, alength = get_clipping(cigar, alength, front_pattern)
                end_clipping, alength = get_clipping(cigar, alength, end_pattern)
                percent_clipped = ((front_clipping + end_clipping) / alength) * 100
                if percent_clipped < max_percent:
                    output.writerow(row)
            else:
                output.writerow(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Filter clipped alignments",
        description="Filter out alignments from minimap2 that are too heavily clipped",
    )
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("-m", "--max_percent", default=50.0, type=float)
    args = parser.parse_args()
    main(args.input_file, args.output_file, args.max_percent)
