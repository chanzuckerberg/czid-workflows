import argparse
import json

from Bio import SeqIO
from scipy import stats
import numpy as np


def read_length_metrics(
    fastq_path: str,
    json_output_path: str,
):
    read_lengths = np.array([len(seq.seq) for seq in SeqIO.parse(fastq_path, "fastq")])
    mean = np.mean(read_lengths)
    # NOTE: int casts to convert from numpy int64 to python int because numpy int64 is not
    #   JSON serializable
    read_length_stats = {
        "read_length_median": np.median(read_lengths),
        "read_length_mode": int(stats.mode(read_lengths).mode[0]),
        "read_length_absolute_deviation": np.sum(np.abs(read_lengths - mean)) / len(read_lengths),
        "read_length_min": int(np.min(read_lengths)),
        "read_length_max": int(np.max(read_lengths)),
        "read_length_mean": mean,
        "read_length_standard_deviation": np.std(read_lengths),
    }

    with open(json_output_path, 'w') as f:
        json.dump(read_length_stats, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq-path")
    parser.add_argument("--json-output-path")
    args = parser.parse_args()
    read_length_metrics(args.fastq_path, args.json_output_path)
