import argparse
import csv
import logging
import os
from os.path import join
from subprocess import run
from typing import List

import marisa_trie
from Bio import SeqIO

# Number of sequences in a taxid to trigger compression
COMPRESSION_THRESHOLD = 12


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def greedy_partition(csv_path: str, threshold: float):
    clusters_dict = {}
    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        labels = next(reader)
        for i, row in enumerate(reader):
            if i in clusters_dict:
                continue

            cluster = [j for j, value in enumerate(row) if float(value) > threshold and j not in clusters_dict]
            for j in cluster:
                clusters_dict[j] = None
            clusters_dict[i] = cluster

    return set(labels[i] for i, cluster in clusters_dict.items() if cluster)


def compress_taxid(dirname: str, taxid: int, k: int, st: int):
    taxid_dir = join(dirname, str(taxid))
    output_path = join(taxid_dir, "compressed.fasta")

    with open(join(taxid_dir, "count"), "r") as f:
        count = int(f.read())
        if count < COMPRESSION_THRESHOLD:
            os.rename(join(taxid_dir, "input.fasta"), output_path)
            return output_path

    run(["sourmash", "sketch", "dna", "-p", f"k={k},scaled={st}", "--singleton", "--output", "sig", "input.fasta"], cwd=taxid_dir, check=True)
    run(["sourmash", "compare", "sig", "--containment", "-o", "cmp.dist", "--csv", "cmp.csv"], cwd=taxid_dir, check=True)

    labels_to_keep = greedy_partition(join(taxid_dir, "cmp.csv"), 0.6)

    with open(join(taxid_dir, "input.fasta"), "r") as f_r, open(output_path, "w") as f_w:
        for record in SeqIO.parse(f_r, "fasta"):
            if record.id in labels_to_keep:
                SeqIO.write(record, f_w, "fasta")

    return output_path


def compress_nt(nt_filepath: str, accession2taxid_path: str, taxids_to_drop: List[int], k: int, st: int, compressed_nt_filepath: str):
    # Split NT by taxid
    accession2taxid = marisa_trie.RecordTrie("L").mmap(accession2taxid_path)
    dirname = "nt_by_taxid"
    with open(nt_filepath) as f_r:
        for record in SeqIO.parse(f_r, "fasta"):
            non_versioned_accession_id = record.id.rsplit(".", 1)[0]
            if non_versioned_accession_id not in accession2taxid:
                continue

            taxid = accession2taxid[non_versioned_accession_id][0][0]
            if taxid in taxids_to_drop:
                continue

            os.makedirs(join(dirname, str(taxid)), exist_ok=True)
            with open(join(dirname, str(taxid), "input.fasta"), "a") as f_w:
                SeqIO.write(record, f_w, "fasta")

            if not os.path.exists(join(dirname, str(taxid), "count")):
                count = 0
            else:
                with open(join(dirname, str(taxid), "count") , "r") as f:
                    count = int(f.read())
            with open(join(dirname, str(taxid), "count"), "w") as f:
                f.write(str(count + 1))


    for taxid in os.listdir(dirname):
        taxid = int(taxid)
        compressed_fasta = compress_taxid(dirname, taxid, k, st)

        logger.info(f"Compressed taxid {taxid} to {compressed_fasta}")
        logger.info(f"{compressed_nt_filepath}")

        with open(compressed_fasta, "r") as f_r, open(compressed_nt_filepath, "a") as f_w:
            f_w.writelines(f_r)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nt-filepath")
    parser.add_argument("--accession2taxid-path")
    parser.add_argument("--taxids-to-drop", nargs="+", type=int)
    parser.add_argument("--k", type=int)
    parser.add_argument("--st", type=int)
    parser.add_argument("--compressed-nt-filepath")
    args = parser.parse_args()

    compress_nt(args.nt_filepath, args.accession2taxid_path, args.taxids_to_drop or [], args.k, args.st, args.compressed_nt_filepath)
