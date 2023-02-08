import argparse
import json
from pathlib import Path
import subprocess
import sys


def process_md5sum_file(filepath: str):
    """using md5sum here so I don't have to read in the files directly, change if it causes too many issues"""
    result = subprocess.run(["md5sum", filepath], capture_output=True).stdout.decode(
        "utf-8"
    )
    return result.split("  ")[0]


def run_md5_dict(output_dict: dict):
    """takes in a dict of key: value where value is either a file or a list of files or none"""
    key_to_hash = {}
    for key, value in output_dict.items():
        if isinstance(value, list):
            new_list = []
            for item in value:
                if Path(item).exists():
                    new_list.append(process_md5sum_file(item))
                else:
                    new_list.append(item)
            key_to_hash[key] = new_list
        elif value and Path(value).exists():
            key_to_hash[key] = process_md5sum_file(value)
        else:
            """if the input is not a file, compare the value directly"""
            key_to_hash[key] = value
    return key_to_hash


def compare_hashes(out_hash_1: dict, out_hash_2: dict, verbose: bool = True):
    keys_in_old_only = out_hash_1.keys() - out_hash_2.keys()
    keys_in_new_only = out_hash_2.keys() - out_hash_1.keys()

    if keys_in_old_only and verbose:
        print(f"outputs removed in new run {keys_in_old_only}\n")
    elif keys_in_old_only:
        print(*keys_in_old_only, sep="\n")

    if keys_in_new_only and verbose:
        print(f"outputs introduced in new run {keys_in_new_only}\n")
    elif keys_in_new_only:
        print(*keys_in_new_only, sep="\n")

    for key, hash_val in out_hash_1.items():
        if new_hash_val := out_hash_2.get(key, None):
            if hash_val != new_hash_val and verbose:
                print(f"{key} differs - old {hash_val} - new {new_hash_val}\n")
            elif hash_val != new_hash_val:
                print(key)


def create_md5_dict(outfile: str):
    with open(outfile) as f:
        output = json.load(f)
    return run_md5_dict(output)


def compare_outputs(outfile1: str, outfile2: str, verbose: bool = True):
    out_hash_1 = create_md5_dict(outfile1)
    out_hash_2 = create_md5_dict(outfile2)

    compare_hashes(out_hash_1, out_hash_2, verbose)


def print_json_hash(outfile: str):
    out_hash = create_md5_dict(outfile)
    print(json.dumps(out_hash))


def main():
    # set up the argument parser
    parser = argparse.ArgumentParser(
        description="A script that reads in outputs.json files from miniwdl and runs md5sum on the files. Prints out a diff of the files. Depends on md5sum. "
    )
    parser.add_argument(
        "miniwdl_output_1",
        metavar="miniwdl-output-1",
        type=str,
        help="first output json. If this is the only argument runs md5sum on the outputs and outputs the k/v file",
    )

    parser.add_argument(
        "miniwdl_output_2",
        metavar="miniwdl-output-2",
        type=str,
        nargs="?",
        help="Second output json. If present compares the first and second output files.",
    )

    parser.add_argument(
        "-v", "--verbose", help="Enable verbose output", action="store_true"
    )

    # parse the command line arguments
    args = parser.parse_args()

    # process the arguments
    if args.miniwdl_output_2:
        compare_outputs(args.miniwdl_output_1, args.miniwdl_output_2, args.verbose)
    else:
        print_json_hash(args.miniwdl_output_1)


if __name__ == "__main__":
    main()
