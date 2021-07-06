#!/usr/bin/env python3
"""
Post-processes the outputs of harvest.py to distill essential metrics to be tracked over time in a
data warehouse
"""

import sys
import argparse
import json
from datetime import datetime


def main(argv):
    parser = argparse.ArgumentParser(
        argv[0], description="harvest benchmark runs into importable JSON"
    )

    parser.add_argument(
        "harvest_json",
        metavar="HARVEST.json",
        nargs="+",
        help="captured output of harvest.py",
    )
    parser.add_argument(
        "--wall",
        metavar="SECONDS",
        type=int,
        help="elapsed processing time, in seconds",
    )
    parser.add_argument(
        "--version", metavar="vX.Y.Z", type=str, help="idseq-workflows release version"
    )

    args = parser.parse_args(argv[1:])
    harvest(**vars(args))


def harvest(harvest_json, wall, version):
    print(",".join(["version", "sample", "timestamp", "wall_seconds", "NT_aupr", "NR_aupr"]))

    timestamp = datetime.utcnow().isoformat()

    for one_harvest_json in harvest_json:
        with open(one_harvest_json, "r") as infile:
            for sample, results in json.load(infile).items():
                print(
                    ",".join(
                        str(v)
                        for v in [
                            version,
                            sample,
                            timestamp,
                            wall,
                            results["NT_aupr"],
                            results["NR_aupr"],
                        ]
                    )
                )


if __name__ == "__main__":
    main(sys.argv)
