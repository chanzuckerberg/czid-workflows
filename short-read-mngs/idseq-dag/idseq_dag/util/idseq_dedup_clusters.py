"""
idseq-dedup outputs a cluster file in the form of a csv
The first column contains the representative read id of
a cluster, and the second column contains the read id.
"""
from csv import DictReader
from typing import Dict, Optional, List
import idseq_dag.util.log as log
from subprocess import run, PIPE

def parse_clusters_file(
    idseq_dedup_clusters_path: str,
) -> Dict[str, Optional[List]]:
    clusters_dict = {}
    log.write(f"opening clusters file {idseq_dedup_clusters_path}")
    lines = run(["wc", "-l", idseq_dedup_clusters_path], stdout=PIPE, check=True).stdout.split()[0]
    log.write(f"lines in cluster file: {lines}")
    with open(idseq_dedup_clusters_path) as f:
        log.write(f"opened clusters file {idseq_dedup_clusters_path}")
        for i, row in enumerate(DictReader(f)):
            if i % 100 == 0:
                log.write(f"parsed row {i}/{lines}")
            r_read_id, read_id = row["representative read id"], row["read id"]
            if r_read_id not in clusters_dict:
                clusters_dict[r_read_id] = [1]
            else:
                count, *others = clusters_dict[r_read_id]
                clusters_dict[r_read_id][0] += 1
                clusters_dict[r_read_id].append(read_id)
    return clusters_dict
