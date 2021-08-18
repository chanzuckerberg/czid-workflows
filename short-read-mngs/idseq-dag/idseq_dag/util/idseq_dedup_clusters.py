"""
idseq-dedup outputs a cluster file in the form of a csv
The first column contains the representative read id of
a cluster, and the second column contains the read id.
"""
from csv import DictReader
from typing import Dict, Optional, List


def parse_clusters_file(
    idseq_dedup_clusters_path: str,
) -> Dict[str, Optional[List]]:
    clusters_dict = {}
    with open(idseq_dedup_clusters_path) as f:
        for row in DictReader(f):
            r_read_id, read_id = row["representative read id"], row["read id"]
            if r_read_id[0] == "'":
                r_read_id = r_read_id[1:]
            if read_id[0] == "'":
                read_id = read_id[1:]
            if r_read_id not in clusters_dict:
                clusters_dict[r_read_id] = [1]
            else:
                clusters_dict[r_read_id][0] += 1
                clusters_dict[r_read_id].append(read_id)
    return clusters_dict
