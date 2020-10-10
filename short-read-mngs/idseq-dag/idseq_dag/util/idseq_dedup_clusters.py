"""
idseq-dedup outputs a cluster file in the form of a csv
The first column contains the representative read id of
a cluster, and the second column contains the read id.
"""
from csv import DictReader
from typing import Dict, Optional, Tuple


def parse_clusters_file(
    idseq_dedup_clusters_path: str,
) -> Dict[str, Optional[Tuple]]:
    clusters_dict = {}
    with open(idseq_dedup_clusters_path) as f:
        for row in DictReader(f):
            r_read_id, read_id = row["representative read id"], row["read id"]
            if r_read_id not in clusters_dict:
                clusters_dict[r_read_id] = (1,)
            else:
                count, *others = clusters_dict[r_read_id]
                clusters_dict[r_read_id] = tuple([count + 1] + others + [read_id])
    return clusters_dict
