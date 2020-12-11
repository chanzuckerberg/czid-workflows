"""
idseq-dedup outputs a cluster file in the form of a csv
The first column contains the representative read id of
a cluster, and the second column contains the read id.
"""
from csv import DictReader
from typing import Dict, Optional, Tuple
import idseq_dag.util.log as log

def parse_clusters_file(
    idseq_dedup_clusters_path: str,
) -> Dict[str, Optional[Tuple]]:
    clusters_dict = {}
    log.write(f"opening clusters file {idseq_dedup_clusters_path}")
    with open(idseq_dedup_clusters_path) as f:
        log.write(f"opened clusters file {idseq_dedup_clusters_path}")
        for i, row in enumerate(DictReader(f)):
            if i % 100:
                log.write(f"parsed row {i}")
            r_read_id, read_id = row["representative read id"], row["read id"]
            if r_read_id not in clusters_dict:
                clusters_dict[r_read_id] = (1,)
            else:
                count, *others = clusters_dict[r_read_id]
                clusters_dict[r_read_id] = tuple([count + 1] + others + [read_id])
    return clusters_dict
