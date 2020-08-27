"""
Example input lines that form a cluster:

   "0       140nt, >M05295:357:000000000-CRPNR:1:2119:16143:8253... *"
   "1       140nt, >M05295:357:000000000-CRPNR:1:1101:22051:10534... at 1:140:1:140/+/100.00%"
   "2       140nt, >M05295:357:000000000-CRPNR:1:1102:15401:7483... at 1:140:1:140/+/100.00%"
   ...
   "2334    140nt, >M05295:357:000000000-CRPNR:1:1102:13405:3483... at 1:140:1:140/+/100.00%"

Corresponding output line for that cluster:

   "2335    140nt, >M05295:357:000000000-CRPNR:1:2119:16143:8253"

Please note that "..." above does not indicate truncation. CD-HIT-DUP appends "..." to the read
IDs even if the read IDs have not been truncated.

Per CD-HIT-DUP docs, when a "-d" argument is not specified, each read ID is obtained by
splitting the corresponding FASTA line on whitespace and taking the first item.  That is also
how we do it throughout this pipeline.  The code below assumes (and relies upon) the read IDs
being free from whitespace.

Furthremore, we expect read IDs to be unique in a sequencing run.  Violating that assumption,
if it does not break cdhit itself, might produce slightly bogus numbers, because of how it
affects subsampling and per-read DCR correction.  Preserving this uniqueness is why we do not
specify a "-d" flag to cdhit, and allow it thus to use the entire read id.

Further note that, although CD-HIT-DUP documentation states the read marked with '*' is the
one chosen as representative, we have found, particularly with unpaired fasta, the emitted
deduplicated read is not the one marked with '*'.  Hence we handle that case correctly below,
and put a lot of assertions to make sure problems with cdhit output are detected.
"""

import idseq_dag.util.log as log
from idseq_dag.util.fasta import iterator
from typing import Dict, Optional, Set, Tuple


def parse_clusters_file(
    cdhit_clusters_path: str,
    deduped_fasta_path: str,
) -> Dict[str, Optional[Tuple]]:
    # First identify the cluster representative reads emitted by cd-hit-dup.  Originally we
    # used the ".clstr" output of cd-hit-dup for this, but turns out that for unpaired reads
    # the actual deduped output of cdhit contains different representatives.
    clusters_dict: Dict[str, Optional[Tuple]] = {}
    for read in iterator(deduped_fasta_path):
        # A read header looks someting like
        #
        #    >M05295:357:000000000-CRPNR:1:1101:22051:10534 OPTIONAL RANDOM STUFF"
        #     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #
        # where the first character on the line is '>' and the read ID (underlined above)
        # extends from '>' to the first whitespace character, not including '>' itself.
        #
        # The fasta iterator already asserts that read.header[0] is '>'.
        #
        read_id = read.header.split(None, 1)[0][1:]
        clusters_dict[read_id] = None  # not yet known

    def record_cluster_size(
        cluster_size: int,
        emitted_reads_from_cluster: set,
        other_reads_from_cluster: set,
        line_number: int,
    ):
        assert emitted_reads_from_cluster, f"""If this assertion fails,
        CD-HIT-DUP has forgotten to emit a read for this cluster.  In that case,
        just use the current read_id as cluster_representative.  Everything will
        work fine, aside from reduced sensitivity. {line_number}"""

        if len(emitted_reads_from_cluster) != 1:
            """
            If this check is true, cd-hit-dup has emitted multiple reads
            from the same cluster. We have observed that this does happen.
            Everything will run fine, but read counts contributed by that
            cluster will be exaggerated. We are in the process of replacing
            cd-hit-dup with a tool that does not have this bug but until
            that time we will just output a warning so we don't block results.
            """
            log.write(f"more than one read emitted from cd-hit-dup cluster: {emitted_reads_from_cluster} reads emitted on line {line_number}", warning=True)

        cluster_representative = emitted_reads_from_cluster.pop()

        assert cluster_representative in clusters_dict, "If this fails it's our bug here."

        assert cluster_size - 1 == len(other_reads_from_cluster), """other_reads_from_cluster should
        contain the number of reads specified by cluster_size minus cluster_representative:
        {}, {}""".format(cluster_size, other_reads_from_cluster)

        clusters_dict[cluster_representative] = (cluster_size,) + tuple(other_reads_from_cluster)
        return

    with open(cdhit_clusters_path, "r") as clusters_file:
        # set of reads in both dedup1.fa and current cluster; cardinality 1!
        emitted_reads_from_cluster: Set[str] = set()
        other_reads_from_cluster: Set[str] = set()
        cluster_size = 0
        read_id = None
        line_number = 0
        for line in clusters_file:
            line_number += 1
            if line.startswith(">"):
                continue
            parts = line.strip().split()
            serial = int(parts[0])
            assert parts[2][0] == ">", line
            assert parts[2].endswith("..."), line
            if serial == 0 and cluster_size > 0:
                # We've just encountered the first read of a new cluster.  Emit
                # all data held for the old cluster.
                record_cluster_size(
                    cluster_size,
                    emitted_reads_from_cluster,
                    other_reads_from_cluster,
                    line_number,
                )
                emitted_reads_from_cluster = set()
                other_reads_from_cluster = set()
                cluster_size = 0
            assert cluster_size == serial, f"{line_number}: {cluster_size}, {serial}, {line}"
            read_id = parts[2][1:-3]
            cluster_size += 1
            if read_id in clusters_dict:
                emitted_reads_from_cluster.add(read_id)
            else:
                other_reads_from_cluster.add(read_id)
        # record last cluster
        if cluster_size > 0:
            record_cluster_size(
                cluster_size,
                emitted_reads_from_cluster,
                other_reads_from_cluster,
                line_number,
            )

    return clusters_dict
