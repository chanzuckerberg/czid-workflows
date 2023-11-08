import json
from contextlib import ExitStack


def read_json(path):
    with open(path) as infile:
        return json.load(infile)


def read_truth_file(path):
    taxon_abundances = {}
    with open(path, "r") as f:
        for line in f.read().splitlines():
            taxid, _, abundance, rank, species_name = line.split("\t")[:5]
            # All current benchmark datasets have truth data on species level only,
            # so we use this as a simplifying assumption downstream
            assert rank == "species"
            taxon_abundances[taxid] = float(abundance)
    return taxon_abundances


def contigs_stats(contig_lengths, ids=None):
    lengths = []
    for id in ids if ids else contig_lengths.keys():
        lengths.append(contig_lengths[id])
    lengths.sort(reverse=True)
    total_nt = sum(lengths)
    cumlen = 0
    N50 = 0
    for i, length_i in enumerate(lengths):
        if cumlen + length_i >= total_nt / 2:
            N50 = length_i
            break
        cumlen += length_i
    return {"contigs": len(lengths), "contigs_nt": total_nt, "contigs_N50": N50}


def get_contig_stats(contig_summary_json, contig_fasta):
    contig_summary = read_json(contig_summary_json)
    contig_summary = {(elt["count_type"], elt["taxid"]): elt for elt in contig_summary}
    contigs_mapped = set()
    for elt in contig_summary.values():
        for key in elt["contig_counts"]:
            if key != "*":
                contigs_mapped.add(key)

    contig_lengths = load_contig_lengths(
        contig_fasta
    )  # TODO replace with a standard function
    return contig_summary, contigs_mapped, contig_lengths


def load_contig_lengths(contig_fasta):
    """
    Generate dict contig id -> length
    """
    lengths = {}
    with ExitStack() as cleanup:
        lines = cleanup.enter_context(open(contig_fasta))
        cur = None
        for line in lines:
            if "ASSEMBLY FAILED" in line:
                return {}

            if line.startswith(">"):
                if cur is not None:
                    lengths[cur[0]] = cur[1]
                cur = (line[1:].strip(), 0)
            else:
                cur = (cur[0], cur[1] + len(line.strip()))
        if cur is not None:
            lengths[cur[0]] = cur[1]
    return lengths


def harvest_sample_taxon_counts(
    taxon_counts_json,
    contig_summary_json,
    contig_fasta,
    dbtype,
    taxadb,
    calculate_bPM=False,
):
    assert dbtype in ("NR", "NT")

    # read in the taxon counts & contig summary JSON files
    taxon_counts = read_json(taxon_counts_json)["pipeline_output"][
        "taxon_counts_attributes"
    ]

    # read in the contig_summary and contig_lengths
    contig_summary, contig_mapped, contig_lengths = get_contig_stats(
        contig_summary_json, contig_fasta
    )

    """
    FIGURE OUT WHERE TO PUT THIS
    ans.update(contigs_stats(contig_lengths))
    for key, value in contigs_stats(contig_lengths, contigs_mapped).items():
        ans["mapped_" + key] = value
    """
    sample_reads = sum(
        [
            rslt["unique_count"]
            for rslt in taxon_counts
            if rslt["count_type"] == dbtype and rslt["tax_level"] == 1
        ]
    )
    if calculate_bPM:
        sample_bases = sum(
            [
                rslt["base_count"]
                for rslt in taxon_counts
                if rslt["count_type"] == dbtype and rslt["tax_level"] == 1
            ]
        )
    # for each species in the taxon counts (excluding genus/family for now)
    ans = {}
    for rslt in taxon_counts:
        if (
            rslt["count_type"] == dbtype
            and rslt["tax_level"] == 1
            and int(rslt["tax_id"]) >= 0
        ):
            # lookup corresponding contig summary
            contig_summary_key = (dbtype, rslt["tax_id"])
            if contig_summary_key in contig_summary:
                rslt_contigs = contig_summary[contig_summary_key]["contig_counts"]
            else:
                rslt_contigs = {}
            assert rslt["tax_id"] not in ans

            # combine info
            info = {
                "tax_level": rslt["tax_level"],
                "reads": rslt["nonunique_count"],
                "reads_dedup": rslt["unique_count"],
                "avg_aln_len": rslt["alignment_length"],
                "rPM": max(1, round(rslt["unique_count"] * 1e6 / sample_reads)),
            }

            if calculate_bPM:
                info["bPM"] = max(1, round(rslt["base_count"] * 1e6 / sample_bases))

            info.update(
                contigs_stats(contig_lengths, (k for k in rslt_contigs if k != "*"))
            )
            info["contigs_reads"] = sum(
                rslt_contigs[k] for k in rslt_contigs if k != "*"
            )
            if taxadb:
                info["tax_name"] = taxadb.sci_name(int(rslt["tax_id"]))
            ans[rslt["tax_id"]] = info

    # sort by abundance
    return {
        k: v
        for k, v in sorted(
            ans.items(),
            key=lambda kv: (kv[1]["reads_dedup"], kv[1]["avg_aln_len"]),
            reverse=True,
        )
    }
