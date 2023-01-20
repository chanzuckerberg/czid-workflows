import csv
import json
from collections import defaultdict

from idseq_dag.util.parsing import HitSummaryMergedReader
from idseq_dag.util.m8 import MIN_CONTIG_SIZE, build_should_keep_filter, generate_taxon_count_json_from_m8
from idseq_dag.util.count import READ_COUNTING_MODE, ReadCountingMode, get_read_cluster_size, load_duplicate_cluster_sizes

def generate_taxon_summary(
    read2contig,
    contig2lineage,
    read_dict,
    added_reads_dict,
    db_type,
    should_keep,
    duplicate_cluster_sizes_path=None,
    use_min_contig_size=True
):
    # Return an array with
    # { taxid: , tax_level:, contig_counts: { 'contig_name': <count>, .... } }
    duplicate_cluster_sizes = None
    if duplicate_cluster_sizes_path:
        duplicate_cluster_sizes = load_duplicate_cluster_sizes(duplicate_cluster_sizes_path)

    def new_summary():
        return defaultdict(lambda: defaultdict(lambda: [0, 0]))

    genus_summary = new_summary()
    species_summary = new_summary()

    def record_read(species_taxid, genus_taxid, contig, read_id):
        if duplicate_cluster_sizes:
            cluster_size = get_read_cluster_size(duplicate_cluster_sizes, read_id)
        else:
            cluster_size = 1

        def increment(counters):
            counters[0] += 1
            counters[1] += cluster_size

        increment(species_summary[species_taxid][contig])
        increment(genus_summary[genus_taxid][contig])

    for read_id, read_info in read_dict.items():
        contig = read2contig.get(read_id, '*')
        lineage = contig2lineage.get(contig)
        if contig != '*' and lineage:
            species_taxid, genus_taxid, _family_taxid = lineage
        else:
            # not mapping to a contig, or missing contig lineage
            species_taxid, genus_taxid = read_info["species_taxid"], read_info["genus_taxid"]
            contig = '*'
        if should_keep((species_taxid, genus_taxid)):
            record_read(species_taxid, genus_taxid, contig, read_id)

    for read_id in added_reads_dict.keys():
        contig = read2contig.get(read_id, "*")
        species_taxid, genus_taxid, _family_taxid = contig2lineage[contig]
        if should_keep((species_taxid, genus_taxid)):
            record_read(species_taxid, genus_taxid, contig, read_id)

    # Filter out contigs that contain too few unique reads.
    # This used to happen in db_loader in idseq-web.  Any code left there that still appears to
    # do this filtering is effectively a no-op and the filtering cannot be done there because
    # the non-unique read counts are no longer output by the pipeline.
    for summary in [species_summary, genus_summary]:
        for taxid in list(summary.keys()):
            contig_counts = summary[taxid]
            for contig in list(contig_counts.keys()):
                unique_count, nonunique_count = contig_counts[contig]
                if unique_count < MIN_CONTIG_SIZE and use_min_contig_size:
                    del contig_counts[contig]
                else:
                    contig_counts[contig] = nonunique_count if READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL else unique_count
            if not contig_counts:
                del summary[taxid]

    # construct the array for output
    output_array = []
    for idx, summary in enumerate([species_summary, genus_summary]):
        tax_level = idx + 1
        for taxid, contig_counts in summary.items():
            entry = {'taxid': taxid, 'tax_level': tax_level,
                     'count_type': db_type.upper(), 'contig_counts': contig_counts}
            output_array.append(entry)

    return output_array


def generate_taxon_summary_from_hit_summary(
    read_to_contig_tsv_path: str,
    reassigned_m8_path: str,
    hitsummary_path: str,
    taxid_to_lineage_path: str,
    deuterostome_db_path: str,
    taxon_whitelist_path: str,
    taxon_blacklist_path: str,
    db_type: str,
    refined_counts_with_dcr_output_path: str,
    contig_summary_output_path: str,
):
    read_to_contig = {}
    contig_to_lineage = {}
    added_reads_dict = {}
    read_to_base_count = {}

    with open(read_to_contig_tsv_path) as f:
        for read_id, contig_id, sequence_length in csv.reader(f, delimiter="\t"):
            read_to_base_count[read_id] = int(sequence_length)  # reads should be theoretically unique
            read_to_contig[read_id] = contig_id

    with open(hitsummary_path) as f:
        for row in HitSummaryMergedReader(f):
            contig_to_lineage[row["contig_id"]] = (
                row["contig_species_taxid"],
                row["contig_genus_taxid"],
                row["contig_family_taxid"],
            )
            added_reads_dict[row["read_id"]] = 0

    generate_taxon_count_json_from_m8(
        reassigned_m8_path,
        hitsummary_path,
        db_type.upper(),
        taxid_to_lineage_path,
        deuterostome_db_path,
        taxon_whitelist_path,
        taxon_blacklist_path,
        duplicate_cluster_sizes_path=None,
        output_json_file=refined_counts_with_dcr_output_path,
        read_to_base_count=read_to_base_count,
    )

    should_keep = build_should_keep_filter(
        deuterostome_db_path,
        taxon_whitelist_path,
        taxon_blacklist_path,
    )

    contig_taxon_summary = generate_taxon_summary(
        read_to_contig,
        contig_to_lineage,
        {},
        added_reads_dict,
        db_type,
        should_keep,
        use_min_contig_size=False
    )

    with open(contig_summary_output_path, 'w') as contig_outf:
        json.dump(contig_taxon_summary, contig_outf)
