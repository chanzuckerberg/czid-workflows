import argparse
import os
import math
import json
from enum import Enum
from .call_hits import build_should_keep_filter
from ..utils.dictfile import open_file_db_by_extension
from ..utils.parsing import BlastnOutput6NTRerankedReader, HitSummaryMergedReader
from ..utils import lineage, log


class ReadCountingMode(Enum):
    COUNT_UNIQUE = "COUNT UNIQUE READS"
    COUNT_ALL = "COUNT ALL READS"


READ_COUNTING_MODE = ReadCountingMode.COUNT_ALL

################
### BEGIN CP ###
################
def generate_taxon_count_json_from_m8(
    blastn_6_path,
    hit_level_path,
    count_type,
    lineage_map_path,
    deuterostome_path,
    taxon_whitelist_path,
    taxon_blacklist_path,
    duplicate_cluster_sizes_path,
    output_json_file,
    read_to_base_count={},
):
    # Parse through hit file and m8 input file and format a JSON file with
    # our desired attributes, including aggregated statistics.

    duplicate_cluster_sizes = None
    if duplicate_cluster_sizes_path:
        duplicate_cluster_sizes = load_duplicate_cluster_sizes(
            duplicate_cluster_sizes_path
        )

    should_keep = build_should_keep_filter(
        deuterostome_path, taxon_whitelist_path, taxon_blacklist_path
    )
    # Setup
    aggregation = {}
    with open(hit_level_path) as hit_level_f, open(
        blastn_6_path
    ) as blastn_6_f, open_file_db_by_extension(lineage_map_path, "lll") as lineage_map:

        num_ranks = len(lineage.NULL_LINEAGE)
        # See https://en.wikipedia.org/wiki/Double-precision_floating-point_format
        MIN_NORMAL_POSITIVE_DOUBLE = 2.0 ** -1022

        log.log_context("generate_taxon_count_json_from_m8", {"substep": "loop_1"})
        # Lines in m8_file and hit_level_file correspond (same read_id)
        for hit_row, blastn_6_row in zip(
            HitSummaryMergedReader(hit_level_f),
            BlastnOutput6NTRerankedReader(blastn_6_f),
        ):
            # Retrieve data values from files
            read_id = hit_row["read_id"]
            hit_level = hit_row["level"]
            hit_taxid = hit_row["taxid"]
            if hit_level < 0:
                log.write("hit_level < 0", debug=True)
            hit_source_count_type = hit_row.get("source_count_type")

            msg = "read_ids in %s and %s do not match: %s vs. %s" % (
                os.path.basename(blastn_6_path),
                os.path.basename(hit_level_path),
                blastn_6_row["qseqid"],
                read_id,
            )
            assert blastn_6_row["qseqid"] == read_id, msg
            percent_identity = blastn_6_row["pident"]
            alignment_length = blastn_6_row["length"]

            if count_type == "merged_NT_NR" and hit_source_count_type == "NR":
                # NOTE: At the moment of the change, applied ONLY in the scope of the prototype of NT/NR consensus project.
                # Protein alignments (NR) are done at amino acid level. Each amino acid is composed of 3 nucleotides.
                # To make alignment length values comparable across NT and NR alignments (for combined statistics),
                # the NR alignment lengths are multiplied by 3.
                alignment_length *= 3
            e_value = blastn_6_row["evalue"]

            # These have been filtered out before the creation of blastn_6_f and hit_level_f
            assert alignment_length > 0
            assert -0.25 < percent_identity < 100.25
            assert e_value == e_value

            # e_value could be 0 when large contigs are mapped
            if e_value <= MIN_NORMAL_POSITIVE_DOUBLE:
                e_value = MIN_NORMAL_POSITIVE_DOUBLE
            e_value = math.log10(e_value)

            # Retrieve the taxon lineage and mark meaningless calls with fake
            # taxids.
            # lineage_map expects string ids
            hit_taxids_all_levels = lineage_map.get(
                str(hit_taxid), lineage.NULL_LINEAGE
            )
            cleaned_hit_taxids_all_levels = lineage.validate_taxid_lineage(
                hit_taxids_all_levels, hit_taxid, hit_level
            )
            assert num_ranks == len(cleaned_hit_taxids_all_levels)

            if should_keep(cleaned_hit_taxids_all_levels):
                # Aggregate each level and collect statistics
                agg_key = tuple(cleaned_hit_taxids_all_levels)
                while agg_key:
                    agg_bucket = aggregation.get(agg_key)
                    if not agg_bucket:
                        agg_bucket = {
                            "nonunique_count": 0,
                            "unique_count": 0,
                            "base_count": 0,
                            "sum_percent_identity": 0.0,
                            "sum_alignment_length": 0.0,
                            "sum_e_value": 0.0,
                        }
                        aggregation[agg_key] = agg_bucket
                    if duplicate_cluster_sizes:
                        agg_bucket["nonunique_count"] += get_read_cluster_size(
                            duplicate_cluster_sizes, read_id
                        )
                    else:
                        agg_bucket["nonunique_count"] += 1
                    agg_bucket["unique_count"] += 1
                    agg_bucket["base_count"] += read_to_base_count.get(read_id, 0)
                    agg_bucket["sum_percent_identity"] += percent_identity
                    agg_bucket["sum_alignment_length"] += alignment_length
                    agg_bucket["sum_e_value"] += e_value
                    if hit_source_count_type:
                        agg_bucket.setdefault("source_count_type", set()).add(
                            hit_source_count_type
                        )
                    # Chomp off the lowest rank as we aggregate up the tree
                    agg_key = agg_key[1:]

    # Produce the final output
    taxon_counts_attributes = []
    log.log_context("generate_taxon_count_json_from_m8", {"substep": "loop_2"})
    for agg_key, agg_bucket in aggregation.items():
        unique_count = agg_bucket["unique_count"]
        nonunique_count = agg_bucket["nonunique_count"]
        tax_level = num_ranks - len(agg_key) + 1
        # TODO: Extend taxonomic ranks as indicated on the commented out lines.
        taxon_counts_row = {
            "tax_id": agg_key[0],
            "tax_level": tax_level,
            # 'species_taxid' : agg_key[tax_level - 1] if tax_level == 1 else "-100",
            "genus_taxid": agg_key[2 - tax_level] if tax_level <= 2 else "-200",
            "family_taxid": agg_key[3 - tax_level] if tax_level <= 3 else "-300",
            # 'order_taxid' : agg_key[4 - tax_level] if tax_level <= 4 else "-400",
            # 'class_taxid' : agg_key[5 - tax_level] if tax_level <= 5 else "-500",
            # 'phyllum_taxid' : agg_key[6 - tax_level] if tax_level <= 6 else "-600",
            # 'kingdom_taxid' : agg_key[7 - tax_level] if tax_level <= 7 else "-700",
            # 'domain_taxid' : agg_key[8 - tax_level] if tax_level <= 8 else "-800",
            "count": nonunique_count  # this field will be consumed by the webapp
            if READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL
            else unique_count,
            "nonunique_count": nonunique_count,
            "unique_count": unique_count,
            "dcr": nonunique_count / unique_count,
            "percent_identity": agg_bucket["sum_percent_identity"] / unique_count,
            "alignment_length": agg_bucket["sum_alignment_length"] / unique_count,
            "e_value": agg_bucket["sum_e_value"] / unique_count,
            "count_type": count_type,
            "base_count": agg_bucket["base_count"],
        }
        if agg_bucket.get("source_count_type"):
            taxon_counts_row["source_count_type"] = list(
                agg_bucket["source_count_type"]
            )

        taxon_counts_attributes.append(taxon_counts_row)
    output_dict = {
        "pipeline_output": {"taxon_counts_attributes": taxon_counts_attributes}
    }

    with log.log_context(
        "generate_taxon_count_json_from_m8",
        {"substep": "json_dump", "output_json_file": output_json_file},
    ):
        with open(output_json_file, "w") as outf:
            json.dump(output_dict, outf)
            outf.flush()


##############
### END CP ###
##############


def main():
    parser = argparse.ArgumentParser(
        description="Generate taxon count JSON from m8 and hitsummary files."
    )

    parser.add_argument(
        "-f", "--m8-file", type=str, help="Path to deduped m8 file", required=True
    )
    parser.add_argument(
        "-s", "--hitsummary", type=str, help="Path to hitsummary file", required=True
    )
    parser.add_argument(
        "-c",
        "--count-type",
        type=str,
        choices=["NT", "NR", "merged_NT_NR"],
        help="Count type",
        required=True,
    )
    parser.add_argument(
        "-l", "--lineage", type=str, help="Path to lineage map file", required=True
    )
    parser.add_argument(
        "--output", type=str, default="output", help="Output prefix for JSON file"
    )
    parser.add_argument(
        "--deuterostome-path", default=None, type=str, help="Path to deuterostome file"
    )
    parser.add_argument(
        "--taxon-whitelist-path",
        default=None,
        type=str,
        help="Path to taxon whitelist file",
    )
    parser.add_argument(
        "--taxon-blacklist-path",
        default=None,
        type=str,
        help="Path to taxon blacklist file",
    )
    parser.add_argument(
        "--duplicate-cluster-sizes-path",
        default=None,
        type=str,
        help="Path to duplicate cluster sizes file",
    )
    args = parser.parse_args()

    suffix = "_with_dcr" if args.duplicate_cluster_sizes_path else ""

    generate_taxon_count_json_from_m8(
        args.m8_file,
        args.hitsummary,
        args.count_type,
        args.lineage,
        args.deuterostome_path,
        args.taxon_whitelist_path,
        args.taxon_blacklist_path,
        args.duplicate_cluster_sizes_path,
        f"{args.output}_counts{suffix}.json",
    )


if __name__ == "__main__":
    main()
