import argparse
from typing import Iterable
import random
import logging
from collections import Counter
from ..utils.dictfile import open_file_db_by_extension
from ..utils.parsing import BlastnOutput6Reader, BlastnOutput6Writer, HitSummaryWriter

log = logging.getLogger(__name__)
log.log_context = lambda x, y: print(x, y)  # monkey patch log function
log.write = lambda x: print(x)  # monkey patch log function

# mock lineage


class lineage:
    NULL_SPECIES_ID = "-100"
    NULL_GENUS_ID = "-200"
    NULL_FAMILY_ID = "-300"
    NULL_LINEAGE = (NULL_SPECIES_ID, NULL_GENUS_ID, NULL_FAMILY_ID)


################
### CP PASTE ###
################
def build_should_keep_filter(
    deuterostome_path, taxon_whitelist_path, taxon_blacklist_path
):

    # See also HOMO_SAPIENS_TAX_IDS in idseq-web
    taxids_to_remove = set(["9605", "9606"])

    if taxon_blacklist_path:
        with log.log_context(
            "generate_taxon_count_json_from_m8", {"substep": "read_blacklist_into_set"}
        ):
            taxids_to_remove.update(read_file_into_set(taxon_blacklist_path))

    if deuterostome_path:
        with log.log_context(
            "generate_taxon_count_json_from_m8", {"substep": "read_file_into_set"}
        ):
            taxids_to_remove.update(read_file_into_set(deuterostome_path))

    if taxon_whitelist_path:
        with log.log_context(
            "generate_taxon_count_json_from_m8", {"substep": "read_whitelist_into_set"}
        ):
            taxids_to_keep = read_file_into_set(taxon_whitelist_path)

    def is_blacklisted(hits: Iterable[str]):
        for taxid in hits:
            if int(taxid) >= 0 and taxid in taxids_to_remove:
                return True
        return False

    def is_whitelisted(hits: Iterable[str]):
        if not taxon_whitelist_path:
            return True
        for taxid in hits:
            if int(taxid) >= 0 and taxid in taxids_to_keep:
                return True
        return False

    def should_keep(hits: Iterable[str]):
        # In some places in the code taxids are ints rather than strings, this would lead
        # to a silent failure here so it is worth the explicit check.
        non_strings = [h for h in hits if not isinstance(h, str)]
        assert not non_strings, f"should_keep recieved non-string inputs {non_strings}"
        return is_whitelisted(hits) and not is_blacklisted(hits)

    return should_keep


def _call_hits_m8_work(
    input_blastn_6_path,
    lineage_map,
    accession2taxid_dict,
    output_blastn_6_path,
    output_summary,
    min_alignment_length,
    deuterostome_path,
    taxon_whitelist_path,
    taxon_blacklist_path,
):
    lineage_cache = {}

    should_keep = build_should_keep_filter(
        deuterostome_path, taxon_whitelist_path, taxon_blacklist_path
    )

    # Helper functions
    def get_lineage(accession_id):
        """Find the lineage of the accession ID and utilize a cache for
        performance by reducing random IOPS, ameliorating a key performance
        bottleneck
        """
        if accession_id in lineage_cache:
            return lineage_cache[accession_id]
        accession_taxid = accession2taxid_dict.get(accession_id.split(".")[0], "NA")
        result = lineage_map.get(accession_taxid, lineage.NULL_LINEAGE)
        lineage_cache[accession_id] = result
        return result

    def accumulate(hits, accession_id):
        """Accumulate hits for summarizing hit information and specificity at
        each taxonomy level.
        """
        lineage_taxids = get_lineage(accession_id)
        for level, taxid_at_level in enumerate(lineage_taxids):
            if int(taxid_at_level) < 0:
                # Skip if we have a negative taxid. When an accession doesn't
                # provide species level info, it doesn't contradict any info
                # provided by other accessions. This occurs a lot and
                # handling it in this way seems to work well.
                continue
            accession_list = hits[level].get(taxid_at_level, []) + [accession_id]
            hits[level][taxid_at_level] = accession_list

    def most_frequent_accession(accession_list):
        counts = Counter(accession_list)
        return counts.most_common(1)[0][0]

    # FIXME: https://jira.czi.team/browse/IDSEQ-2738
    #  We want to move towards a general randomness solution in which
    #  all randomness is seeded based on the content of the original input.
    #  This is currently introducing non-determinism and hard coding
    #  an arbitrary seed here shouldn't impact correctness. This is only used
    #  to break ties.
    randgen = random.Random(x=4)  # chosen by fair dice role, guaranteed to be random

    def call_hit_level_v2(hits):
        """ Always call hit at the species level with the taxid with most matches """
        species_level_hits = hits[0]
        max_match = 0
        taxid_candidates = []
        for taxid, accession_list in species_level_hits.items():
            accession_len = len(accession_list)
            if accession_len > max_match:
                taxid_candidates = [taxid]
                max_match = accession_len
            elif accession_len == max_match:
                taxid_candidates.append(taxid)
        if max_match > 0:
            selected_taxid = taxid_candidates[0]
            if len(taxid_candidates) > 1:
                selected_taxid = randgen.sample(taxid_candidates, 1)[0]
            accession_id = most_frequent_accession(species_level_hits[selected_taxid])
            return 1, selected_taxid, accession_id
        return -1, "-1", None

    # Deduplicate m8 and summarize hits
    summary = {}
    count = 0
    LOG_INCREMENT = 50000
    log.write(f"Starting to summarize hits from {input_blastn_6_path}.")
    with open(input_blastn_6_path) as input_blastn_6_f:
        for row in BlastnOutput6Reader(
            input_blastn_6_f,
            filter_invalid=True,
            min_alignment_length=min_alignment_length,
        ):
            read_id, accession_id, bitscore = (
                row["qseqid"],
                row["sseqid"],
                row["bitscore"],
            )
            # The Expect value (E) is a parameter that describes the number of
            # hits one can 'expect' to see by chance when searching a database of
            # a particular size. It decreases exponentially as the Score (S) of
            # the match increases. Essentially, the E value describes the random
            # background noise. https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web
            # &PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
            # We have since moved to using the bitscore rather than the e-value
            my_best_bitscore, hits, _ = summary.get(
                read_id, (float("-inf"), [{}, {}, {}], None)
            )
            if my_best_bitscore < bitscore:
                # If we find a new better bitscore we want to start accumulation over
                hits = [{}, {}, {}]
                accumulate(hits, accession_id)
                my_best_bitscore = bitscore
            elif my_best_bitscore == bitscore:
                # If we find another accession with the same bitscore we want to accumulate it
                accumulate(hits, accession_id)
            summary[read_id] = my_best_bitscore, hits, call_hit_level_v2(hits)
            count += 1
            if count % LOG_INCREMENT == 0:
                log.write(
                    f"Summarized hits for {count} read ids from {input_blastn_6_path}, and counting."
                )

    log.write(f"Summarized hits for all {count} read ids from {input_blastn_6_path}.")

    # Generate output files. outf is the main output_m8 file and outf_sum is
    # the summary level info.
    emitted = set()
    with open(output_blastn_6_path, "w") as blastn_6_out_f, open(
        output_summary, "w"
    ) as hit_summary_out_f, open(input_blastn_6_path) as input_blastn_6_f:
        blastn_6_writer = BlastnOutput6Writer(blastn_6_out_f)
        hit_summary_writer = HitSummaryWriter(hit_summary_out_f)
        # Iterator over the lines of the m8 file. Emit the hit with the
        # best value that provides the most specific taxonomy
        # information. If there are multiple hits (also called multiple
        # accession IDs) for a given read that all have the same e-value,
        # some may provide species information and some may only provide
        # genus information. We want to emit the one that provides the
        # species information because from that we can infer the rest of
        # the lineage. If we accidentally emitted the one that provided
        # only genus info, downstream steps may have difficulty
        # recovering the species.

        # TODO: Consider all hits within a fixed margin of the best e-value.
        # This change may need to be accompanied by a change to
        # GSNAP/RAPSearch2 parameters.
        for row in BlastnOutput6Reader(
            input_blastn_6_f,
            filter_invalid=True,
            min_alignment_length=min_alignment_length,
        ):
            read_id, accession_id, bitscore = (
                row["qseqid"],
                row["sseqid"],
                row["bitscore"],
            )
            if read_id in emitted:
                continue

            # Read the fields from the summary level info
            best_bitscore, _, (hit_level, taxid, best_accession_id) = summary[read_id]
            if (
                best_bitscore == bitscore
                and best_accession_id in (None, accession_id)
                and should_keep([taxid])
            ):
                # Read out the hit with the best value that provides the
                # most specific taxonomy information.
                emitted.add(read_id)
                blastn_6_writer.writerow(row)
                species_taxid = -1
                genus_taxid = -1
                family_taxid = -1
                if best_accession_id != None:
                    (species_taxid, genus_taxid, family_taxid) = get_lineage(
                        best_accession_id
                    )

                hit_summary_writer.writerow(
                    {
                        "read_id": read_id,
                        "level": hit_level,
                        "taxid": taxid,
                        "accession_id": best_accession_id,
                        "species_taxid": species_taxid,
                        "genus_taxid": genus_taxid,
                        "family_taxid": family_taxid,
                    }
                )


##############
### END CP ###
##############

# ---------------------------------------------------------------------------------#
def call_hits_m8(
    input_m8,
    lineage_map_path,
    accession2taxid_dict_path,
    min_alignment_length,
    deuterostome_path,
    taxon_whitelist_path,
    taxon_blacklist_path,
    output_prefix,
):
    with open_file_db_by_extension(
        lineage_map_path, "lll"
    ) as lineage_map, open_file_db_by_extension(
        accession2taxid_dict_path, "L"
    ) as accession2taxid_dict:
        # With lineage_map and accession2taxid_dict open

        _call_hits_m8_work(
            input_m8,
            lineage_map,
            accession2taxid_dict,
            output_blastn_6_path=f"{output_prefix}.deduped.m8",
            output_summary=f"{output_prefix}.hitsummary.tab",
            min_alignment_length=min_alignment_length,
            deuterostome_path=deuterostome_path,
            taxon_whitelist_path=taxon_whitelist_path,
            taxon_blacklist_path=taxon_blacklist_path,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Given an m8 file, call hits for each query"
    )
    parser.add_argument("-f", "--file", help="Path to the input m8 file", required=True)
    parser.add_argument(
        "-l", "--lineage", help="Path to the lineage file", required=True
    )
    parser.add_argument(
        "-a",
        "--accession2taxid",
        help="Path to the accession2taxid file",
        required=True,
    )
    parser.add_argument("-o", "--output", help="Output prefix", required=True)

    parser.add_argument(
        "--min-alignment-length",
        type=int,
        default=0,
        help="Minimum alignment length to consider",
    )
    parser.add_argument(
        "--deuterostome-path", default=None, help="Path to the deuterostome file"
    )
    parser.add_argument(
        "--taxon-whitelist-path", default=None, help="Path to the taxon whitelist"
    )
    parser.add_argument(
        "--taxon-blacklist-path", default=None, help="Path to the taxon blacklist"
    )

    args = parser.parse_args()

    call_hits_m8(
        args.file,
        args.lineage,
        args.accession2taxid,
        args.min_alignment_length,
        args.deuterostome_path,
        args.taxon_whitelist_path,
        args.taxon_blacklist_path,
        args.output,
    )


if __name__ == "__main__":
    main()
