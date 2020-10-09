import json
import math
import os
import random

from collections import Counter
from collections import defaultdict

import idseq_dag.util.command as command
import idseq_dag.util.lineage as lineage
import idseq_dag.util.log as log

from idseq_dag.util.count import READ_COUNTING_MODE, ReadCountingMode, get_read_cluster_size, load_duplicate_cluster_sizes
from idseq_dag.util.dict import open_file_db_by_extension

# NT alginments with shorter length are associated with a high rate of false positives.
# NR doesn't have this problem because Rapsearch2 contains an equivalent filter.
# Nevertheless, it may be useful to re-filter blastx results.
NT_MIN_ALIGNMENT_LEN = 36

# Alignments with e-values greater than 1 are low-quality alignments and associated with
# a high rate of false-positives. These should be filtered at all alignment steps.
MAX_EVALUE_THRESHOLD = 1

# blastn output format 6 as documented in
# http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# it's also the format of our GSNAP and RAPSEARCH2 output
BLAST_OUTPUT_SCHEMA = {
    "qseqid": str,
    "sseqid": str,
    "pident": float,
    "length": int,
    "mismatch": int,
    "gapopen": int,
    "qstart": int,
    "qend": int,
    "sstart": int,
    "send": int,
    "evalue": float,
    "bitscore": float,
}


# Additional blastn output columns.
BLAST_OUTPUT_NT_SCHEMA = dict(BLAST_OUTPUT_SCHEMA, **{
    "qlen": int,      # query sequence length, helpful for computing qcov
    "slen": int,      # subject sequence length, so far unused in IDseq
})


# Re-ranked output of blastn.  One row per query.  Two additional columns.
RERANKED_BLAST_OUTPUT_NT_SCHEMA = dict(BLAST_OUTPUT_NT_SCHEMA, **{
    "qcov": float,     # fraction of query covered by the optimal set of HSPs
    "hsp_count": int   # cardinality of optimal fragment cover;  see BlastCandidate
})


RERANKED_BLAST_OUTPUT_SCHEMA = {
    'nt': {
        'contig_level': RERANKED_BLAST_OUTPUT_NT_SCHEMA,   # only NT contigs are reranked
        'read_level': BLAST_OUTPUT_SCHEMA
    },
    'nr': {
        'contig_level': BLAST_OUTPUT_SCHEMA,
        'read_level': BLAST_OUTPUT_SCHEMA,
    }
}

# The minimum read count for a valid contig. We ignore contigs below this read count in most downstream analyses.
# This constant is hardcoded in at least 4 places in idseq-web.  TODO: Make it a DAG parameter.
MIN_CONTIG_SIZE = 4


def parse_tsv(path, schema, expect_headers=False, raw_lines=False):
    '''Parse TSV file with given schema, yielding a dict per line.  See BLAST_OUTPUT_SCHEMA, for example.  When expect_headers=True, treat the first line as column headers.'''
    assert expect_headers == False, "Headers not yet implemented."
    schema_items = schema.items()
    with open(path, "r") as stream:
        for line_number, raw_line in enumerate(stream):
            try:
                row = raw_line.rstrip("\n").split("\t")
                assert len(row) == len(schema)
                row_dict = {cname: ctype(vstr) for vstr,
                            (cname, ctype) in zip(row, schema_items)}
            except:
                msg = f"{path}:{line_number + 1}:  Parse error.  Input does not conform to schema: {schema}"
                log.write(msg, warning=True)
                raise
            if raw_lines:
                yield row_dict, raw_line
            else:
                yield row_dict


def unparse_tsv(path, schema, rows_generator):
    schema_keys = schema.keys()
    with open(path, 'w') as stream:
        for row in rows_generator:
            stream.write("\t".join(str(row[k]) for k in schema_keys) + "\n")


def log_corrupt(m8_file, line):
    msg = m8_file + " is corrupt at line:\n" + line + \
        "\n----> delete it and its corrupt ancestors before restarting run"
    log.write(msg)
    return msg


def summarize_hits(hit_summary_file, min_reads_per_genus=0):
    ''' Parse the hit summary file from alignment and get the relevant into'''
    read_dict = {}  # read_id => line
    accession_dict = {}  # accession => (species, genus)
    genus_read_counts = defaultdict(int)  # genus => read_counts
    genus_species = defaultdict(set)  # genus => list of species
    genus_accessions = defaultdict(set)  # genus => list of accessions
    total_reads = 0
    with open(hit_summary_file, 'r') as hsf:
        for line in hsf:
            read = line.rstrip().split("\t")
            accession_id, species_taxid, genus_taxid, family_taxid = read[3:7]
            read_dict[read[0]] = read
            total_reads += 1
            if accession_id == 'None' or accession_id == "":
                continue
            accession_dict[accession_id] = (
                species_taxid, genus_taxid, family_taxid)
            if int(genus_taxid) > 0:
                genus_read_counts[genus_taxid] += 1
                genus_species[genus_taxid].add(species_taxid)
                genus_accessions[genus_taxid].add(accession_id)
    selected_genera = {}  # genus => accession_list
    for genus_taxid, reads in genus_read_counts.items():
        if reads >= min_reads_per_genus and len(genus_species[genus_taxid]) > 1:
            selected_genera[genus_taxid] = list(genus_accessions[genus_taxid])

    return (read_dict, accession_dict, selected_genera)


# TODO:  Deprecate this iterate_m8() function, particularly its debug switches and modes,
# in favor of the better encapsulated and more flexible parse_tsv() above.
def iterate_m8(m8_file, min_alignment_length=0, debug_caller=None, logging_interval=25000000, full_line=False):
    """Generate an iterator over the m8 file and return values for each line.
    Work around and warn about any invalid hits detected.
    Return a subset of values (read_id, accession_id, percent_id, alignment_length,
    e_value, line) by default.
    """
    invalid_hits = 0
    last_invalid_line = None
    with open(m8_file, 'r', encoding='utf-8') as m8f:
        line_count = 0
        for line in m8f:
            line_count += 1
            if line and line[0] == '#':
                # skip comment lines
                continue
            parts = line.split("\t")
            # Must have at least 12 parts per line
            assert len(parts) >= 12, log_corrupt(m8_file, line)

            read_id = parts[0]
            accession_id = parts[1]
            percent_id = float(parts[2])
            alignment_length = int(parts[3])
            num_mismatches = int(parts[4])
            num_gaps = int(parts[5])
            query_start = int(parts[6])
            query_end = int(parts[7])
            subject_start = int(parts[8])
            subject_end = int(parts[9])
            e_value = float(parts[10])
            bitscore = float(parts[11])

            # GSNAP outputs bogus alignments (non-positive length /
            # impossible percent identity / NaN e-value) sometimes,
            # and usually they are not the only assignment, so rather than
            # killing the job, we just skip them. If we don't filter these
            # out here, they will override the good data when computing min(
            # evalue), pollute averages computed in the json, and cause the
            # webapp loader to crash as the Rails JSON parser cannot handle
            # NaNs. Test if e_value != e_value to test if e_value is NaN
            # because NaN != NaN.
            if alignment_length <= 0 or not -0.25 < percent_id < 100.25 or e_value != e_value:
                invalid_hits += 1
                last_invalid_line = line
                continue

            # *** Alignment Length Filter ***
            # Alignments with length <=35 bp are associated with false-positives in NT
            # when the .m8 file is generated via service = "gsnap" or using blast db "nt",
            # the min_alignment_length filter is passed to require a minimum .m8 length
            # note:
            #    - this is already taken care of for RAPSEARCH2 .m8 files via the -l parameter
            #    - the -l parameter used in rapsearch2 is shorter due to shorter overall alignment lengths
            #      produced by RAPSEARCH2
            ###
            if alignment_length < min_alignment_length:
                continue

            # *** E-value Filter ***
            # Alignments with e-value > 1 are low-quality and associated with false-positives in
            # all alignments steps (NT and NR). When the e-value is greater than 1, ignore the
            # alignment
            ###
            if e_value > MAX_EVALUE_THRESHOLD:
                continue

            if debug_caller and line_count % logging_interval == 0:
                msg = "Scanned {} m8 lines from {} for {}, and going.".format(
                    line_count, m8_file, debug_caller)
                log.write(msg)
            if full_line:
                yield (read_id, accession_id, percent_id, alignment_length, num_mismatches, num_gaps,
                       query_start, query_end, subject_start, subject_end, e_value, bitscore, line)
            else:
                yield (read_id, accession_id, percent_id, alignment_length,
                       e_value, bitscore, line)

    # Warn about any invalid hits outputted by GSNAP
    if invalid_hits:
        msg = "Found {} invalid hits in {};  last invalid hit line: {}".format(
            invalid_hits, m8_file, last_invalid_line)
        log.write(msg, warning=True)
    if debug_caller:
        msg = "Scanned all {} m8 lines from {} for {}.".format(
            line_count, m8_file, debug_caller)
        log.write(msg)


def read_file_into_set(file_name):
    with open(file_name, 'r') as f:
        S = set(x.rstrip() for x in f)
        S.discard('')
    return S


@command.run_in_subprocess
def call_hits_m8(input_m8, lineage_map_path, accession2taxid_dict_path,
                 output_m8, output_summary, min_alignment_length):
    """
    Determine the optimal taxon assignment for each read from the alignment
    results. When a read aligns to multiple distinct references, we need to
    assess at which level in the taxonomic hierarchy the multiple alignments
    reach consensus. We refer to this process of controlling for specificity
    as 'hit calling'.

    Input:
    - m8 file of multiple alignments per read

    Outputs:
    - cleaned m8 file with a single, optimal alignment per read
    - file with summary information, including taxonomy level at which
    specificity is reached

    Details:
    - A taxon is a group of any rank (e.g. species, genus, family, etc.).

    - A hit is a match of a read to a known reference labeled with an
    accession ID. We use NCBI's mapping of accession IDs to taxonomy IDs in
    order to retrieve the full taxonomic hierarchy for the accession ID.

    - The full taxonomy hierarchy for a hit is called its "lineage" (species,
    genus, family, etc.). A hit will normally have (positive) NCBI taxon IDs
    at all levels of the hierarchy, but there are some exceptions:

        - We use an artificial negative taxon ID if we have determined that
        the alignment is not specific at the taxonomy level under
        consideration. This happens when a read's multiple reference matches
        do not agree on taxon ID at the given level.

        For example, a read may match 5 references that all belong to
        different species (e.g. Escherichia albertii, Escherichia vulneris,
        Escherichia coli, ...), but to the same genus (Escherichia). In this
        case, we use the taxon ID for the genus (Escherichia) at the
        genus-level, but we populate the species-level with an artificial
        negative ID. The artificial ID is defined based on a negative base (
        INVALID_CALL_BASE_ID), the taxon level (e.g. 2 for genus), and the
        valid parent ID (e.g. genus Escherichia's taxon ID): see helper
        function cleaned_taxid_lineage for the precise formula.

        - Certain entries in NCBI may not have a full lineage classification;
        for example species and family will be defined but genus will be
        undefined. In this case, we populate the undefined taxonomic level
        with an artificial negative ID defined in the same manner as above.

    - m8 files correspond to BLAST tabular output format 6:
        Columns: read_id | _ref_id | percent_identity | alignment_length...

        * read_id = query (e.g., gene) sequence id
        * ref_id = subject (e.g., reference genome) sequence id
        * percent_identity = percentage of identical matches
        * alignment_length = length of the alignments
        * e_value = the expect value

        See:
        * http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
        * http://www.metagenomics.wiki/tools/blast/evalue
    """
    with open_file_db_by_extension(lineage_map_path) as lineage_map, \
         open_file_db_by_extension(accession2taxid_dict_path) as accession2taxid_dict:  # noqa
        _call_hits_m8_work(input_m8, lineage_map, accession2taxid_dict,
                           output_m8, output_summary, min_alignment_length)


def _call_hits_m8_work(input_m8, lineage_map, accession2taxid_dict,
                       output_m8, output_summary, min_alignment_length):
    lineage_cache = {}

    # Helper functions
    def get_lineage(accession_id):
        """Find the lineage of the accession ID and utilize a cache for
        performance by reducing random IOPS, ameliorating a key performance
        bottleneck
        """
        if accession_id in lineage_cache:
            return lineage_cache[accession_id]
        accession_taxid = accession2taxid_dict.get(
            accession_id.split(".")[0], "NA")
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
            accession_list = hits[level].get(
                taxid_at_level, []) + [accession_id]
            hits[level][taxid_at_level] = accession_list

    def most_frequent_accession(accession_list):
        counts = Counter(accession_list)
        return counts.most_common(1)[0][0]

    def _call_hit_level(hits):  # TODO: Remove unused function
        for level, hits_at_level in enumerate(hits):
            if len(hits_at_level) == 1:
                taxid, accession_list = hits_at_level.popitem()
                accession_id = most_frequent_accession(accession_list)
                return level + 1, taxid, accession_id
        return -1, "-1", None

    # FIXME: https://jira.czi.team/browse/IDSEQ-2738
    #  We want to move towards a general randomness solution in which
    #  all randomness is seeded based on the content of the original input.
    #  This is currently introducing non-determinism and hard coding
    #  an arbitrary seed here shouldn't impact correctness. This is only used
    #  to break ties.
    randgen = random.Random(x=4)  # chosen by fair dice role, guaranteed to be random

    def call_hit_level_v2(hits):
        ''' Always call hit at the species level with the taxid with most matches '''
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
            accession_id = most_frequent_accession(
                species_level_hits[selected_taxid])
            return 1, selected_taxid, accession_id
        return -1, "-1", None

    # Deduplicate m8 and summarize hits
    summary = {}
    count = 0
    LOG_INCREMENT = 50000
    log.write("Starting to summarize hits from {}.".format(input_m8))
    for read_id, accession_id, _percent_id, _alignment_length, e_value, _bitscore, _line in iterate_m8(
            input_m8, min_alignment_length, "call_hits_m8_scan"):
        # The Expect value (E) is a parameter that describes the number of
        # hits one can 'expect' to see by chance when searching a database of
        # a particular size. It decreases exponentially as the Score (S) of
        # the match increases. Essentially, the E value describes the random
        # background noise. https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web
        # &PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
        my_best_evalue, hits, _ = summary.get(read_id, (float("inf"), [{}, {}, {}], None))
        if my_best_evalue > e_value:
            # If we find a new better e value we want to start accumulation over
            hits = [{}, {}, {}]
            accumulate(hits, accession_id)
            my_best_evalue = e_value
        elif my_best_evalue == e_value:
            # If we find another accession with the same e value we want to accumulate it
            accumulate(hits, accession_id)
        summary[read_id] = my_best_evalue, hits, call_hit_level_v2(hits)
        count += 1
        if count % LOG_INCREMENT == 0:
            msg = "Summarized hits for {} read ids from {}, and counting.".format(
                count, input_m8)
            log.write(msg)

    log.write("Summarized hits for all {} read ids from {}.".format(
        count, input_m8))

    # Generate output files. outf is the main output_m8 file and outf_sum is
    # the summary level info.
    emitted = set()
    with open(output_m8, "w") as outf:
        with open(output_summary, "w") as outf_sum:
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
            for read_id, accession_id, _percent_id, _alignment_length, e_value, _bitscore, line in iterate_m8(
                    input_m8, min_alignment_length, "call_hits_m8_emit_deduped_and_summarized_hits"):
                if read_id in emitted:
                    continue

                # Read the fields from the summary level info
                best_e_value, _, (hit_level, taxid,
                                  best_accession_id) = summary[read_id]
                if best_e_value == e_value and best_accession_id in (
                        None, accession_id):
                    # Read out the hit with the best value that provides the
                    # most specific taxonomy information.
                    emitted.add(read_id)
                    outf.write(line)
                    species_taxid = -1
                    genus_taxid = -1
                    family_taxid = -1
                    if best_accession_id != None:
                        (species_taxid, genus_taxid, family_taxid) = get_lineage(
                            best_accession_id)

                    msg = f"{read_id}\t{hit_level}\t{taxid}\t{best_accession_id}"
                    msg += f"\t{species_taxid}\t{genus_taxid}\t{family_taxid}\n"
                    outf_sum.write(msg)


@command.run_in_subprocess
def generate_taxon_count_json_from_m8(
        m8_file, hit_level_file, e_value_type, count_type, lineage_map_path,
        deuterostome_path, taxon_whitelist_path, taxon_blacklist_path,
        duplicate_cluster_sizes_path, output_json_file):
    # Parse through hit file and m8 input file and format a JSON file with
    # our desired attributes, including aggregated statistics.

    duplicate_cluster_sizes = load_duplicate_cluster_sizes(duplicate_cluster_sizes_path)

    should_keep = build_should_keep_filter(
        deuterostome_path, taxon_whitelist_path, taxon_blacklist_path)
    # Setup
    aggregation = {}
    with open(hit_level_file, 'r', encoding='utf-8') as hit_f, \
         open(m8_file, 'r', encoding='utf-8') as m8_f, \
         open_file_db_by_extension(lineage_map_path) as lineage_map:  # noqa
        # Lines in m8_file and hit_level_file correspond (same read_id)
        hit_line = hit_f.readline()
        m8_line = m8_f.readline()
        num_ranks = len(lineage.NULL_LINEAGE)
        # See https://en.wikipedia.org/wiki/Double-precision_floating-point_format
        MIN_NORMAL_POSITIVE_DOUBLE = 2.0**-1022

        with log.log_context("generate_taxon_count_json_from_m8", {"substep": "loop_1"}):
            while hit_line and m8_line:
                # Retrieve data values from files
                hit_line_columns = hit_line.rstrip("\n").split("\t")
                read_id = hit_line_columns[0]
                hit_level = hit_line_columns[1]
                hit_taxid = hit_line_columns[2]
                if int(hit_level) < 0:
                    log.write('int(hit_level) < 0', debug=True)

                # m8 files correspond to BLAST tabular output format 6:
                # Columns: read_id | _ref_id | percent_identity | alignment_length...
                #
                # * read_id = query (e.g., gene) sequence id
                # * _ref_id = subject (e.g., reference genome) sequence id
                # * percent_identity = percentage of identical matches
                # * alignment_length = length of the alignments
                # * e_value = the expect value
                #
                # See:
                # * http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
                # * http://www.metagenomics.wiki/tools/blast/evalue

                m8_line_columns = m8_line.split("\t")
                msg = "read_ids in %s and %s do not match: %s vs. %s" % (
                    os.path.basename(m8_file), os.path.basename(hit_level_file),
                    m8_line_columns[0], hit_line_columns[0])
                assert m8_line_columns[0] == hit_line_columns[0], msg
                percent_identity = float(m8_line_columns[2])
                alignment_length = float(m8_line_columns[3])
                e_value = float(m8_line_columns[10])

                # These have been filtered out before the creation of m8_f and hit_f
                assert alignment_length > 0
                assert -0.25 < percent_identity < 100.25
                assert e_value == e_value
                if e_value_type != 'log10':
                    # e_value could be 0 when large contigs are mapped
                    if e_value <= MIN_NORMAL_POSITIVE_DOUBLE:
                        e_value = MIN_NORMAL_POSITIVE_DOUBLE
                    e_value = math.log10(e_value)

                # Retrieve the taxon lineage and mark meaningless calls with fake
                # taxids.
                hit_taxids_all_levels = lineage_map.get(
                    hit_taxid, lineage.NULL_LINEAGE)
                cleaned_hit_taxids_all_levels = lineage.validate_taxid_lineage(
                    hit_taxids_all_levels, hit_taxid, hit_level)
                assert num_ranks == len(cleaned_hit_taxids_all_levels)

                if should_keep(cleaned_hit_taxids_all_levels):
                    # Aggregate each level and collect statistics
                    agg_key = tuple(cleaned_hit_taxids_all_levels)
                    while agg_key:
                        agg_bucket = aggregation.get(agg_key)
                        if not agg_bucket:
                            agg_bucket = {
                                'nonunique_count': 0,
                                'unique_count': 0,
                                'sum_percent_identity': 0.0,
                                'sum_alignment_length': 0.0,
                                'sum_e_value': 0.0
                            }
                            aggregation[agg_key] = agg_bucket
                        agg_bucket['nonunique_count'] += get_read_cluster_size(
                            duplicate_cluster_sizes, read_id)
                        agg_bucket['unique_count'] += 1
                        agg_bucket['sum_percent_identity'] += percent_identity
                        agg_bucket['sum_alignment_length'] += alignment_length
                        agg_bucket['sum_e_value'] += e_value
                        # Chomp off the lowest rank as we aggregate up the tree
                        agg_key = agg_key[1:]

                hit_line = hit_f.readline()
                m8_line = m8_f.readline()

    # Produce the final output
    taxon_counts_attributes = []
    with log.log_context("generate_taxon_count_json_from_m8", {"substep": "loop_2"}):
        for agg_key, agg_bucket in aggregation.items():
            unique_count = agg_bucket['unique_count']
            nonunique_count = agg_bucket['nonunique_count']
            tax_level = num_ranks - len(agg_key) + 1
            # TODO: Extend taxonomic ranks as indicated on the commented out lines.
            taxon_counts_attributes.append({
                "tax_id":
                agg_key[0],
                "tax_level":
                tax_level,
                # 'species_taxid' : agg_key[tax_level - 1] if tax_level == 1 else "-100",
                'genus_taxid':
                agg_key[2 - tax_level] if tax_level <= 2 else "-200",
                'family_taxid':
                agg_key[3 - tax_level] if tax_level <= 3 else "-300",
                # 'order_taxid' : agg_key[4 - tax_level] if tax_level <= 4 else "-400",
                # 'class_taxid' : agg_key[5 - tax_level] if tax_level <= 5 else "-500",
                # 'phyllum_taxid' : agg_key[6 - tax_level] if tax_level <= 6 else "-600",
                # 'kingdom_taxid' : agg_key[7 - tax_level] if tax_level <= 7 else "-700",
                # 'domain_taxid' : agg_key[8 - tax_level] if tax_level <= 8 else "-800",
                "count":  # this field will be consumed by the webapp
                nonunique_count if READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL else unique_count,
                "nonunique_count":
                nonunique_count,
                "unique_count":
                unique_count,
                "dcr":
                nonunique_count / unique_count,
                "percent_identity":
                agg_bucket['sum_percent_identity'] / unique_count,
                "alignment_length":
                agg_bucket['sum_alignment_length'] / unique_count,
                "e_value":
                agg_bucket['sum_e_value'] / unique_count,
                "count_type":
                count_type
            })
        output_dict = {
            "pipeline_output": {
                "taxon_counts_attributes": taxon_counts_attributes
            }
        }

    with log.log_context(
        "generate_taxon_count_json_from_m8",
        {"substep": "json_dump", "output_json_file": output_json_file}
    ):
        with open(output_json_file, 'w') as outf:
            json.dump(output_dict, outf)
            outf.flush()


def build_should_keep_filter(
    deuterostome_path,
    taxon_whitelist_path,
    taxon_blacklist_path
):

    # See also HOMO_SAPIENS_TAX_IDS in idseq-web
    taxids_to_remove = set(['9605', '9606'])

    if taxon_blacklist_path:
        with log.log_context("generate_taxon_count_json_from_m8", {"substep": "read_blacklist_into_set"}):
            taxids_to_remove.update(read_file_into_set(taxon_blacklist_path))

    if deuterostome_path:
        with log.log_context("generate_taxon_count_json_from_m8", {"substep": "read_file_into_set"}):
            taxids_to_remove.update(read_file_into_set(deuterostome_path))

    if taxon_whitelist_path:
        with log.log_context("generate_taxon_count_json_from_m8", {"substep": "read_whitelist_into_set"}):
            taxids_to_keep = read_file_into_set(taxon_whitelist_path)

    def is_blacklisted(hits):
        for taxid in hits:
            if int(taxid) >= 0 and taxid in taxids_to_remove:
                return True
        return False

    def is_whitelisted(hits):
        if not taxon_whitelist_path:
            return True
        for taxid in hits:
            if int(taxid) >= 0 and taxid in taxids_to_keep:
                return True
        return False

    def should_keep(hits):
        return is_whitelisted(hits) and not is_blacklisted(hits)

    return should_keep
