import json
import math
import os
import random
import csv

from abc import ABC
from collections import Counter
from collections import defaultdict
from typing import Any, Callable, Dict, Iterable, List, TextIO, Tuple

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
_BLAST_OUTPUT_SCHEMA = [
    ("qseqid", str),
    ("sseqid", str),
    ("pident", float),
    ("length", int),
    ("mismatch", int),
    ("gapopen", int),
    ("qstart", int),
    ("qend", int),
    ("sstart", int),
    ("send", int),
    ("evalue", float),
    ("bitscore", float),
]


# Additional blastn output columns.
_BLAST_OUTPUT_NT_SCHEMA = [
    ("qlen", int),      # query sequence length, helpful for computing qcov
    ("slen", int),      # subject sequence length, so far unused in IDseq
]


# Re-ranked output of blastn.  One row per query.  Two additional columns.
_RERANKED_BLAST_OUTPUT_NT_SCHEMA = [
    ("qcov", float),     # fraction of query covered by the optimal set of HSPs
    ("hsp_count", int),   # cardihnality of optimal fragment cover;  see BlastCandidate
]


# The minimum read count for a valid contig. We ignore contigs below this read count in most downstream analyses.
# This constant is hardcoded in at least 4 places in idseq-web.  TODO: Make it a DAG parameter.
MIN_CONTIG_SIZE = 4


class _TSVWithSchemaReader(ABC):
    """
    The _TSVWithSchemaReader class is a base class for reading TSVs without headers with schemas provided by subclasses

    This class takes a stream for reading as well as a list of partial schemas. These schemas are used to construct a
    map from the number of fields to the appropriate schema. The first partial schema becomes the first complete schema
    and each subsequent schema is appended to the previous complete schema to make a complete schema. This guaruntees that
    there is exactly one schema for each field count. This schema recognition is necessary because previously we were using
    optional fields to capture a finite set of different variants.
    """

    def __init__(self, tsv_stream: TextIO, *schemas: List[Tuple[str, Callable[[str], Any]]]) -> None:
        self._tsv_stream = tsv_stream
        assert schemas, "_TSVWithSchemaReader requires at least one schema"
        n = len(schemas[0])
        self._schema_map = {n: schemas[0]}
        for schema in schemas[1:]:
            assert schema, "_TSVWithSchemaReader does not support empty schemas"
            self._schema_map[n + len(schema)] = schemas[n] + schema
            n += len(schema)
        self._generator = self._read_all()

    def _read_all(self) -> Iterable[Dict[str, Any]]:
        """
        Parse TSV file with given schema, yielding a dict per line.
        See _BLAST_OUTPUT_SCHEMA, for example.
        When expect_headers=True, treat the first line as column headers.
        When strict mode is True, all columns in schema are required. When strict mode is False, will return None values for column not found.
        """
        for row in csv.reader(self._tsv_stream, delimiter="\t"):
            # The output of rapsearch2 contains comments that start with '#', these should be skipped
            if row and row[0][0] == "#":
                continue
            if len(row) not in self._schema_map:
                raise Exception(f"Parse error. Input line: \"{row}\" has {len(row)} columns, no associated schema found in {self._schema_map}")
            yield {
                key: _type(value) for ((key, _type), value) in zip(self._schema_map[len(row)], row)
            }

    def fields(self, columns) -> List[str]:
        return [k for k, _ in self._schema_map[columns]]

    def __iter__(self):
        return self

    def __next__(self, *args) -> Dict[str, Any]:
        return next(self._generator, *args)


class _TSVWithSchemaWriter(ABC):
    def __init__(self, tsv_stream: TextIO, *schemas: List[Tuple[str, Callable[[str], Any]]]) -> None:
        self._schema_map = {len(schema): schema for schema in schemas}
        self._writer = csv.writer(tsv_stream, delimiter="\t")

    def _dict_row_to_list(self, row: Dict[str, Any]) -> List[Any]:
        if len(row) not in self._schema_map:
            raise Exception(f"{self.path}: Write error. Input line: \"{row}\" has {len(row)} columns, no associated schema found in {self._schema_map}")
        try:
            return [row[key] for (key, _) in self._schema_map[len(row)]]
        except KeyError as e:
            raise Exception(f"{self.path}: Write error. Input line: \"{row}\" has {len(row)} columns, associated schema: {self._schema_map[len(row)]} does not have key: {e}")

    def write_all(self, rows: Iterable[Dict[str, Any]]) -> None:
        self._writer.writerows(self._dict_row_to_list(row) for row in rows)

    def write(self, row: Dict[str, Any]) -> None:
        self._writer.writerow(self._dict_row_to_list(row))

class BlastnOutput6Reader(_TSVWithSchemaReader):
    def __init__(self, tsv_stream: TextIO, filter_invalid: bool = False, min_alignment_length: int = 0):
        super().__init__(tsv_stream, _BLAST_OUTPUT_SCHEMA, _BLAST_OUTPUT_NT_SCHEMA, _RERANKED_BLAST_OUTPUT_NT_SCHEMA)
        if filter_invalid:
            self._generator = (row for row in self._generator if self.row_is_valid(row, min_alignment_length))

    @staticmethod
    def row_is_valid(row, min_alignment_length) -> bool:
        # GSNAP outputs bogus alignments (non-positive length /
        # impossible percent identity / NaN e-value) sometimes,
        # and usually they are not the only assignment, so rather than
        # killing the job, we just skip them. If we don't filter these
        # killing the job, we just skip them. If we don't filter these
        # out here, they will override the good data when computing min(
        # evalue), pollute averages computed in the json, and cause the
        # webapp loader to crash as the Rails JSON parser cannot handle
        # NaNs. Test if e_value != e_value to test if e_value is NaN
        # because NaN != NaN.
        # *** E-value Filter ***
        # Alignments with e-value > 1 are low-quality and associated with false-positives in
        # all alignments steps (NT and NR). When the e-value is greater than 1, ignore the
        # alignment
        ###
        return all([
            row["length"] > min_alignment_length,
            -0.25 < row["pident"] < 100.25,
            row["evalue"] != row["evalue"],
            row["evalue"] <= MAX_EVALUE_THRESHOLD,
        ])


class BlastnOutput6Writer(_TSVWithSchemaWriter):
    def __init__(self, tsv_stream: TextIO):
        super().__init__(tsv_stream, _BLAST_OUTPUT_SCHEMA, _BLAST_OUTPUT_NT_SCHEMA, _RERANKED_BLAST_OUTPUT_NT_SCHEMA)


_HIT_SUMMARY_SCHEMA = [
    ("read_id", str),
    ("level", int),
    ("taxid", int),
    ("accession_id", str),
    ("species_taxid", int),
    ("genus_taxid", int),
    ("family_taxid", int),
]

_HIT_SUMMARY_SCHEMA_WITH_CONTIG = [
    ("contig_id", str),
    ("contig_accession_id", str),
    ("contig_species_taxid", int),
    ("contig_genus_taxid", int),
    ("contig_family_taxid", int),
]

_HIT_SUMMARY_SCHEMA_WITH_ASSEMBLY_SOURCE = [
    ("from_assembly", str),
]

_HIT_SUMMARY_SCHEMA_MERGED = [
    ("source_count_type", str),
]

class HitSummaryReader(_TSVWithSchemaReader):
    def __init__(self, tsv_stream: TextIO) -> None:
        super().__init__(
            tsv_stream,
            _HIT_SUMMARY_SCHEMA,
            _HIT_SUMMARY_SCHEMA_WITH_CONTIG,
            _HIT_SUMMARY_SCHEMA_WITH_ASSEMBLY_SOURCE,
            _HIT_SUMMARY_SCHEMA_MERGED,
        )

class HitSummaryWriter(_TSVWithSchemaWriter):
    def __init__(self, tsv_stream: TextIO) -> None:
        super().__init__(
            tsv_stream,
            _HIT_SUMMARY_SCHEMA,
            _HIT_SUMMARY_SCHEMA_WITH_CONTIG,
            _HIT_SUMMARY_SCHEMA_WITH_ASSEMBLY_SOURCE,
            _HIT_SUMMARY_SCHEMA_MERGED,
        )


def summarize_hits(hit_summary_file_path: str, min_reads_per_genus=0):
    ''' Parse the hit summary file from alignment and get the relevant into'''
    read_dict = {}  # read_id => line
    accession_dict = {}  # accession => (species, genus)
    genus_read_counts = defaultdict(int)  # genus => read_counts
    genus_species = defaultdict(set)  # genus => list of species
    genus_accessions = defaultdict(set)  # genus => list of accessions
    total_reads = 0
    with open(hit_summary_file_path) as hit_summary_f:
        for read in HitSummaryReader(hit_summary_f):
            read_id = read["read_id"]
            accession_id = read["accession_id"]
            species_taxid = read["species_taxid"]
            genus_taxid = read["genus_taxid"]
            family_taxid = read["family_taxid"]

            read_dict[read_id] = read
            total_reads += 1
            if accession_id == "None" or accession_id == "":
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


def _call_hits_m8_work(input_blastn_6_path, lineage_map, accession2taxid_dict,
                       output_blastn_6_path, output_summary, min_alignment_length):
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
    log.write(f"Starting to summarize hits from {input_blastn_6_path}.")
    with open(input_blastn_6_path) as input_blastn_6_f:
        for row in BlastnOutput6Reader(input_blastn_6_f, filter_invalid=True, min_alignment_length=min_alignment_length):
            read_id, accession_id, e_value = row["qseqid"], row["sseqid"], row["evalue"]
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
                log.write(f"Summarized hits for {count} read ids from {input_blastn_6_path}, and counting.")

    log.write(f"Summarized hits for all {count} read ids from {input_blastn_6_path}.")

    # Generate output files. outf is the main output_m8 file and outf_sum is
    # the summary level info.
    emitted = set()
    with open(output_blastn_6_path, "w") as blastn_6_out_f, open(output_summary, "w") as outf_sum, open(input_blastn_6_path) as input_blastn_6_f:
        blastn_6_writer = BlastnOutput6Writer(blastn_6_out_f)
        # TODO: (tmorse) create sum file parser
        outf_sum_writer = csv.writer(outf_sum, delimiter="\t")
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
        for row in BlastnOutput6Reader(input_blastn_6_f, filter_invalid=True, min_alignment_length=min_alignment_length):
            read_id, accession_id, e_value = row["qseqid"], row["sseqid"], row["evalue"]
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
                blastn_6_writer.write(row)
                species_taxid = -1
                genus_taxid = -1
                family_taxid = -1
                if best_accession_id != None:
                    (species_taxid, genus_taxid, family_taxid) = get_lineage(
                        best_accession_id)

                outf_sum_writer.writerow([
                    read_id,
                    hit_level,
                    taxid,
                    best_accession_id,
                    species_taxid,
                    genus_taxid,
                    family_taxid,
                ])


@command.run_in_subprocess
def generate_taxon_count_json_from_m8(
        blastn_6_path, hit_level_file, count_type, lineage_map_path,
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
         open_file_db_by_extension(lineage_map_path) as lineage_map, \
         open(blastn_6_path) as blastn_6_f:

        # Lines in m8_file and hit_level_file correspond (same read_id)
        blastn_6_reader = BlastnOutput6Reader(blastn_6_f)
        # TODO (tmorse): make schema for hit reader
        hit_reader = csv.reader(hit_f, delimiter="\t")
        hit_row = next(hit_reader, None)
        blastn_6_row = next(blastn_6_reader, None)
        num_ranks = len(lineage.NULL_LINEAGE)
        # See https://en.wikipedia.org/wiki/Double-precision_floating-point_format
        MIN_NORMAL_POSITIVE_DOUBLE = 2.0**-1022

        with log.log_context("generate_taxon_count_json_from_m8", {"substep": "loop_1"}):
            while hit_row and blastn_6_row:
                # Retrieve data values from files
                read_id = hit_row[0]
                hit_level = hit_row[1]
                hit_taxid = hit_row[2]
                if int(hit_level) < 0:
                    log.write('int(hit_level) < 0', debug=True)
                hit_source_count_type = hit_row[13] if len(hit_row) >= 14 else None

                msg = "read_ids in %s and %s do not match: %s vs. %s" % (
                    os.path.basename(m8_file), os.path.basename(hit_level_file),
                    blastn_6_row["qseqid"], hit_row[0])
                assert blastn_6_row["qseqid"] == hit_row[0], msg
                percent_identity = blastn_6_row["pident"]
                alignment_length = blastn_6_row["length"]

                if count_type == 'merged_NT_NR' and hit_source_count_type == 'NR':
                    # NOTE: At the moment of the change, applied ONLY in the scope of the prototype of NT/NR consensus project.
                    # Protein alignments (NR) are done at amino acid level. Each amino acid is composed of 3 nucleotides.
                    # To make alignment length values comparable across NT and NR alignments (for combined statistics),
                    # the NR alignment lengths are multiplied by 3.
                    alignment_length *= 3
                e_value = blastn_6_row["evalue"]

                # These have been filtered out before the creation of m8_f and hit_f
                assert alignment_length > 0
                assert -0.25 < percent_identity < 100.25
                assert e_value == e_value

                if count_type == "NT" or hit_source_count_type == "NT":
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
                        if hit_source_count_type:
                            agg_bucket.setdefault('source_count_type', set()).add(hit_source_count_type)
                        # Chomp off the lowest rank as we aggregate up the tree
                        agg_key = agg_key[1:]

                hit_row = next(hit_reader, None)
                blastn_6_row = next(blastn_6_reader, None)

    # Produce the final output
    taxon_counts_attributes = []
    with log.log_context("generate_taxon_count_json_from_m8", {"substep": "loop_2"}):
        for agg_key, agg_bucket in aggregation.items():
            unique_count = agg_bucket['unique_count']
            nonunique_count = agg_bucket['nonunique_count']
            tax_level = num_ranks - len(agg_key) + 1
            # TODO: Extend taxonomic ranks as indicated on the commented out lines.
            taxon_counts_row = {
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
            }
            if agg_bucket.get('source_count_type'):
                taxon_counts_row['source_count_type'] = list(agg_bucket['source_count_type'])

            taxon_counts_attributes.append(taxon_counts_row)
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
