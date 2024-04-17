from collections import defaultdict

import idseq_dag.util.lineage as lineage

from idseq_dag.util.dict import open_file_db_by_extension
from idseq_dag.util.m8 import NT_MIN_ALIGNMENT_LEN
from idseq_dag.util.parsing import BlastnOutput6Reader, BlastnOutput6Writer
from idseq_dag.util.parsing import BlastnOutput6NTReader
from idseq_dag.util.parsing import BlastnOutput6NTRerankedWriter
from idseq_dag.util.parsing import MAX_EVALUE_THRESHOLD


# When composing a query cover from non-overlapping fragments, consider fragments
# that overlap less than this fraction to be disjoint.
_NT_MIN_OVERLAP_FRACTION = 0.1

# Ignore NT local alignments (in blastn) with sequence similarity below 80%.
#
# Considerations:
#
#   Should not exceed the equivalent threshold for GSNAP.  Otherwise, contigs that
#   fail the higher BLAST standard would fall back on inferior results from GSNAP.
#
#   Experts agree that any NT match with less than 80% identity can be safely ignored.
#   For such highly divergent sequences, we should hope protein search finds a better
#   match than nucleotide search.
#
_NT_MIN_PIDENT = 80


def _intervals_overlap(p, q):
    '''Return True iff the intersection of p and q covers more than NT_MIN_OVERLAP_FRACTION of either p or q.'''
    p_len = p[1] - p[0]
    q_len = q[1] - q[0]
    shorter_len = min(p_len, q_len)
    intersection_len = max(0, min(p[1], q[1]) - max(p[0], q[0]))
    return intersection_len > (_NT_MIN_OVERLAP_FRACTION * shorter_len)


# Note:  HSP is a BLAST term that stands for "highest-scoring pair", i.e., a local
# alignment with no gaps that achieves one of the highest alignment scores
# in a given search.  Each HSP spans an interval in the query sequence as well
# as an interval in the subject sequence.


def _query_interval(row):
    # decode HSP query interval
    # depending on strand, qstart and qend may be reversed
    q = (row["qstart"], row["qend"])
    return min(*q), max(*q)


def _hsp_overlap(hsp_1, hsp_2):
    # Should we take into account subject_interval overlap?
    return _intervals_overlap(_query_interval(hsp_1), _query_interval(hsp_2))


def _intersects(needle, haystack):
    ''' Return True iff needle intersects haystack.  Ignore overlap < NT_MIN_OVERLAP_FRACTION. '''
    return any(_hsp_overlap(needle, hay) for hay in haystack)

def _get_lineage(accession_id, lineage_map, accession2taxid_dict, lineage_cache={}):
    if accession_id in lineage_cache:
        return lineage_cache[accession_id]
    accession_taxid = accession2taxid_dict.get(
        accession_id.split(".")[0], "NA")
    result = lineage_map.get(accession_taxid, lineage.NULL_LINEAGE)
    lineage_cache[accession_id] = result
    return result

class BlastCandidate:
    '''Documented in function get_top_m8_nt() below.'''

    def __init__(self, hsps):
        self.hsps = hsps
        self.optimal_cover = None
        self.total_score = None

        if hsps:
            # All HSPs here must be for the same query and subject sequence.
            h_0 = hsps[0]
            for h_i in self.hsps[1:]:
                assert h_i["qseqid"] == h_0["qseqid"]
                assert h_i["sseqid"] == h_0["sseqid"]
            self.qseqid = h_0["qseqid"]
            self.sseqid = h_0["sseqid"]

    def optimize(self):
        # Find a subset of disjoint HSPs with maximum sum of bitscores.
        # Initial implementation:  Super greedy.  Takes advantage of the fact
        # that blast results are sorted by bitscore, highest first.
        self.optimal_cover = [self.hsps[0]]
        for next_hsp in self.hsps[1:]:
            if not _intersects(next_hsp, self.optimal_cover):
                self.optimal_cover.append(next_hsp)
        # total_score
        self.total_score = sum(hsp["bitscore"] for hsp in self.optimal_cover)

    def summary(self):
        '''Summary stats are used later in generate_coverage_viz.'''
        r = dict(self.optimal_cover[0])
        # aggregate across the optimal cover's HSPs
        r["length"] = sum(hsp["length"] for hsp in self.optimal_cover)
        if r["length"] == 0:
            # it's never zero, but...
            r["length"] = 1
        r["pident"] = sum(hsp["pident"] * hsp["length"] for hsp in self.optimal_cover) / r["length"]
        r["bitscore"] = sum(hsp["bitscore"] for hsp in self.optimal_cover)
        # min(min(... and max(max(... below because, depending on strand, qstart and qend may be reversed
        r["qstart"] = min(min(hsp["qstart"], hsp["qend"]) for hsp in self.optimal_cover)
        r["qend"] = max(max(hsp["qstart"], hsp["qend"]) for hsp in self.optimal_cover)
        r["sstart"] = min(min(hsp["sstart"], hsp["send"]) for hsp in self.optimal_cover)
        r["send"] = max(max(hsp["sstart"], hsp["send"]) for hsp in self.optimal_cover)
        # add these two
        r["qcov"] = r["length"] / r["qlen"]
        r["hsp_count"] = len(self.optimal_cover)
        return r


def _optimal_hit_for_each_query_nr(blast_output_path, max_evalue):
    contigs_to_best_alignments = defaultdict(list)
    accession_counts = defaultdict(lambda: 0)

    with open(blast_output_path) as blastn_6_f:
        # For each contig, get the alignments that have the best total score (may be multiple if there are ties).
        for alignment in BlastnOutput6Reader(blastn_6_f):
            if alignment["evalue"] > max_evalue:
                continue
            query = alignment["qseqid"]
            best_alignments = contigs_to_best_alignments[query]

            if len(best_alignments) == 0 or best_alignments[0]["bitscore"] < alignment["bitscore"]:
                contigs_to_best_alignments[query] = [alignment]
            # Add all ties to best_hits.
            elif len(best_alignments) > 0 and best_alignments[0]["bitscore"] == alignment["bitscore"]:
                contigs_to_best_alignments[query].append(alignment)

        # Create a map of accession to best alignment count.
        for _contig_id, alignments in contigs_to_best_alignments.items():
            for alignment in alignments:
                accession_counts[alignment["sseqid"]] += 1

        # For each contig, pick the optimal alignment based on the accession that has the most best alignments.
        # If there is still a tie, arbitrarily pick the first one (later we could factor in which taxid has the most blast candidates)
        for contig_id, alignments in contigs_to_best_alignments.items():
            optimal_alignment = None
            for alignment in alignments:
                if not optimal_alignment or accession_counts[optimal_alignment["sseqid"]] < accession_counts[alignment["sseqid"]]:
                    optimal_alignment = alignment

            yield optimal_alignment


def _filter_and_group_hits_by_query(blast_output_path, min_alignment_length, min_pident, max_evalue):
    # Filter and group results by query, yielding one result group at a time.
    # A result group consists of all hits for a query, grouped by subject.
    # Please see comment in get_top_m8_nt for more context.
    current_query = None
    current_query_hits = None
    previously_seen_queries = set()
    with open(blast_output_path) as blastn_6_f:
        # Please see comments explaining the definition of "hsp" elsewhere in this file.
        for hsp in BlastnOutput6NTReader(blastn_6_f):
            # filter local alignment HSPs based on minimum length and sequence similarity
            if hsp["length"] < min_alignment_length:
                continue
            if hsp["pident"] < min_pident:
                continue
            if hsp["evalue"] > max_evalue:
                continue
            query, subject = hsp["qseqid"], hsp["sseqid"]
            if query != current_query:
                assert query not in previously_seen_queries, "blastn output appears out of order, please resort by (qseqid, sseqid, score)"
                previously_seen_queries.add(query)
                if current_query_hits:
                    yield current_query_hits
                current_query = query
                current_query_hits = defaultdict(list)
            current_query_hits[subject].append(hsp)
        if current_query_hits:
            yield current_query_hits


# An iterator that, for contig, yields to optimal hit for the contig in the nt blast_output.
def _optimal_hit_for_each_query_nt(blast_output, lineage_map, accession2taxid_dict,
                                   min_alignment_length, min_pident, max_evalue, summary=True):
    contigs_to_blast_candidates = {}
    accession_counts = defaultdict(lambda: 0)

    # For each contig, get the collection of blast candidates that have the best total score (may be multiple if there are ties).
    for query_hits in _filter_and_group_hits_by_query(blast_output, min_alignment_length, min_pident, max_evalue):
        best_hits = {}
        for _subject, hit_fragments in query_hits.items():
            # For NT, we take a special approach where we try to find a subset of disjoint HSPs
            # with maximum sum of bitscores.
            hit = BlastCandidate(hit_fragments)
            hit.optimize()

            # We prioritize the specificity of the hit; hits with species taxids are taken before hits without
            # Specificity is just the index of the tuple returned by _get_lineage(); 0 for species, 1 for genus, etc.
            lineage_taxids = _get_lineage(hit.sseqid, lineage_map, accession2taxid_dict)
            specificity = next(level for level, taxid_at_level in enumerate(lineage_taxids) if int(taxid_at_level) > 0)

            if (specificity not in best_hits) or best_hits[specificity][0].total_score < hit.total_score:
                best_hits[specificity] = [hit]
            # Add all ties to best_hits[specificity].
            elif len(best_hits[specificity]) > 0 and best_hits[specificity][0].total_score == hit.total_score:
                best_hits[specificity].append(hit)

        specific_best_hits = next(hits for specificity, hits in best_hits.items() if len(hits) > 0)
        contigs_to_blast_candidates[specific_best_hits[0].qseqid] = specific_best_hits

    # Create a map of accession to blast candidate count.
    for _contig_id, blast_candidates in contigs_to_blast_candidates.items():
        for blast_candidate in blast_candidates:
            accession_counts[blast_candidate.sseqid] += 1

    # For each contig, pick the optimal hit based on the accession that has the most total blast candidates.
    # If there is still a tie, arbitrarily pick the first one (later we could factor in which taxid has the most blast candidates)
    for contig_id, blast_candidates in contigs_to_blast_candidates.items():
        optimal_hit = None
        for blast_candidate in blast_candidates:
            if not optimal_hit or accession_counts[optimal_hit.sseqid] < accession_counts[blast_candidate.sseqid]:
                optimal_hit = blast_candidate

        if summary:
            yield optimal_hit.summary()
        else:
            for hsp in optimal_hit.optimal_cover:
                yield hsp


def get_top_m8_nr(
    blast_output,
    blast_top_blastn_6_path,
    max_evalue=MAX_EVALUE_THRESHOLD,
):
    ''' Get top m8 file entry for each contig from blast_output and output to blast_top_m8 '''
    with open(blast_top_blastn_6_path, "w") as blast_top_blastn_6_f:
        BlastnOutput6Writer(blast_top_blastn_6_f).writerows(
            _optimal_hit_for_each_query_nr(blast_output, max_evalue)
        )


def get_top_m8_nt(
    blast_output,
    lineage_map_path,
    accession2taxid_dict_path,
    blast_top_blastn_6_path,
    min_alignment_length=NT_MIN_ALIGNMENT_LEN,
    min_pident=_NT_MIN_PIDENT,
    max_evalue=MAX_EVALUE_THRESHOLD,
):
    '''
    For each contig Q (query) and reference S (subject), extend the highest-scoring
    fragment alignment of Q to S with other *non-overlapping* fragments as far as
    possible;  then output the S with the highest cumulative bitscore for Q.

    We assume cumulative bitscores can be compared across subjects for the same
    query, since that appears to be one of the ranking methods in online-blast.
    '''

    # HSP is a BLAST term that stands for "highest-scoring pair", i.e., a local
    # alignment with no gaps that achieves one of the highest alignment scores
    # in a given search.
    #
    # The output of BLAST consists of one HSP per line, where lines corresponding
    # to the same (query, subject) pair are ordered by decreasing bitscore.
    #
    # See http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    # for documentation of Blast output.

    # Output the optimal hit for each query.

    with open(blast_top_blastn_6_path, "w") as blast_top_blastn_6_f, \
        open_file_db_by_extension(lineage_map_path, "lll") as lineage_map, \
        open_file_db_by_extension(accession2taxid_dict_path, "L") as accession2taxid_dict:  # noqa
        BlastnOutput6NTRerankedWriter(blast_top_blastn_6_f).writerows(
            _optimal_hit_for_each_query_nt(blast_output, lineage_map, accession2taxid_dict, min_alignment_length, min_pident, max_evalue, False)
        )
