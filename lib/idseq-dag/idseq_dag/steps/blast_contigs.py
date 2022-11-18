import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log
import idseq_dag.util.m8 as m8
import idseq_dag.util.s3 as s3
import csv
import json
import multiprocessing
import os

from collections import defaultdict
from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.trace_lock import TraceLock

from idseq_dag.steps.run_assembly import PipelineStepRunAssembly
from idseq_dag.util.count import READ_COUNTING_MODE, ReadCountingMode, get_read_cluster_size, load_duplicate_cluster_sizes
from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3, DEFAULT_WHITELIST_S3
from idseq_dag.util.m8 import MIN_CONTIG_SIZE, NT_MIN_ALIGNMENT_LEN
from idseq_dag.util.parsing import BlastnOutput6Reader, BlastnOutput6Writer
from idseq_dag.util.parsing import BlastnOutput6NTReader
from idseq_dag.util.parsing import BlastnOutput6NTRerankedReader, BlastnOutput6NTRerankedWriter
from idseq_dag.util.parsing import HitSummaryReader, HitSummaryMergedWriter
from idseq_dag.util.parsing import MAX_EVALUE_THRESHOLD


MIN_REF_FASTA_SIZE = 25
MIN_ASSEMBLED_CONTIG_SIZE = 25

# When composing a query cover from non-overlapping fragments, consider fragments
# that overlap less than this fraction to be disjoint.
NT_MIN_OVERLAP_FRACTION = 0.1

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
NT_MIN_PIDENT = 80


def intervals_overlap(p, q):
    '''Return True iff the intersection of p and q covers more than NT_MIN_OVERLAP_FRACTION of either p or q.'''
    p_len = p[1] - p[0]
    q_len = q[1] - q[0]
    shorter_len = min(p_len, q_len)
    intersection_len = max(0, min(p[1], q[1]) - max(p[0], q[0]))
    return intersection_len > (NT_MIN_OVERLAP_FRACTION * shorter_len)


# Note:  HSP is a BLAST term that stands for "highest-scoring pair", i.e., a local
# alignment with no gaps that achieves one of the highest alignment scores
# in a given search.  Each HSP spans an interval in the query sequence as well
# as an interval in the subject sequence.


def query_interval(row):
    # decode HSP query interval
    # depending on strand, qstart and qend may be reversed
    q = (row["qstart"], row["qend"])
    return min(*q), max(*q)


def hsp_overlap(hsp_1, hsp_2):
    # Should we take into account subject_interval overlap?
    return intervals_overlap(query_interval(hsp_1), query_interval(hsp_2))


def intersects(needle, haystack):
    ''' Return True iff needle intersects haystack.  Ignore overlap < NT_MIN_OVERLAP_FRACTION. '''
    return any(hsp_overlap(needle, hay) for hay in haystack)


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
            if not intersects(next_hsp, self.optimal_cover):
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


def _optimal_hit_for_each_query_nt(blast_output, min_alignment_length, min_pident, max_evalue, summary=True):
    contigs_to_blast_candidates = {}
    accession_counts = defaultdict(lambda: 0)

    # For each contig, get the collection of blast candidates that have the best total score (may be multiple if there are ties).
    for query_hits in _filter_and_group_hits_by_query(blast_output, min_alignment_length, min_pident, max_evalue):
        best_hits = []
        for _subject, hit_fragments in query_hits.items():
            # For NT, we take a special approach where we try to find a subset of disjoint HSPs
            # with maximum sum of bitscores.
            hit = BlastCandidate(hit_fragments)
            hit.optimize()

            if len(best_hits) == 0 or best_hits[0].total_score < hit.total_score:
                best_hits = [hit]
            # Add all ties to best_hits.
            elif len(best_hits) > 0 and best_hits[0].total_score == hit.total_score:
                best_hits.append(hit)

        contigs_to_blast_candidates[best_hits[0].qseqid] = best_hits

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


def _get_top_m8_nt(
    blast_output,
    blast_top_blastn_6_path,
    min_alignment_length=NT_MIN_ALIGNMENT_LEN,
    min_pident=NT_MIN_PIDENT,
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

    with open(blast_top_blastn_6_path, "w") as blast_top_blastn_6_f:
        BlastnOutput6NTRerankedWriter(blast_top_blastn_6_f).writerows(
            _optimal_hit_for_each_query_nt(
                blast_output,
                min_alignment_length,
                min_pident,
                max_evalue,
                False,
            )
        )

def _get_top_m8_nr(
    blast_output,
    blast_top_blastn_6_path,
    max_evalue=MAX_EVALUE_THRESHOLD,
):
    ''' Get top m8 file entry for each contig from blast_output and output to blast_top_m8'''
    with open(blast_top_blastn_6_path, "w") as blast_top_blastn_6_f:
        BlastnOutput6Writer(blast_top_blastn_6_f).writerows(
            _optimal_hit_for_each_query_nr(blast_output, max_evalue)
        )


def _update_read_dict(read2contig, blast_top_blastn_6_path, read_dict, accession_dict, db_type):
    consolidated_dict = read_dict
    read2blastm8 = {}
    contig2accession = {}
    contig2lineage = {}
    added_reads = {}

    with open(blast_top_blastn_6_path) as blast_top_blastn_6_f:
        blastn_6_reader = BlastnOutput6NTRerankedReader(blast_top_blastn_6_f) if db_type == 'nt' else BlastnOutput6Reader(blast_top_blastn_6_f)
        for row in blastn_6_reader:
            contig_id = row["qseqid"]
            accession_id = row["sseqid"]
            if contig2accession.get(contig_id, None):
                # if contig already exists, synthesize subcontigs
                curr = contig2accession[contig_id][1]

                # take the weighted mean
                curr["pident"] = (curr["pident"] * curr["length"] + row["pident"] * row["length"]) / (curr["length"] + row["length"])
                curr["evalue"] = (curr["evalue"] * curr["length"] + row["evalue"] * row["length"]) / (curr["length"] + row["length"])

                # take the sum for these values
                curr["mismatch"] += row["mismatch"]
                curr["gapopen"] += row["gapopen"]
                curr["bitscore"] += row["bitscore"]
                curr["length"] += row["length"]

                # take the min/max for these values
                curr["qstart"] = min(curr["qstart"], row["qstart"])
                curr["qend"] = max(curr["qend"], row["qend"])
                curr["sstart"] = min(curr["sstart"], row["sstart"])
                curr["send"] = max(curr["send"], row["send"])
            else:
                # else add the contig
                contig2accession[contig_id] = (accession_id, row)
            if accession_id in accession_dict:
                contig2lineage[contig_id] = accession_dict[accession_id]

        for read_id, contig_id in read2contig.items():
            (accession, m8_row) = contig2accession.get(contig_id, (None, None))
            # accession_dict comes from hit_summary, which comes from alignment and is filtered for taxids
            # this means that we don't need to filter here because we will never get unfiltered taxids from
            # accession_dict, however it may be missing accessions so we must handle that case.
            if accession and accession in accession_dict:
                (species_taxid, genus_taxid, family_taxid) = accession_dict[accession]
                if read_id in consolidated_dict:
                    consolidated_dict[read_id]["taxid"] = species_taxid
                    consolidated_dict[read_id]["contig_id"] = contig_id
                    consolidated_dict[read_id]["contig_accession_id"] = accession
                    consolidated_dict[read_id]["contig_species_taxid"] = species_taxid
                    consolidated_dict[read_id]["contig_genus_taxid"] = genus_taxid
                    consolidated_dict[read_id]["contig_family_taxid"] = family_taxid
                else:
                    added_reads[read_id] = {
                        "read_id": read_id,
                        "level": 1,
                        "taxid": species_taxid,
                        "accession_id": accession,
                        "species_taxid": species_taxid,
                        "genus_taxid": genus_taxid,
                        "family_taxid": family_taxid,
                        "contig_id": contig_id,
                        "contig_accession_id": accession,
                        "contig_species_taxid": species_taxid,
                        "contig_genus_taxid": genus_taxid,
                        "contig_family_taxid": family_taxid,
                        "from_assembly": "from_assembly",
                    }
            if m8_row:
                read2blastm8[read_id] = m8_row
        return (consolidated_dict, read2blastm8, contig2lineage, added_reads)


def _generate_taxon_summary(
    read2contig,
    contig2lineage,
    read_dict,
    added_reads_dict,
    db_type,
    duplicate_cluster_sizes_path,
    should_keep,
    read2contig_bases={},
):
    # Return an array with
    # { taxid: , tax_level:, contig_counts: { 'contig_name': <count>, .... } }
    duplicate_cluster_sizes = None
    if duplicate_cluster_sizes_path:
        duplicate_cluster_sizes = load_duplicate_cluster_sizes(duplicate_cluster_sizes_path)

    def new_summary():
        return defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))

    genus_summary = new_summary()
    species_summary = new_summary()

    def record_read(species_taxid, genus_taxid, contig, read_id, base_count=None):
        cluster_size = get_read_cluster_size(duplicate_cluster_sizes, read_id) if duplicate_cluster_sizes else 1

        def increment(counters):
            counters[0] += 1
            counters[1] += cluster_size
            if base_count:
                counters[2] += base_count

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
            record_read(species_taxid, genus_taxid, contig, read_id, read2contig_bases.get((read_id, contig), None))

    for read_id in added_reads_dict.keys():
        contig = read2contig[read_id]
        species_taxid, genus_taxid, _family_taxid = contig2lineage[contig]
        if should_keep((species_taxid, genus_taxid)):
            record_read(species_taxid, genus_taxid, contig, read_id, read2contig_bases.get((read_id, contig), None))

    # Filter out contigs that contain too few unique reads.
    # This used to happen in db_loader in idseq-web.  Any code left there that still appears to
    # do this filtering is effectively a no-op and the filtering cannot be done there because
    # the non-unique read counts are no longer output by the pipeline.
    for summary in [species_summary, genus_summary]:
        for taxid in list(summary.keys()):
            contig_counts = summary[taxid]
            for contig in list(contig_counts.keys()):
                unique_count, nonunique_count, _ = contig_counts[contig]
                if unique_count < MIN_CONTIG_SIZE:
                    del contig_counts[contig]
                elif not read2contig_bases:  # if there are bases to count, output everything
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


def _generate_m8_and_hit_summary(
    consolidated_dict,
    added_reads,
    read2blastm8,
    hit_summary_path,
    deduped_blastn_6_path,
    refined_hit_summary_path,
    refined_blastn_6_path,
):
    ''' generate new m8 and hit_summary based on consolidated_dict and read2blastm8 '''
    # Generate new hit summary
    new_read_ids = added_reads.keys()
    with open(hit_summary_path) as hit_summary_f, open(refined_hit_summary_path, "w") as refined_hit_summary_f:
        refined_hit_summary_writer = HitSummaryMergedWriter(refined_hit_summary_f)
        for read in HitSummaryReader(hit_summary_f):
            refined_hit_summary_writer.writerow(consolidated_dict[read["read_id"]])
        # add the reads that are newly blasted
        for read_id in new_read_ids:
            refined_hit_summary_writer.writerow(added_reads[read_id])
    # Generate new M8
    with open(deduped_blastn_6_path) as deduped_blastn_6_f, open(refined_blastn_6_path, "w") as refined_blastn_6_f:
        refined_blastn_6_writer = BlastnOutput6NTRerankedWriter(refined_blastn_6_f)
        for row in BlastnOutput6Reader(deduped_blastn_6_f):
            new_row = read2blastm8.get(row["qseqid"], row)
            new_row["qseqid"] = row["qseqid"]
            refined_blastn_6_writer.writerow(new_row)

        # add the reads that are newly blasted
        for read_id in new_read_ids:
            new_row = read2blastm8.get(read_id)
            new_row["qseqid"] = read_id
            refined_blastn_6_writer.writerow(new_row)


def generate_contig_taxon_summary(
    read_to_contig_tsv_path,
    m8_file,
    hit_summary,
    db_type,
    deuterostome_db_path,
    taxon_whitelist_path,
    taxon_blacklist_path,
    lineage_db_path,
    top_m8_output,
    refined_hit_summary_output,
    contig_summary_json_output,
    refined_counts_with_dcr_output,
    duplicate_cluster_sizes_path=None,
):
    if db_type.lower() == "nt":
        _get_top_m8_nt(m8_file, top_m8_output)
    elif db_type.lower() == "nr":
        _get_top_m8_nr(m8_file, top_m8_output)

    read_dict, accession_dict, _ = m8.summarize_hits(hit_summary)
    read2contig = {}
    read2contig_bases = {}
    with open(read_to_contig_tsv_path) as f:
        for read_id, contig_id, sequence in csv.reader(f, delimiter="\t"):
            read2contig[read_id] = contig_id
            read2contig_bases[(read_id, contig_id)] = len(sequence)  # these should be roughly equal since we're only using primary hits

    updated_read_dict, read2blastm8, contig2lineage, added_reads = _update_read_dict(
        read2contig,
        top_m8_output,
        read_dict,
        accession_dict,
        db_type,
    )

    should_keep = m8.build_should_keep_filter(
        deuterostome_db_path,
        taxon_whitelist_path,
        taxon_blacklist_path,
    )
    contig_taxon_summary = _generate_taxon_summary(
        read2contig,
        contig2lineage,
        updated_read_dict,
        added_reads,
        db_type,
        duplicate_cluster_sizes_path,
        should_keep,
        read2contig_bases,
    )

    with open(contig_summary_json_output, 'w') as f:
        json.dump(contig_taxon_summary, f)

    refined_m8 = "refined.m8"
    _generate_m8_and_hit_summary(
        updated_read_dict,
        added_reads,
        read2blastm8,
        hit_summary,
        m8_file,
        refined_hit_summary_output,
        refined_m8,
    )

    m8.generate_taxon_count_json_from_m8(
        refined_m8,
        refined_hit_summary_output,
        db_type.upper(),
        lineage_db_path,
        deuterostome_db_path,
        taxon_whitelist_path,
        taxon_blacklist_path,
        duplicate_cluster_sizes_path,
        refined_counts_with_dcr_output,
    )


class PipelineStepBlastContigs(PipelineStep):  # pylint: disable=abstract-method
    """ The BLAST step is run independently for the contigs. First against the NT-BLAST database
    constructed from putative taxa identified from short read alignments to NCBI NT using GSNAP.
    Then, against the NR-BLAST database constructed from putative taxa identified from short read
    alignments to NCBI NR with Rapsearch2.

    For NT:
    ```
    blast_command
    -query {assembled_contig}
    -db {blast_index_path}
    -out {blast_m8}
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
    -evalue 1e-10
    -max_target_seqs 5000
    -num_threads 16
    ```

    For NR:

    ```
    blast_command
    -query {assembled_contig}
    -db {blast_index_path}
    -out {blast_m8}
    -outfmt 6
    -num_alignments 5
    -num_threads 16
    ```
    """

    # Opening lineage sqlite in parallel for NT and NR sometimes hangs.
    # In fact it's generally not helpful to run NT and NR in parallel in this step,
    # so we may broaden the scope of this mutex.
    cya_lock = multiprocessing.RLock()

    def run(self):
        '''
            1. summarize hits
            2. built blast index
            3. blast assembled contigs to the index
            4. update the summary
        '''
        _align_m8, deduped_m8, hit_summary, orig_counts_with_dcr = self.input_files_local[0]
        assembled_contig, _assembled_scaffold, bowtie_sam, _contig_stats = self.input_files_local[1]
        reference_fasta, = self.input_files_local[2]
        duplicate_cluster_sizes_path, = self.input_files_local[3]

        blast_m8, refined_m8, refined_hit_summary, refined_counts_with_dcr, contig_summary_json, blast_top_m8 = self.output_files_local()

        assert refined_counts_with_dcr.endswith("with_dcr.json"), self.output_files_local()
        assert orig_counts_with_dcr.endswith("with_dcr.json"), self.output_files_local()

        db_type = self.additional_attributes["db_type"]
        no_assembled_results = (
            os.path.getsize(assembled_contig) < MIN_ASSEMBLED_CONTIG_SIZE or
            os.path.getsize(reference_fasta) < MIN_REF_FASTA_SIZE)

        if no_assembled_results:
            # No assembled results or refseq fasta available.
            # Create empty output files.
            command.write_text_to_file(' ', blast_m8)
            command.write_text_to_file(' ', blast_top_m8)
            command.copy_file(deduped_m8, refined_m8)
            command.copy_file(hit_summary, refined_hit_summary)
            command.copy_file(orig_counts_with_dcr, refined_counts_with_dcr)
            command.write_text_to_file('[]', contig_summary_json)
            return  # return in the middle of the function

        (read_dict, accession_dict, _selected_genera) = m8.summarize_hits(hit_summary)
        PipelineStepBlastContigs.run_blast(db_type, blast_m8, assembled_contig, reference_fasta, blast_top_m8)
        read2contig = {}
        PipelineStepRunAssembly.generate_info_from_sam(bowtie_sam, read2contig, duplicate_cluster_sizes_path)

        (updated_read_dict, read2blastm8, contig2lineage, added_reads) = _update_read_dict(
            read2contig, blast_top_m8, read_dict, accession_dict, db_type)
        self.generate_m8_and_hit_summary(updated_read_dict, added_reads, read2blastm8,
                                         hit_summary, deduped_m8,
                                         refined_hit_summary, refined_m8)

        # Generating taxon counts based on updated results
        lineage_db = s3.fetch_reference(
            self.additional_files["lineage_db"],
            self.ref_dir_local,
            allow_s3mi=False)  # Too small to waste s3mi

        deuterostome_db = None
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = s3.fetch_reference(self.additional_files["deuterostome_db"],
                                                 self.ref_dir_local, allow_s3mi=False)  # Too small for s3mi

        blacklist_s3_file = self.additional_files.get('taxon_blacklist', DEFAULT_BLACKLIST_S3)
        taxon_blacklist = s3.fetch_reference(blacklist_s3_file, self.ref_dir_local)

        taxon_whitelist = None
        if self.additional_attributes.get("use_taxon_whitelist"):
            taxon_whitelist = s3.fetch_reference(self.additional_files.get("taxon_whitelist", DEFAULT_WHITELIST_S3),
                                                 self.ref_dir_local)

        with TraceLock("PipelineStepBlastContigs-CYA", PipelineStepBlastContigs.cya_lock, debug=False):
            with log.log_context("PipelineStepBlastContigs", {"substep": "generate_taxon_count_json_from_m8", "db_type": db_type, "refined_counts": refined_counts_with_dcr}):
                m8.generate_taxon_count_json_from_m8(refined_m8, refined_hit_summary, db_type.upper(),
                                                     lineage_db, deuterostome_db, taxon_whitelist, taxon_blacklist,
                                                     duplicate_cluster_sizes_path, refined_counts_with_dcr)

        # generate contig stats at genus/species level
        with log.log_context("PipelineStepBlastContigs", {"substep": "generate_taxon_summary"}):
            contig_taxon_summary = _generate_taxon_summary(
                read2contig,
                contig2lineage,
                updated_read_dict,
                added_reads,
                db_type,
                duplicate_cluster_sizes_path,
                # same filter as applied in generate_taxon_count_json_from_m8
                m8.build_should_keep_filter(deuterostome_db, taxon_whitelist, taxon_blacklist)
            )

        with log.log_context("PipelineStepBlastContigs", {"substep": "generate_taxon_summary_json", "contig_summary_json": contig_summary_json}):
            with open(contig_summary_json, 'w') as contig_outf:
                json.dump(contig_taxon_summary, contig_outf)

        # Upload additional file
        contig2lineage_json = os.path.join(os.path.dirname(contig_summary_json), f"contig2lineage.{db_type}.json")
        with log.log_context("PipelineStepBlastContigs", {"substep": "contig2lineage_json", "contig2lineage_json": contig2lineage_json}):
            with open(contig2lineage_json, 'w') as c2lf:
                json.dump(contig2lineage, c2lf)

        self.additional_output_files_hidden.append(contig2lineage_json)

    @staticmethod
    def run_blast(db_type, blast_m8, *args):
        blast_index_path = os.path.join(os.path.dirname(blast_m8), f"{db_type}_blastindex")
        if db_type == 'nt':
            PipelineStepBlastContigs.run_blast_nt(blast_index_path, blast_m8, *args)
        else:
            assert db_type == 'nr'
            PipelineStepBlastContigs.run_blast_nr(blast_index_path, blast_m8, *args)

    @staticmethod
    def run_blast_nt(blast_index_path, blast_m8, assembled_contig, reference_fasta, blast_top_m8):
        blast_type = 'nucl'
        blast_command = 'blastn'
        command.execute(
            command_patterns.SingleCommand(
                cmd="makeblastdb",
                args=[
                    "-in",
                    reference_fasta,
                    "-dbtype",
                    blast_type,
                    "-out",
                    blast_index_path
                ],
            )
        )
        command.execute(
            command_patterns.SingleCommand(
                cmd=blast_command,
                args=[
                    "-query",
                    assembled_contig,
                    "-db",
                    blast_index_path,
                    "-out",
                    blast_m8,
                    "-outfmt",
                    '6 ' + ' '.join(field for field, _ in BlastnOutput6NTReader.SCHEMA),
                    '-evalue',
                    1e-10,
                    '-max_target_seqs',
                    5000,
                    "-num_threads",
                    16
                ],  # type: ignore
                # We can only pass BATCH_SIZE as an env var.  The default is 100,000 for blastn;  10,000 for blastp.
                # Blast concatenates input queries until they exceed this size, then runs them together, for efficiency.
                # Unfortunately if too many short and low complexity queries are in the input, this can expand too
                # much the memory required.  We have found empirically 10,000 to be a better default.  It is also the
                # value used as default for remote blast.
                env=dict(os.environ, BATCH_SIZE="10000")
            )
        )
        # further processing of getting the top m8 entry for each contig.
        _get_top_m8_nt(blast_m8, blast_top_m8)

    @staticmethod
    def run_blast_nr(blast_index_path, blast_m8, assembled_contig, reference_fasta, blast_top_m8):
        blast_type = 'prot'
        blast_command = 'blastx'
        command.execute(
            command_patterns.SingleCommand(
                cmd="makeblastdb",
                args=[
                    "-in",
                    reference_fasta,
                    "-dbtype",
                    blast_type,
                    "-out",
                    blast_index_path
                ]
            )
        )
        command.execute(
            command_patterns.SingleCommand(
                cmd=blast_command,
                args=[
                    "-query",
                    assembled_contig,
                    "-db",
                    blast_index_path,
                    "-out",
                    blast_m8,
                    "-outfmt",
                    6,
                    "-num_alignments",
                    5,
                    "-num_threads",
                    16
                ]
            )
        )
        # further processing of getting the top m8 entry for each contig.
        _get_top_m8_nr(blast_m8, blast_top_m8)
