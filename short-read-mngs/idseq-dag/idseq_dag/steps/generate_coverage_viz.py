import json
import math
import os

from collections import defaultdict

import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.dict import open_file_db_by_extension
from idseq_dag.util.m8 import MIN_CONTIG_SIZE
from idseq_dag.util.parsing import BlastnOutput6NTRerankedReader

# These constants can be overridden with the additional_attributes dict:
# The maximum number of bins to divide the accession length into when computing coverage.
MAX_NUM_BINS_COVERAGE = 500
# The number of accessions to show per taxon.
# For each taxon, we show all accessions with a contig (even if there are more than num_accessions_per_taxon of them).
# We then add accessions with only reads until we reach num_accessions_per_taxon.
NUM_ACCESSIONS_PER_TAXON = 10

class PipelineStepGenerateCoverageViz(PipelineStep):  # pylint: disable=abstract-method
    """Pipeline step to generate JSON files for coverage viz to
    be consumed by the web app.
    """

    def run(self):
        """
        Extract data from input files.
        Generate coverage viz data.
        Output JSON output files.
        """
        max_num_bins_coverage = self.additional_attributes.get("max_num_bins_coverage", MAX_NUM_BINS_COVERAGE)
        num_accessions_per_taxon = self.additional_attributes.get("num_accessions_per_taxon", NUM_ACCESSIONS_PER_TAXON)
        min_contig_size = self.additional_attributes.get("min_contig_size", MIN_CONTIG_SIZE)

        info_db = s3.fetch_reference(
            self.additional_files["info_db"],
            self.ref_dir_local,
            allow_s3mi=True)
        with open_file_db_by_extension(info_db) as info_dict:
            # Extract data from input files.
            (taxon_data, accession_data, contig_data, read_data) = self.prepare_data(
                self.input_files_local, info_dict, min_contig_size, num_accessions_per_taxon
            )

        # Generate the coverage viz data for each accession.
        coverage_viz_data = self.generate_coverage_viz_data(accession_data, contig_data, read_data, max_num_bins_coverage)

        # Generate the summary data, which contains a dict of all taxons for which coverage viz data is available.
        # For each taxon, summary data for the best accessions, plus the number of total accessions, is included.
        coverage_viz_summary_data = self.generate_coverage_viz_summary_data(taxon_data, accession_data, coverage_viz_data)

        coverage_viz_summary = self.output_files_local()[0]
        # Write the summary JSON file which is initially loaded on the report page.
        with open(coverage_viz_summary, 'w') as cvs:
            json.dump(coverage_viz_summary_data, cvs)

        # Create a separate coverage viz JSON file for each accession.
        # This file will be passed to the front-end when the user views that particular accession.
        coverage_viz_dir = os.path.join(self.output_dir_local, "coverage_viz")
        command.make_dirs(coverage_viz_dir)
        for accession_id in coverage_viz_data:
            upload_file = os.path.join(coverage_viz_dir, f"{accession_id}_coverage_viz.json")

            with open(upload_file, 'w') as uf:
                json.dump(coverage_viz_data[accession_id], uf)

        self.additional_output_folders_hidden.append(coverage_viz_dir)

    @staticmethod
    def prepare_data(input_files_local, info_dict, min_contig_size, num_accessions_per_taxon):
        """
        Extract taxon, accession, contig, and read data from input files.
        Remove taxons with no contigs in any accession.
        Select the best accessions for each taxon.
        """
        (_reassigned_m8, hit_summary, blast_top_m8) = input_files_local[0]
        (contig_coverage_json, contig_stats_json, contigs_fasta) = input_files_local[1]
        (gsnap_deduped_m8, ) = input_files_local[2]

        # Get a map from valid contigs (based on min_contig_size) to number of reads in that contig.
        valid_contigs_with_read_counts = PipelineStepGenerateCoverageViz.get_valid_contigs_with_read_counts(contig_stats_json, min_contig_size)

        # Get a map from accession name to hits (reads and contigs)
        # Also get a map from taxons to accession names.
        (accession_data, taxon_data) = PipelineStepGenerateCoverageViz.generate_accession_data(hit_summary, valid_contigs_with_read_counts)
        # Remove taxons with no contigs in any of their accessions.
        PipelineStepGenerateCoverageViz.remove_taxons_with_no_contigs(accession_data, taxon_data)

        # Add the total_length and name of the accession to the accession data.
        PipelineStepGenerateCoverageViz.augment_accession_data_with_info(info_dict, accession_data)

        # Get unassigned reads. Use a set for performance.
        unassigned_reads_set = PipelineStepGenerateCoverageViz.get_unassigned_reads_set(accession_data)

        # Extract information about contigs and reads.
        contig_data = PipelineStepGenerateCoverageViz.generate_contig_data(blast_top_m8, valid_contigs_with_read_counts)
        read_data = PipelineStepGenerateCoverageViz.generate_read_data(gsnap_deduped_m8, unassigned_reads_set)

        # Add coverage to the contig data.
        PipelineStepGenerateCoverageViz.augment_contig_data_with_coverage(contig_coverage_json, contig_data)

        # Add byteranges to the contig data.
        PipelineStepGenerateCoverageViz.augment_contig_data_with_byteranges(contigs_fasta, contig_data)

        # Select the best accessions for each taxon.
        (taxon_data, accession_data) = PipelineStepGenerateCoverageViz.select_best_accessions_per_taxon(
            taxon_data, accession_data, contig_data, read_data, num_accessions_per_taxon
        )

        return (
            taxon_data,
            accession_data,
            contig_data,
            read_data
        )

    @staticmethod
    def generate_coverage_viz_data(accession_data, contig_data, read_data, max_num_bins):
        """
        Generate coverage viz data for each accession in accession_data.
        """
        coverage_viz_data = {}

        for accession_id, accession_obj in accession_data.items():
            total_length = accession_obj["total_length"]

            # Number of bins to calculate coverage for.
            num_bins = min(max_num_bins, total_length)

            # Aggregate the reads and contigs into "hit groups".
            # Divide the accession up into a number of bins, and group together small reads and contigs that fall in the same bin.
            hit_groups = PipelineStepGenerateCoverageViz.generate_hit_group_json(accession_obj, accession_id, contig_data, read_data, num_bins)

            # Calculate the coverage for the accession, based on the reads and contigs.
            (coverage, coverage_bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
                accession_id, accession_obj, contig_data, read_data, num_bins
            )

            # Calculate statistics for the accession.
            accession_stats = PipelineStepGenerateCoverageViz.calculate_accession_stats(
                accession_obj,
                contig_data,
                read_data
            )

            coverage_viz_data[accession_id] = {
                "total_length": total_length,
                "name": accession_obj["name"],
                "hit_groups": hit_groups,
                "coverage": coverage,
                "coverage_bin_size": coverage_bin_size,
                "max_aligned_length": accession_stats["max_aligned_length"],
                "coverage_depth": _format_number(accession_stats["coverage_depth"]),
                "coverage_breadth": _format_percent(accession_stats["coverage_breadth"]),
                "avg_prop_mismatch": _format_percent(accession_stats["avg_prop_mismatch"]),
            }

        return coverage_viz_data

    @staticmethod
    def generate_coverage_viz_summary_data(taxon_data, accession_data, coverage_viz_obj):
        """
        Generate a summary dict that maps taxons to accession data.
        """
        coverage_viz_summary_data = {}

        for taxon, data in taxon_data.items():
            coverage_viz_summary_data[taxon] = {
                "best_accessions": list(map(lambda accession_id: {
                    "id": accession_id,
                    "name": accession_data[accession_id]["name"],
                    "num_contigs": len(accession_data[accession_id]["contigs"]),
                    "num_reads": len(accession_data[accession_id]["reads"]),
                    "score": _format_number(accession_data[accession_id]["score"]),
                    "coverage_depth": coverage_viz_obj[accession_id]["coverage_depth"]
                }, data["accessions"])),
                "num_accessions": data["num_total_accessions"]
            }

        return coverage_viz_summary_data

    @staticmethod
    def get_valid_contigs_with_read_counts(contig_stats_json, min_contig_size):
        """
        Get a dict that maps valid contigs to their read count.
        A contig is valid if it is larger than min_contig_size.
        """
        valid_contigs_with_read_counts = {}
        with open(contig_stats_json, 'r') as csj:
            contig_read_counts = json.load(csj)

            for contig in contig_read_counts:
                if contig != "*" and contig_read_counts[contig] >= min_contig_size:
                    valid_contigs_with_read_counts[contig] = contig_read_counts[contig]

        return valid_contigs_with_read_counts

    @staticmethod
    def generate_accession_data(hit_summary, valid_contigs_with_read_counts):
        """
        Generate a dict that maps accessions to the reads and contigs that were assigned to them.
        Also generate a dict that maps taxons to accessions.
        """
        # Use a set for contigs, since multiple lines in the hitsummary can map to the same contig.
        accession_data = defaultdict(lambda: {'reads': [], 'contigs': set()})
        # Use a set for accessions, since multiple lines in the hitsummary can map to the same accession.
        taxon_data = defaultdict(lambda: {'accessions': set(), 'num_total_accessions': 0})

        line_count = 0
        with open(hit_summary, 'r') as hs:
            for line in hs:
                line_count += 1
                if line_count % 100000 == 0:
                    log.write(f"{line_count} lines in the hitsummary file processed.")

                values = line.rstrip().split("\t")

                # Only count the contig if the contig is valid, i.e. it is larger than min_contig_size.
                # For a particular line in hit_summary, if the line only has 7 columns, this means the read wasn't assembled.
                # If the line has 12 or 13 columns, then the read was assembled into a contig.
                if len(values) >= 12 and values[7] in valid_contigs_with_read_counts:
                    taxon_data[values[9]]["accessions"].add(values[8])
                    accession_data[values[8]]["contigs"].add(values[7])
                else:
                    taxon_data[values[4]]["accessions"].add(values[3])
                    accession_data[values[3]]["reads"].append(values[0])

        # Convert the contig set to a list.
        for accession_id, data in accession_data.items():
            accession_data[accession_id]["contigs"] = list(data["contigs"])

        # Add the number of total accessions to taxon_data.
        # Later, accessions will be filtered to just the best accessions.
        for taxon, data in taxon_data.items():
            taxon_data[taxon]["num_total_accessions"] = len(data["accessions"])

        return (accession_data, taxon_data)

    @staticmethod
    def remove_taxons_with_no_contigs(accession_data, taxon_data):
        """
        Removes taxons with no assembled contigs from taxon_data.
        Also removes all accessions corresponding to these taxons from accession_data
        Mutates input parameters.
        """
        taxons_to_remove = []

        # Find all taxons with no mapped contigs.
        for taxon, data in taxon_data.items():
            total_contigs = sum([len(accession_data[accession_id]["contigs"]) for accession_id in data["accessions"]])

            if total_contigs == 0:
                taxons_to_remove.append(taxon)

        # Delete all taxons and corresponding accessions.
        for taxon in taxons_to_remove:
            accessions = list(taxon_data[taxon]["accessions"])
            del taxon_data[taxon]
            for accession_id in accessions:
                del accession_data[accession_id]

    @staticmethod
    def augment_accession_data_with_info(info_dict, accession_data):
        """
        Augment the accession data dictionary with accession data pulled from the info_db.
        """
        for accession_id in accession_data:
            entry = info_dict.get(accession_id)

            if entry:
                # The name of the accession.
                accession_data[accession_id]["name"] = entry[0]
                # The total length of the accession in base pairs.
                accession_data[accession_id]["total_length"] = int(entry[1])
            else:
                accession_data[accession_id]["name"] = "Unknown accession"
                accession_data[accession_id]["total_length"] = 0

    @staticmethod
    def get_unassigned_reads_set(accession_data):
        """
        Get the set of unassigned reads from the accession_data.
        """
        unassigned_reads_set = set()
        for _, accession_obj in accession_data.items():
            for read in accession_obj["reads"]:
                unassigned_reads_set.add(read)

        return unassigned_reads_set

    @staticmethod
    def generate_hit_data_from_m8(blastn_6_path, valid_hits):
        """
        Generate hit data from an m8 file.
        Only include hits whose name appears in the valid_hits collection.
        """
        # M8 file should have at least a single line.
        # Anything less than this is considered an empty file.
        MIN_M8_FILE_SIZE = 25
        hits = {}

        # File is empty.
        if os.path.getsize(blastn_6_path) < MIN_M8_FILE_SIZE:
            return hits

        with open(blastn_6_path) as blastn_6_f:
            for hit in BlastnOutput6NTRerankedReader(blastn_6_f):

                if hit["qseqid"] in valid_hits:
                    # Blast output is per HSP, yet the hit represents a set of HSPs,
                    # so these fields have been aggregated across that set by
                    # function summary_row() in class CandidateHit.
                    hits[hit["qseqid"]] = {
                        "accession": hit["sseqid"],
                        "percent_id": hit["pident"],
                        "alignment_length": hit["length"],
                        "num_mismatches": hit["mismatch"],
                        "num_gaps": hit["gapopen"],
                        "query_start": hit["qstart"],
                        "query_end": hit["qend"],
                        "subject_start": hit["sstart"],
                        "subject_end": hit["send"],
                        "prop_mismatch": hit["mismatch"] / max(1, hit["length"]),
                    }

            return hits

    @staticmethod
    def generate_contig_data(blast_top_m8, valid_contigs_with_read_counts):
        """
        Generate contig data from blast_top_m8.
        """
        contigs = PipelineStepGenerateCoverageViz.generate_hit_data_from_m8(blast_top_m8, valid_contigs_with_read_counts)

        # Include some additional data.
        for contig_id, contig_obj in contigs.items():
            name_parts = contig_id.split("_")
            # Total length of the contig. We extract this from the contig name.
            contig_obj["total_length"] = int(name_parts[3])
            # The contig read count.
            contig_obj["num_reads"] = valid_contigs_with_read_counts[contig_id]

        return contigs

    @staticmethod
    def generate_read_data(gsnap_deduped_m8, unassigned_reads_set):
        """
        Generate read data from gnap_deduped_m8.
        We process gsnap.deduped.m8 instead of gsnap.reassigned.m8 because we ignore contigs with read_count < 4.
        However, these contigs still get reassigned in gsnap.reassigned.m8,
        and overwrite the original read alignment to the accession, which we need. So we can't use gsnap.reassigned.m8.
        """
        return PipelineStepGenerateCoverageViz.generate_hit_data_from_m8(gsnap_deduped_m8, unassigned_reads_set)

    @staticmethod
    def augment_contig_data_with_coverage(contig_coverage_json, contig_data):
        """
        Augment contig data with the contig coverage array.
        """
        with open(contig_coverage_json, 'r') as ccj:
            contig_coverage = json.load(ccj)

            for contig_name in contig_coverage:
                if contig_name in contig_data:
                    contig_data[contig_name]["coverage"] = contig_coverage[contig_name]["coverage"]

    @staticmethod
    def augment_contig_data_with_byteranges(contigs_fasta, contig_data):
        """
        Augment contig data with the byterange location in the contigs.fasta file for each contig.
        """
        with open(contigs_fasta, 'r') as cf:
            seq_offset = 0
            seq_len = 0
            contig_name = ""

            # Process each line in contigs.fasta.
            for line in cf:
                # If the line is a header file, process the contig we just traversed.
                if line[0] == '>':  # header line
                    if seq_len > 0 and contig_name in contig_data:
                        contig_data[contig_name]["byterange"] = [seq_offset, seq_len]

                    seq_offset = seq_offset + seq_len
                    seq_len = len(line)
                    contig_name = line[1:].strip()
                else:
                    seq_len += len(line)

            # Process the last contig once we reach the end of the file.
            if seq_len > 0 and contig_name in contig_data:
                contig_data[contig_name]["byterange"] = [seq_offset, seq_len]

    @staticmethod
    def select_best_accessions_per_taxon(taxon_data, accession_data, contig_data, _read_data, num_accessions_per_taxon):
        """
        Select the accessions that we will generate coverage viz data for.
        We select all accessions that have a contig.
        If there are too few such accessions, we also select the best accessions that only have reads until we have enough.
        """
        def get_score(accession_id):
            """
            Calculate a score for each accession.
            """
            accession_obj = accession_data[accession_id]

            contig_lengths = [contig_data[contig_name]["alignment_length"] for contig_name in accession_obj["contigs"]]

            max_contig_length = max(contig_lengths) if len(contig_lengths) > 0 else 0
            total_contig_length = sum(contig_lengths) if len(contig_lengths) > 0 else 0
            num_reads = len(accession_obj["reads"])

            # We add max_contig_length and total_contig_length in order to give some additional weight to long contigs.
            # We use the length of contigs, but only the number of reads, so contigs should dominate this score whenever they are present.
            return max_contig_length + total_contig_length + num_reads

        filtered_taxon_data = {}
        filtered_accession_data = {}

        for taxon, data in taxon_data.items():
            for accession_id in data["accessions"]:
                accession_data[accession_id]["score"] = get_score(accession_id)

            # Sort the accessions
            sorted_accessions = sorted(data["accessions"], key=lambda accession_id: accession_data[accession_id]["score"], reverse=True)

            # Take ALL accessions with one or more contigs.
            accessions_with_contigs = list(filter(lambda accession_id: len(accession_data[accession_id]["contigs"]) >= 1, sorted_accessions))

            # If this is enough accessions, we are done.
            if len(accessions_with_contigs) >= num_accessions_per_taxon:
                best_accessions = accessions_with_contigs
            # Otherwise, fill the remaining spots with accessions without contigs.
            else:
                # Get all accessions without contigs.
                accessions_without_contigs = list(filter(lambda accession_id: len(accession_data[accession_id]["contigs"]) == 0, sorted_accessions))
                # Fill in the remaining spots
                filtered_accessions = accessions_with_contigs + accessions_without_contigs[0: num_accessions_per_taxon - len(accessions_with_contigs)]
                # Re-sort the accessions, as it's possible that an accession with no contigs has a higher score than an accession with contigs.
                best_accessions = sorted(filtered_accessions, key=lambda accession_id: accession_data[accession_id]["score"], reverse=True)

            filtered_taxon_data[taxon] = {
                "accessions": best_accessions,
                "num_total_accessions": taxon_data[taxon]["num_total_accessions"]
            }

            # Populate filtered_accession_data with only the best accessions.
            for accession_id in filtered_taxon_data[taxon]["accessions"]:
                filtered_accession_data[accession_id] = accession_data[accession_id]

        return (filtered_taxon_data, filtered_accession_data)

    @staticmethod
    def generate_hit_group_json(accession_data, accession_id, contig_data, read_data, num_bins):
        """
        In order to display numerous tiny hits in the coverage viz, we aggregate them into groups.
        """
        # If a read of contig is larger than the bin size, treat it as an individual hit instead of aggregating it.
        individual_reads = []
        individual_contigs = []

        # Initialize bins.
        bin_size = accession_data["total_length"] / num_bins
        read_bins = [[] for i in range(num_bins)]
        contig_bins = [[] for i in range(num_bins)]

        # Add each hit to the appropriate array: individual_reads, individual_contigs, read_bins, or contig_bins.
        def process_hit(hit_type, hit_name):
            hit_data = read_data if hit_type == "read" else contig_data

            if hit_name not in hit_data:
                log.write(f"Could not find {hit_type} in map: {hit_name}")
                return

            hit_obj = hit_data[hit_name]

            # hitsummary is more strict than reassigned.
            # Sometimes reassigned will have a value for accession, but hitsummary won't.
            if hit_obj["accession"] != accession_id:
                log.write(f"Mismatched accession for {hit_name}: {hit_obj['accession']} (reassigned) versus {accession_id} (hitsummary)")
                return

            (accession_start, accession_end) = _align_interval(_decrement_lower_bound((hit_obj["subject_start"], hit_obj["subject_end"])))

            # If the hit is larger than the bin size, treat it as an individual hit.
            if accession_end - accession_start >= bin_size:
                if hit_type == "read":
                    individual_reads.append(hit_name)
                else:
                    individual_contigs.append(hit_name)

            # Otherwise, put the hit into a bin based on its midpoint
            else:
                hit_midpoint = (accession_end + accession_start) / 2
                hit_bin_index = _floor_with_min(hit_midpoint / bin_size, 0)

                if hit_type == "read":
                    read_bins[hit_bin_index].append(hit_name)
                else:
                    contig_bins[hit_bin_index].append(hit_name)

        for read_name in accession_data["reads"]:
            process_hit("read", read_name)

        for contig_name in accession_data["contigs"]:
            process_hit("contig", contig_name)

        # Generate the hit group JSON for individual hits.
        hit_groups = []
        for read_name in individual_reads:
            read_obj = read_data[read_name]
            hit_groups.append(PipelineStepGenerateCoverageViz.get_hit_group_json([], [read_obj], bin_size))

        for contig_name in individual_contigs:
            contig_obj = contig_data[contig_name]
            hit_groups.append(PipelineStepGenerateCoverageViz.get_hit_group_json([contig_obj], [], bin_size))

        # Generate the hit group JSON for aggregated hits.
        for i in range(num_bins):
            reads = read_bins[i]
            contigs = contig_bins[i]

            # Ignore empty bins.
            if len(reads) + len(contigs) == 0:
                continue
            else:
                read_objs = list(map(lambda read_name: read_data[read_name], reads))
                contig_objs = list(map(lambda contig_name: contig_data[contig_name], contigs))
                hit_groups.append(PipelineStepGenerateCoverageViz.get_hit_group_json(contig_objs, read_objs, bin_size))

        return hit_groups

    @staticmethod
    def get_hit_group_json(contig_objs, read_objs, bin_size):
        """
        Generate the JSON for a group of hits (reads and contigs) that are being aggregated into a single bin.
        """
        num_contigs = len(contig_objs)
        num_reads = len(read_objs)

        # Calculate some stats that only apply to contig_objs.
        contig_r = sum([contig_obj["num_reads"] for contig_obj in contig_objs])
        contig_byteranges = [contig_obj["byterange"] for contig_obj in contig_objs]

        # Treat read_objs and contig_objs the same from here onwards. They share many of the same fields.
        hit_objs = contig_objs + read_objs
        num_hits = num_contigs + num_reads

        # Averaging a particular field across all hits.
        def avg_field(field):
            return sum(map(lambda hit_obj: hit_obj[field], hit_objs)) / num_hits

        # Calculate the bounds of this hit group.
        endpoints = []

        for hit_obj in hit_objs:
            endpoints.append(hit_obj["subject_start"])
            endpoints.append(hit_obj["subject_end"])

        hit_group_start = min(endpoints)
        hit_group_end = max(endpoints)

        # Also calculate the bin index based on the midpoint of the hit group.
        # If the bounds of the hit group is too narrow,
        # the front-end will display the bounds of the bin instead.
        hit_group_midpoint = ((hit_group_start - 1) + hit_group_end) / 2
        hit_group_bin_index = _floor_with_min(hit_group_midpoint / bin_size, 0)

        # Use an array of numbers instead of a dict with field names to save space in the JSON file.
        # There will be many of these hit group arrays.
        return [
            num_contigs,
            num_reads,
            # Total number of contig reads.
            contig_r,
            # Alignment range
            hit_group_start,
            hit_group_end,
            # Alignment length. Can be different from alignment range.
            _format_number(avg_field("alignment_length")),
            # Percent identity
            _format_percent(avg_field("percent_id") / 100),
            # Number of mismatches
            _format_number(avg_field("num_mismatches")),
            # Number of gaps
            _format_number(avg_field("num_gaps")),
            # Bin index of midpoint.
            hit_group_bin_index,
            # Byteranges in the contigs.fasta file for each contig.
            contig_byteranges
        ]

    @staticmethod
    def calculate_accession_coverage(accession_id, accession_data, contig_data, read_data, num_bins):
        """
        Divide the accession length into a number of bins, and calculate the average coverage for each bin.
        """
        bin_size = accession_data["total_length"] / num_bins
        coverage = [{
            "depth": 0,
            "endpoints": [],
            "num_reads": 0,
            "num_contigs": 0
        } for i in range(num_bins)]

        # First, process each contig.
        # For each contig, figure out which bins the contig overlaps.
        # For each bin, figure out how much coverage the contig contribtes to that bin.
        for contig_name in accession_data["contigs"]:
            if contig_name not in contig_data:
                log.write(f"Could not find contig in contig data: {contig_name}")
                continue

            contig_obj = contig_data[contig_name]
            # Ignore contigs with accession mismatch
            if contig_obj["accession"] != accession_id:
                continue

            # The bins and coverage array are 0-indexed, but subject start/end and coverage start/end are 1-indexed.
            # We convert everything to 0-index here and stay in 0-index for the rest of the function.
            # NOTE: We decrement only the lower bound here so that we can treat the discrete integer indices as a continuous interval
            # while converting from accession interval to contig interval to contig coverage interval. This makes the math easier.
            # These conversions are necessary because the accession interval, contig interval, and contig coverage interval
            # might all be different sizes.
            # We convert back to integer indices when we calculate coverage_arr_start/_end.
            (subject_start, subject_end) = _decrement_lower_bound((contig_obj["subject_start"], contig_obj["subject_end"]))
            (query_start, query_end) = _decrement_lower_bound((contig_obj["query_start"], contig_obj["query_end"]))

            # Find all bins that this contig overlaps, and calculate average coverage for each bin separately.
            (bin_start, bin_end) = _align_interval((subject_start / bin_size, subject_end / bin_size))

            for i in range(_floor_with_min(bin_start, 0), _ceil_with_max(bin_end, num_bins)):
                # Our goal is to figure out which part of the contig coverage array corresponds to this bin.
                # Get the section of the accession that corresponds to the current bin and overlaps with the contig.
                accession_interval = [bin_size * max(bin_start, i), bin_size * min(bin_end, i + 1)]

                # Convert the accession interval to a section of the contig by using the alignment data.
                contig_interval = _transform_interval(accession_interval, subject_start, subject_end, query_start, query_end)

                # The contig coverage array should be the same length as the contig length.
                # If not, convert to the appropriate range in the coverage array.
                if contig_obj["total_length"] == len(contig_obj["coverage"]):
                    coverage_interval = _align_interval((contig_interval[0], contig_interval[1]))
                else:
                    coverage_interval = _transform_interval(contig_interval, 0, contig_obj["total_length"], 0, len(contig_obj["coverage"]))
                    coverage_interval = _align_interval((coverage_interval[0], coverage_interval[1]))

                # Convert back to integer indices.
                # This is the range of values in the contig coverage array that corresponds to the section of the contig that overlaps with this bin.
                (coverage_arr_start, coverage_arr_end) = (_floor_with_min(coverage_interval[0], 0), _ceil_with_max(coverage_interval[1], len(contig_obj["coverage"])))

                # Guard against a division-by-zero bug caused a rounding error.
                # There are circumstances where a contig might have (bin_start, bin_end) = (200, 477.06) but with rounding errors this becomes (199.9999997, 477.06).
                # This causes us to attempt to process bin 199 for the interval (199.99999997, 200). This interval is so small that
                # the coverage interval for the contig ends up being (coverage_arr_start, coverage_arr_end) = (322.0, 322.0) and having length 0.
                # In normal cases, this interval should have at least length 1 because of the floor and ceil.
                # We should just disregard this edge case, because the contig doesn't really overlap this bin (it's a rounding error)
                if coverage_arr_end - coverage_arr_start > 0:
                    # Get the average coverage for the section of the contig that overlaps with this bin.
                    avg_coverage_for_coverage_interval = sum(contig_obj["coverage"][coverage_arr_start: coverage_arr_end]) / (coverage_arr_end - coverage_arr_start)

                    # Multiply by the proportion of the bin that the contig covers.
                    avg_coverage_for_bin = avg_coverage_for_coverage_interval * (abs(accession_interval[1] - accession_interval[0]) / bin_size)

                    coverage[i]["depth"] += avg_coverage_for_bin
                    coverage[i]["endpoints"].append([max(i * bin_size, accession_interval[0]), 1])
                    coverage[i]["endpoints"].append([min((i + 1) * bin_size, accession_interval[1]), -1])
                    coverage[i]["num_contigs"] += 1

        # The logic for processing reads is very similar to contigs above, but the avg coverage on the read is simply 1.
        for read_name in accession_data["reads"]:
            if read_name not in read_data:
                log.write(f"Could not find read in read data: {read_name}")
                continue

            read_obj = read_data[read_name]
            # Ignore reads with accession mismatch
            if read_obj["accession"] != accession_id:
                continue

            (subject_start, subject_end) = _decrement_lower_bound((read_obj["subject_start"], read_obj["subject_end"]))

            # Find all bins that this read overlaps, and calculate average coverage for each bin separately.
            (bin_start, bin_end) = _align_interval((subject_start / bin_size, subject_end / bin_size))

            for i in range(_floor_with_min(bin_start, 0), _ceil_with_max(bin_end, num_bins)):
                # Get the section of the accession that corresponds to the current bin and overlaps with the read.
                accession_range = [bin_size * max(bin_start, i), bin_size * min(bin_end, i + 1)]

                # The read coverage is 1. Multiply by the proportion of the bin that the read covers.
                avg_coverage_for_bin = (abs(accession_range[1] - accession_range[0]) / bin_size)

                coverage[i]["depth"] += avg_coverage_for_bin
                coverage[i]["endpoints"].append([max(i * bin_size, accession_range[0]), 1])
                coverage[i]["endpoints"].append([min((i + 1) * bin_size, accession_range[1]), -1])
                coverage[i]["num_reads"] += 1

        final_coverage = []

        # For each index, an array of numbers is generated.
        # The array is sparse. Only bins with nonzero coverage are included.
        for index in range(len(coverage)):
            coverage_obj = coverage[index]
            # Ignore all bins with no coverage.
            if coverage_obj["depth"] > 0:
                # Use an array of numbers instead of a dict with field names to save space in the JSON file.
                # There will be many of these coverage arrays.
                final_coverage.append([
                    index,  # bin index
                    _format_number(coverage_obj["depth"]),  # average coverage depth
                    _format_percent(PipelineStepGenerateCoverageViz.calculate_covered_length(coverage_obj["endpoints"]) / bin_size),  # coverage breadth
                    coverage_obj["num_contigs"],  # number of contigs
                    coverage_obj["num_reads"],  # number of reads
                ])

        return (final_coverage, bin_size)

    @staticmethod
    def calculate_covered_length(endpoints):
        """
        Calculate the total distance covered by a set of overlapping interval endpoints.
        The endpoints were collected earlier while processing reads and contigs.
        This is used for calculating coverage breadth.
        Each interval endpoint is [bp, depth_change], where bp is the base pair where the endpoint is located,
        and depth_change is 1 for an interval start and -1 for an interval end.
        """
        total_covered_length = 0
        cur_depth = 0
        last_start_point = 0

        endpoints.sort(key=lambda endpoint: endpoint[0])

        for endpoint in endpoints:
            cur_depth += endpoint[1]

            if cur_depth < 0:
                raise ValueError("coverage depth of %s is invalid. Malformed endpoints" % cur_depth)

            # Covered length is starting. Set the last start point.
            if endpoint[1] == 1 and cur_depth == 1:
                last_start_point = endpoint[0]
            # Covered length is ending. Add the covered distance to the total.
            elif endpoint[1] == -1 and cur_depth == 0:
                total_covered_length += endpoint[0] - last_start_point

        if cur_depth != 0:
            raise ValueError("coverage depth is %s after traversal. 0 is expected" % cur_depth)

        return total_covered_length

    @staticmethod
    def calculate_accession_stats(accession_data, contig_data, read_data):
        """
        Calculate various statistics for the accession.
        """
        max_aligned_length = 0
        coverage_sum = 0
        endpoints = []
        prop_total_mismatch = 0

        for contig_name in accession_data["contigs"]:
            if contig_name not in contig_data:
                log.write(f"Could not find contig in contig data: {contig_name}")
                continue

            contig_obj = contig_data[contig_name]

            (accession_start, accession_end) = _align_interval(_decrement_lower_bound((contig_obj["subject_start"], contig_obj["subject_end"])))
            (contig_start, contig_end) = _align_interval(_decrement_lower_bound((contig_obj["query_start"], contig_obj["query_end"])))

            # For max_aligned_length
            accession_alignment_length = accession_end - accession_start
            if accession_alignment_length > max_aligned_length:
                max_aligned_length = accession_alignment_length

            # For coverage_depth
            # Restrict to the part of the coverage that corresponds to the alignment.
            coverage_sum += sum(contig_obj["coverage"][contig_start: contig_end])

            # For avg_prop_mismatch
            prop_total_mismatch += contig_obj["prop_mismatch"]

            # For coverage_breadth
            endpoints.append([accession_start, 1])
            endpoints.append([accession_end, -1])

        for read_name in accession_data["reads"]:
            if read_name not in read_data:
                log.write(f"Could not find read in read data: {read_name}")
                continue

            read_obj = read_data[read_name]

            (accession_start, accession_end) = _align_interval(_decrement_lower_bound((read_obj["subject_start"], read_obj["subject_end"])))

            # For max_aligned_length
            read_length = accession_end - accession_start
            if read_length > max_aligned_length:
                max_aligned_length = read_length

            # For coverage_depth
            coverage_sum += read_length

            # For avg_prop_mismatch
            prop_total_mismatch += read_obj["prop_mismatch"]

            # For coverage_breadth
            endpoints.append([accession_start, 1])
            endpoints.append([accession_end, -1])

        return {
            "max_aligned_length": max_aligned_length,
            # Divide the coverage sum by the total length.
            "coverage_depth": coverage_sum / accession_data["total_length"],
            # Calculate the total length covered by hits, and divide by the total length of the accession.
            "coverage_breadth": PipelineStepGenerateCoverageViz.calculate_covered_length(endpoints) / accession_data["total_length"],
            # Sum up the total prop mismatch and divide by number of contigs and reads.
            "avg_prop_mismatch": prop_total_mismatch / (len(accession_data["contigs"]) + len(accession_data["reads"]))
        }


# Private utility methods used by the pipeline step.

def _decrement_lower_bound(interval):
    """
    Decrement the lower bound of a possibly inverted interval.
    i.e. turn (1, 5) into (0, 5) and (5, 1) into (5, 0).
    Intervals in m8 alignment files can be inverted.
    """
    (bound_one, bound_two) = interval
    if bound_one < bound_two:
        return (bound_one - 1, bound_two)
    else:
        return (bound_one, bound_two - 1)

def _align_interval(interval):
    """
    Flip inverted intervals so the lower number is first.
    """
    (bound_one, bound_two) = interval
    return (min(bound_one, bound_two), max(bound_one, bound_two))

def _transform_interval(
        interval, first_domain_start, first_domain_end, second_domain_start, second_domain_end
):
    """
    Transform an interval from one domain to another.
    The interval should be within the first domain [first_domain_start, first_domain_end]
    For example, _transform_interval((3, 5), 0, 10, 100, 200) would return (130, 150)
    """
    def transform_value(value):
        position = (value - first_domain_start) / (first_domain_end - first_domain_start)
        return position * (second_domain_end - second_domain_start) + second_domain_start

    return [transform_value(value) for value in interval]

def _format_number(number):
    """
    Formatter that adds a couple more sig-figs to small numbers so they don't show up as "0.0"
    """
    if number < 0.1:
        return round(number, 3)
    if number < 1:
        return round(number, 2)
    return round(number, 1)

def _format_percent(number):
    return round(number, 3)

def _round_if_within_epsilon(num, epsilon=0.001):
    r_num = round(num)

    if abs(r_num - num) < epsilon:
        return r_num
    else:
        return num

def _floor_with_min(num, min_value):
    """
    Floor function that also takes a min value, to deal with possible rounding errors.
    Rounding errors in floats can cause out of index errors when we index into an array.
    """
    return max(math.floor(num), min_value)

def _ceil_with_max(num, max_value):
    """
    Ceil function that also takes a max value, to deal with possible rounding errors.
    Rounding errors in floats can cause out of index errors when we index into an array.
    """
    return min(math.ceil(num), max_value)
