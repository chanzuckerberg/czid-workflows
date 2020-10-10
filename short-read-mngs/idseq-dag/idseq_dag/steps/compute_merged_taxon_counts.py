import json
import csv

from collections import namedtuple

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3, DEFAULT_WHITELIST_S3
from idseq_dag.util.m8 import BLAST_OUTPUT_SCHEMA, generate_taxon_count_json_from_m8, parse_tsv
from idseq_dag.util.schemas import TAB_SCHEMA_MERGED
from idseq_dag.util.s3 import fetch_reference


SpeciesAlignmentResults = namedtuple('SpeciesAlignmentResults', ['contig', 'read'])

class ComputeMergedTaxonCounts(PipelineStep):
    '''
    Merge taxon results from NT database and from NR database
    '''

    def run(self):
        # Retrieve inputs and outputs
        nt_m8, nt_hitsummary2_tab, nr_m8, nr_hitsummary2_tab, cluster_sizes_filename = self.input_files_local
        merged_m8_filename, merged_hit_filename, merged_taxon_count_filename = self.output_files_local()

        # Create new merged m8 and hit summary files
        nr_alignment_per_read = {}
        # TODO: consider going through nt first, since that will save us from storing results in memory
        # for all the reads that get they hit from NT contigs
        for nr_hit_dict in parse_tsv(nr_hitsummary2_tab, TAB_SCHEMA_MERGED, strict=False):
            nr_alignment_per_read[nr_hit_dict["read_id"]] = SpeciesAlignmentResults(
                contig=nr_hit_dict.get("contig_species_taxid"),
                read=nr_hit_dict.get("species_taxid"),
            )

        with open(merged_m8_filename, 'w') as output_m8, open(merged_hit_filename, 'w') as output_hit:
            # first pass for NR and output to m8 files if assignment should come from NT
            for nt_hit_dict, [nt_m8_dict, nt_m8_row] in zip(
                parse_tsv(nt_hitsummary2_tab, TAB_SCHEMA_MERGED, strict=False),
                parse_tsv(nt_m8, BLAST_OUTPUT_SCHEMA, raw_lines=True, strict=False)
            ):
                # assert files match
                assert nt_hit_dict['read_id'] == nt_m8_dict["qseqid"], f"Mismatched m8 and hit summary files for nt [{nt_hit_dict['read_id']} != {nt_m8_dict['qseqid']}]"

                nr_alignment = nr_alignment_per_read.get(nt_hit_dict["read_id"])
                nt_alignment = SpeciesAlignmentResults(
                    contig=nt_hit_dict.get("contig_species_taxid"),
                    read=nt_hit_dict.get("species_taxid")
                )
                if nt_alignment.contig:
                    output_m8.write(nt_m8_row)
                    nt_hit_dict["source_count_type"] = "NT"
                    self._write_tsv_row(nt_hit_dict, TAB_SCHEMA_MERGED, output_hit)
                    if nr_alignment:
                        del nr_alignment_per_read[nt_hit_dict["read_id"]]
                elif nr_alignment and nr_alignment.contig:
                    continue
                elif nt_alignment.read:
                    output_m8.write(nt_m8_row)
                    nt_hit_dict["source_count_type"] = "NR"
                    self._write_tsv_row(nt_hit_dict, TAB_SCHEMA_MERGED, output_hit)
                    if nr_alignment:
                        del nr_alignment_per_read[nt_hit_dict["read_id"]]
                elif nr_alignment and nr_alignment.read:
                    continue
                else:
                    raise Exception("NO ALIGNMENTS FOUND - SHOULD NOT BE HERE")


            # dump remaining reads from NR
            for nr_hit_dict, [nr_m8_dict, nr_m8_row] in zip(
                parse_tsv(nr_hitsummary2_tab, TAB_SCHEMA_MERGED, strict=False),
                parse_tsv(nr_m8, BLAST_OUTPUT_SCHEMA, raw_lines=True, strict=False)
            ):
                # assert files match
                assert nr_hit_dict['read_id'] == nr_m8_dict["qseqid"], f"Mismatched m8 and hit summary files for NR [{nr_hit_dict['read_id']} {nr_m8_dict['qseqid']}]."

                nr_alignment = nr_alignment_per_read.get(nr_hit_dict["read_id"])

                if nr_alignment:
                    output_m8.write(nr_m8_row)
                    nr_hit_dict["source_count_type"] = "NR"
                    self._write_tsv_row(nr_hit_dict, TAB_SCHEMA_MERGED, output_hit)

        # Create new merged m8 and hit summary files
        self.create_taxon_count_file(merged_m8_filename, merged_hit_filename, cluster_sizes_filename, merged_taxon_count_filename)

    def create_taxon_count_file(self, merged_m8_filename, merged_hit_filename, cluster_sizes_filename, merged_taxon_count_filename):
        # TOOO: Can this be consolidated throughout the pipeline?
        # This setup is mostly repeated in three steps. The list of taxa do not seem to change.
        count_type = 'merged_NT_NR'
        lineage_db = fetch_reference(
            self.additional_files["lineage_db"],
            self.ref_dir_local,
            allow_s3mi=False)
        deuterostome_db = None
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = fetch_reference(self.additional_files["deuterostome_db"],
                                                 self.ref_dir_local, allow_s3mi=False)  # Too small for s3mi
        taxon_whitelist = None
        if self.additional_attributes.get("use_taxon_whitelist"):
            taxon_whitelist = fetch_reference(self.additional_files.get("taxon_whitelist", DEFAULT_WHITELIST_S3),
                                                 self.ref_dir_local)
        blacklist_s3_file = self.additional_files.get('taxon_blacklist', DEFAULT_BLACKLIST_S3)
        taxon_blacklist = fetch_reference(blacklist_s3_file, self.ref_dir_local)
        cdhit_cluster_sizes_path = cluster_sizes_filename


        generate_taxon_count_json_from_m8(
            merged_m8_filename,
            merged_hit_filename,
            count_type,
            lineage_db,
            deuterostome_db,
            taxon_whitelist,
            taxon_blacklist,
            cdhit_cluster_sizes_path,
            merged_taxon_count_filename
        )

    def _write_tsv_row(self, data, schema, output):
        output.write("\t".join(str(data.get(k, "")) for k in schema.keys()) + "\n")

    def count_reads(self):
        pass
