import json

from collections import namedtuple

import idseq_dag.util.log as log

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3
from idseq_dag.util.m8 import generate_taxon_count_json_from_m8
from idseq_dag.util.parsing import HitSummaryMergedReader, HitSummaryMergedWriter, BlastnOutput6NTRerankedReader, BlastnOutput6NTRerankedWriter


SpeciesAlignmentResults = namedtuple('SpeciesAlignmentResults', ['contig', 'read'])


def _create_taxon_count_file(
    merged_m8_filename,
    merged_hit_filename,
    merged_taxon_count_filename,
    lineage_db,
    deuterostome_db,
    taxon_whitelist,
    taxon_blacklist,
    cluster_sizes_path,
):
    # TOOO: Can this be consolidated throughout the pipeline?
    # This setup is mostly repeated in three steps. The list of taxa do not seem to change.
    count_type = 'merged_NT_NR'
    generate_taxon_count_json_from_m8(
        merged_m8_filename,
        merged_hit_filename,
        count_type,
        lineage_db,
        deuterostome_db,
        taxon_whitelist,
        taxon_blacklist,
        duplicate_cluster_sizes_path=cluster_sizes_path,
        output_json_file=merged_taxon_count_filename
    )


def _merge_taxon_counts(
    nt_m8,
    nt_hitsummary2_tab,
    nr_m8,
    nr_hitsummary2_tab,
    lineage_db,
    merged_m8_filename,
    merged_hit_filename,
    merged_taxon_count_filename,
    deuterostome_db,
    taxon_whitelist,
    taxon_blacklist,
    cluster_sizes_path,
):
    # Create new merged m8 and hit summary files
    nr_alignment_per_read = {}
    with open(nr_hitsummary2_tab) as nr_hit_summary_f:
        # if this is a bottleneck, consider
        # (1) if processing time bottleneck, load all the data to memory
        # (2) if memory bottleneck, going through nt first, since that will save us from storing
        #     results in memory for all the reads that get their hit from NT contigs
        for nr_hit_dict in HitSummaryMergedReader(nr_hit_summary_f):
            nr_alignment_per_read[nr_hit_dict["read_id"]] = SpeciesAlignmentResults(
                contig=nr_hit_dict.get("contig_species_taxid"),
                read=nr_hit_dict.get("species_taxid"),
            )

    with open(merged_m8_filename, "w") as output_blastn_6_f, open(merged_hit_filename, "w") as output_hit_summary_f:
        output_blastn_6_writer = BlastnOutput6NTRerankedWriter(output_blastn_6_f)
        output_hit_summary_writer = HitSummaryMergedWriter(output_hit_summary_f)

        with open(nt_m8) as input_nt_blastn_6_f, open(nt_hitsummary2_tab) as input_nt_hit_summary_f:
            # first pass for NR and output to m8 files if assignment should come from NT
            for nt_hit_dict, nt_m8_dict in zip(
                HitSummaryMergedReader(input_nt_hit_summary_f),
                BlastnOutput6NTRerankedReader(input_nt_blastn_6_f),
            ):
                # assert files match
                assert nt_hit_dict['read_id'] == nt_m8_dict["qseqid"], f"Mismatched m8 and hit summary files for nt [{nt_hit_dict['read_id']} != {nt_m8_dict['qseqid']}]"

                nr_alignment = nr_alignment_per_read.get(nt_hit_dict["read_id"])
                nt_alignment = SpeciesAlignmentResults(
                    contig=nt_hit_dict.get("contig_species_taxid"),
                    read=nt_hit_dict.get("species_taxid")
                )
                has_nt_contig_hit = nt_alignment.contig
                has_nr_contig_hit = nr_alignment and nr_alignment.contig
                has_nt_read_hit = nt_alignment.read
                has_nr_read_hit = nr_alignment and nr_alignment.read
                if has_nt_contig_hit or (not has_nr_contig_hit and has_nt_read_hit):
                    output_blastn_6_writer.writerow(nt_m8_dict)
                    nt_hit_dict["source_count_type"] = "NT"
                    output_hit_summary_writer.writerow(nt_hit_dict)
                    if nr_alignment:
                        del nr_alignment_per_read[nt_hit_dict["read_id"]]
                elif has_nr_contig_hit or has_nr_read_hit:
                    continue
                else:
                    raise Exception("NO ALIGNMENTS FOUND - Should not be here")

        with open(nr_m8) as input_nr_blastn_6_f, open(nr_hitsummary2_tab) as input_nr_hit_summary_f:
            # dump remaining reads from NR
            for nr_hit_dict, nr_m8_dict in zip(
                HitSummaryMergedReader(input_nr_hit_summary_f),
                BlastnOutput6NTRerankedReader(input_nr_blastn_6_f),
            ):
                # assert files match
                assert nr_hit_dict['read_id'] == nr_m8_dict["qseqid"], f"Mismatched m8 and hit summary files for NR [{nr_hit_dict['read_id']} {nr_m8_dict['qseqid']}]."

                nr_alignment = nr_alignment_per_read.get(nr_hit_dict["read_id"])

                if nr_alignment:
                    output_blastn_6_writer.writerow(nr_m8_dict)
                    nr_hit_dict["source_count_type"] = "NR"
                    output_hit_summary_writer.writerow(nr_hit_dict)

    # Create new merged m8 and hit summary files
    _create_taxon_count_file(
        merged_m8_filename,
        merged_hit_filename,
        merged_taxon_count_filename,
        lineage_db,
        deuterostome_db,
        taxon_whitelist,
        taxon_blacklist,
        cluster_sizes_path,
    )


def _merge_contigs(nt_contig_summary_json, nr_contig_summary_json, merged_contig_summary_json):
    with open(nt_contig_summary_json) as f:
        nt_contigs = json.load(f)
    with open(nr_contig_summary_json) as f:
        nr_contigs = json.load(f)

    merged_contigs = {}
    nt_contigs_name_set = set()
    for taxid_contigs in nt_contigs:
        taxid_contigs['count_type'] = 'merged_NT_NR'
        merged_contigs[taxid_contigs['taxid']] = taxid_contigs
        # save to set to keep track of contigs that aligned to NT
        nt_contigs_name_set |= set(taxid_contigs['contig_counts'].keys())

    for contigs_per_taxid in nr_contigs:
        # remove contigs that aligned to NT
        taxid = contigs_per_taxid['taxid']
        contigs_per_taxid['contig_counts'] = {
            contig_name: contig_counts
            for contig_name, contig_counts in contigs_per_taxid['contig_counts'].items()
            if contig_name not in nt_contigs_name_set
        }
        if len(contigs_per_taxid['contig_counts']) > 0:
            if (taxid in merged_contigs):
                merged_contigs[taxid]['contig_counts'].update(contigs_per_taxid['contig_counts'])
            else:
                contigs_per_taxid['count_type'] = 'merged_NT_NR'
                merged_contigs[taxid] = contigs_per_taxid

    with open(merged_contig_summary_json, 'w') as output_contig_json:
        json.dump(list(merged_contigs.values()), output_contig_json)


def compute_merged_taxon_counts(
    nt_m8,
    nt_hitsummary2_tab,
    nt_contig_summary_json,
    nr_m8,
    nr_hitsummary2_tab,
    nr_contig_summary_json,
    lineage_db,
    deuterostome_db,
    taxon_whitelist,
    taxon_blacklist,
    merged_m8_filename,
    merged_hit_filename,
    merged_taxon_count_filename,
    merged_contig_summary_json,
    cluster_sizes_path=None,
):
    _merge_taxon_counts(
        nt_m8,
        nt_hitsummary2_tab,
        nr_m8,
        nr_hitsummary2_tab,
        lineage_db,
        merged_m8_filename,
        merged_hit_filename,
        merged_taxon_count_filename,
        deuterostome_db,
        taxon_whitelist,
        taxon_blacklist,
        cluster_sizes_path,
    )
    _merge_contigs(nt_contig_summary_json, nr_contig_summary_json, merged_contig_summary_json)


ComputeMergedTaxonCountsInputs = namedtuple('ComputeMergedTaxonCountsInputs', ['nt_m8', 'nt_hitsummary2_tab', 'nt_contig_summary_json', 'nr_m8', 'nr_hitsummary2_tab', 'nr_contig_summary_json', 'cluster_sizes_filename'])
ComputeMergedTaxonCountsOutputs = namedtuple('ComputeMergedTaxonCountsOutputs', ['merged_m8_filename', 'merged_hit_filename', 'merged_taxon_count_filename', 'merged_contig_summary_json'])


class ComputeMergedTaxonCounts(PipelineStep):
    '''
    Merge taxon results from NT database and from NR database.
    Create m8 and hit summary files for the new merged classifier and taxon count and contig json files.
    '''

    def run(self):
        # Unpack inputs and outputs
        self.inputs = ComputeMergedTaxonCountsInputs(*self.input_files_local)
        self.outputs = ComputeMergedTaxonCountsOutputs(*self.output_files_local())

        _merge_taxon_counts(
            self.inputs.nt_m8,
            self.inputs.nt_hitsummary2_tab,
            self.inputs.nr_m8,
            self.inputs.nr_hitsummary2_tab,
            self.additional_files["lineage_db"],
            self.outputs.merged_m8_filename,
            self.outputs.merged_hit_filename,
            self.outputs.merged_taxon_count_filename,
            self.additional_files.get("deuterostome_db", None),
            self.additional_files.get("taxon_whitelist", None),
            self.additional_files.get("taxon_blacklist", DEFAULT_BLACKLIST_S3),
            self.inputs.cluster_sizes_filename,
        )
        _merge_contigs(
            self.inputs.nt_contig_summary_json,
            self.inputs.nr_contig_summary_json,
            self.outputs.merged_contig_summary_json,
        )

    def count_reads(self):
        pass
