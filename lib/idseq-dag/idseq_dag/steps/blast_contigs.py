import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log
import idseq_dag.util.m8 as m8
import idseq_dag.util.s3 as s3
import json
import multiprocessing
import os

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.trace_lock import TraceLock

from idseq_dag.steps.run_assembly import generate_info_from_sam
from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3, DEFAULT_WHITELIST_S3
from idseq_dag.util.parsing import BlastnOutput6Reader
from idseq_dag.util.parsing import BlastnOutput6NTReader
from idseq_dag.util.parsing import BlastnOutput6NTRerankedReader, BlastnOutput6NTRerankedWriter
from idseq_dag.util.parsing import HitSummaryReader, HitSummaryMergedWriter
from idseq_dag.util.top_hits import get_top_m8_nr, get_top_m8_nt
from idseq_dag.util.taxon_summary import generate_taxon_summary


_MIN_REF_FASTA_SIZE = 25
_MIN_ASSEMBLED_CONTIG_SIZE = 25


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
    -evalue 1e-10
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
            os.path.getsize(assembled_contig) < _MIN_ASSEMBLED_CONTIG_SIZE or
            os.path.getsize(reference_fasta) < _MIN_REF_FASTA_SIZE)

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
        generate_info_from_sam(bowtie_sam, read2contig, duplicate_cluster_sizes_path=duplicate_cluster_sizes_path)

        (updated_read_dict, read2blastm8, contig2lineage, added_reads) = self.update_read_dict(
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
                                                     duplicate_cluster_sizes_path=duplicate_cluster_sizes_path, output_json_file=refined_counts_with_dcr)

        # generate contig stats at genus/species level
        with log.log_context("PipelineStepBlastContigs", {"substep": "generate_taxon_summary"}):
            contig_taxon_summary = generate_taxon_summary(
                read2contig,
                contig2lineage,
                updated_read_dict,
                added_reads,
                db_type,
                # same filter as applied in generate_taxon_count_json_from_m8
                m8.build_should_keep_filter(deuterostome_db, taxon_whitelist, taxon_blacklist),
                duplicate_cluster_sizes_path
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
    def generate_m8_and_hit_summary(consolidated_dict, added_reads, read2blastm8,
                                    hit_summary_path, deduped_blastn_6_path,
                                    refined_hit_summary_path, refined_blastn_6_path):
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

    @staticmethod
    def update_read_dict(read2contig, blast_top_blastn_6_path, read_dict, accession_dict, db_type):
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
        get_top_m8_nt(blast_m8, blast_top_m8)

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
                    "-evalue",
                    1e-10,
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
        get_top_m8_nr(blast_m8, blast_top_m8)
