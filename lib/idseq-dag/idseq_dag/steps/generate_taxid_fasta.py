from typing import List
from shutil import copy as file_copy
from idseq_dag.util.parsing import HitSummaryMergedReader
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.lineage as lineage

from idseq_dag.util.dict import open_file_db_by_extension
from idseq_dag.util import fasta


def _get_valid_lineage(hits_by_read_id, lineage_map, read_id: str) -> List[str]:
    # If the read aligned to something, then it would be present in the
    # summary file for count type, and correspondingly in hits_by_read_id
    # even if the hits disagree so much that the
    # "valid_hits" entry is just ("-1", "-1"). If the read didn't align
    # to anything, we also represent that with ("-1", "-1"). This ("-1",
    # "-1") gets translated to NULL_LINEAGE.
    hit_taxid, hit_level = hits_by_read_id.get(read_id, (-1, -1))
    hit_lineage = lineage_map.get(str(hit_taxid), lineage.NULL_LINEAGE)
    return lineage.validate_taxid_lineage(hit_lineage, hit_taxid, hit_level)


def generate_taxid_fasta(
    input_fa_name: str,
    nt_hit_summary_path: str,
    nr_hit_summary_path: str,
    lineage_db_path: str,
    output_fasta_path: str,
):
    with open(nt_hit_summary_path) as nt_hit_summary_f, open(nr_hit_summary_path) as nr_hit_summary_f:
        nr_hits_by_read_id = {row["read_id"]: (row["taxid"], row["level"]) for row in HitSummaryMergedReader(nr_hit_summary_f)}
        nt_hits_by_read_id = {row["read_id"]: (row["taxid"], row["level"]) for row in HitSummaryMergedReader(nt_hit_summary_f)}

    with open(output_fasta_path, "w") as output_fa, \
            open_file_db_by_extension(lineage_db_path, "lll") as lineage_map:  # noqa
        for read in fasta.iterator(input_fa_name):
            # Example read_id: "NR::NT:CP010376.2:NB501961:14:HM7TLBGX2:1:23109
            # :12720:8743/2"
            # Translate the read information into our custom format with fake
            # taxids at non-specific hit levels.
            # TODO: (tmorse) fasta parsing
            annotated_read_id = read.header.lstrip('>')
            read_id = annotated_read_id.split(":", 4)[-1]

            nr_taxid_species, nr_taxid_genus, nr_taxid_family = _get_valid_lineage(
                nr_hits_by_read_id, lineage_map, read_id)
            nt_taxid_species, nt_taxid_genus, nt_taxid_family = _get_valid_lineage(
                nt_hits_by_read_id, lineage_map, read_id)

            fields = ["family_nr", nr_taxid_family, "family_nt", nt_taxid_family]
            fields += ["genus_nr", nr_taxid_genus, "genus_nt", nt_taxid_genus]
            fields += ["species_nr", nr_taxid_species, "species_nt", nt_taxid_species]
            fields += [annotated_read_id]
            new_read_name = ('>' + ':'.join(fields) + '\n')

            output_fa.write(new_read_name)
            output_fa.write(read.sequence + "\n")


CONFORMING_PREAMBLE = ">family_nr:-300:family_nt:-300:genus_nr:-200:genus_nt:-200:species_nr:-100:species_nt:-100:"
def conform_unmapped_read_header(header: str):
    """Converts fasta headers from unmapped reads to match structure of mapped.

    Implementation is very, very tightly coupled with how we output unmapped
    (AKA, unidentified) reads from GenerateAnnotatedFasta task. This function
    takes the headers of those unmapped reads and converts them to match the
    structure we expect for all the mapped reads, eg where `family_nr` etc
    are all explicitly given. See above `CONFORMING_PREAMBLE` for what that
    looks like. Result is a header that looks like this:
    >family_nr:-300:...:species_nt:-100:NR::NT::ORIGINAL_RAW_ID_GOES_HERE

    See `generate_fasta_with_unmapped_included` below for why we do this."""
    return CONFORMING_PREAMBLE + header.lstrip('>')


def generate_fasta_with_unmapped_included(
    mapped_fa_path: str,
    unmapped_fa_path: str,
    output_with_unmapped_path: str,
):
    """Creates additional fasta that has both mapped and unmapped reads.

    This takes an already existing fasta (expected to be the mapped reads with
    taxids in place), copies it, then appends the unmapped (AKA, unidentified)
    reads to it while conforming the headers of those unmapped reads so they
    match up to the same structure used for the taxid+accession mapped reads.

    Intent: while some steps need only the mapped reads to work properly, other
    steps want to have all the non-host reads, both mapped and unmapped. That
    is also what our users generally expect when downloading a sample's reads.
    This provides another fasta that has **all** the non-host reads, but with
    the headers all in a consistent structure so that any downstream parsing
    won't blow up because of varying header structures (and so long as those
    downstream steps aren't very strongly assuming just getting mapped reads).
    """
    file_copy(mapped_fa_path, output_with_unmapped_path)
    with open(output_with_unmapped_path, "a") as output_fa:
        for read in fasta.iterator(unmapped_fa_path):
            conformed_header = conform_unmapped_read_header(read.header)
            output_fa.write(conformed_header + "\n")
            output_fa.write(read.sequence + "\n")


class PipelineStepGenerateTaxidFasta(PipelineStep):
    """Generate taxid FASTA from hit summaries. Intermediate conversion step
    that includes handling of non-specific hits with artificial tax_ids.
    """

    def run(self):
        input_fa_name = self.input_files_local[0][0]
        # We determine if intended for use in `short-read-mngs/postprocess.wdl`
        # or `short-read-mngs/experimental.wdl` by shape of input files list.
        is_structured_as_postprocess_call = (len(self.input_files_local) > 1)
        if is_structured_as_postprocess_call:  # short-read-mngs/postprocess.wdl
            input_fa_name = self.input_files_local[0][0]
            nt_hit_summary_path, nr_hit_summary_path = self.input_files_local[1][2], self.input_files_local[2][2]
        else:  # short-read-mngs/experimental.wdl
            input_fa_name = self.input_files_local[0][0]
            nt_hit_summary_path, nr_hit_summary_path = self.input_files_local[0][1], self.input_files_local[0][2]

        generate_taxid_fasta(
            input_fa_name,
            nt_hit_summary_path,
            nr_hit_summary_path,
            self.additional_files["lineage_db"],
            self.output_files_local()[0],
        )

        if is_structured_as_postprocess_call:
            # For short-read-mngs/postprocess.wdl ONLY we generate additional
            # FASTA that has both the mapped reads we just did and unmappeds.
            input_unidentified_fa_name = self.input_files_local[0][1]
            generate_fasta_with_unmapped_included(
                self.output_files_local()[0],
                input_unidentified_fa_name,
                self.output_files_local()[1],
            )
