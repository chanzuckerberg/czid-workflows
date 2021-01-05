
from csv import DictReader, DictWriter
from typing import Any, Iterable, Optional, Sequence, Text, TextIO, Tuple

# Alignments with e-values greater than 1 are low-quality alignments and associated with
# a high rate of false-positives. These should be filtered at all alignment steps.
MAX_EVALUE_THRESHOLD = 1


class _TypedDictTSVReader(DictReader):
    """
    Similar to DictReader but instead of a sequence of fieldnames it takes a sequence of
    tuples, the first element being the field name and the second being the field's type.
    After reading a row, this class will convert each element to it's associated type.
    """
    def __init__(self, f: Iterable[Text], schema: Sequence[Tuple[str, type]]) -> None:
        fieldnames = [field for field, _ in schema]
        self._types = {field: _type for field, _type in schema}
        super().__init__(f, fieldnames=fieldnames, delimiter="\t")

    def __next__(self):
        row = super().__next__()
        assert len(row) <= len(self._types), f"row {row} contains fields not in schema {self._types}"
        for key, value in row.items():
            if value:
                row[key] = self._types[key](value)
        return row


class _TypedDictTSVWriter(DictWriter):
    """
    This is just a convenience class so you don't need to pull out the field names
    """
    def __init__(self, f: Any, schema: Sequence[Tuple[str]]) -> None:
        fieldnames = [field for field, _ in schema]
        super().__init__(f, fieldnames, delimiter="\t")


# blastn output format 6 as documented in
# http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# it's also the format of our GSNAP and RAPSEARCH2 output
_BLASTN_OUTPUT_6_SCHEMA = [
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
_BLASTN_OUTPUT_6_NT_SCHEMA = _BLASTN_OUTPUT_6_SCHEMA + [
    ("qlen", int),      # query sequence length, helpful for computing qcov
    ("slen", int),      # subject sequence length, so far unused in IDseq
]

# This is needed to pass the output fields to blastn on NT (idseq-dag/idseq_dag/steps/blast_contigs.py:441)
BLASTN_OUTPUT_6_NT_FIELDS = [field for field, _ in _BLASTN_OUTPUT_6_NT_SCHEMA]

# Re-ranked output of blastn.  One row per query.  Two additional columns.
_BLASTN_OUTPUT_6_NT_RERANKED_SCHEMA = _BLASTN_OUTPUT_6_NT_SCHEMA + [
    ("qcov", float),     # fraction of query covered by the optimal set of HSPs
    ("hsp_count", int),   # cardihnality of optimal fragment cover;  see BlastCandidate
]


class _BlastnOutput6ReaderBase(_TypedDictTSVReader):
    """
    This class is a bit of an oddball due to some compatibility concerns. In addition
    to the normal parsing stuff it also:
    1. Ignores comments (lines starting with '#') in tsv files, rapsearch2 adds them
    2. Supports filtering rows that we consider invalid
    """
    def __init__(self, f: Iterable[Text], schema: Sequence[Tuple[str, type]], filter_invalid: bool = False, min_alignment_length: int = 0):
        self._filter_invalid = filter_invalid
        self._min_alignment_length = min_alignment_length

        # The output of rapsearch2 contains comments that start with '#', these should be skipped
        filtered_stream = (line for line in f if not line.startswith("#"))
        super().__init__(filtered_stream, schema,)

    def __next__(self):
        if not self._filter_invalid:
            return super().__next__()
        row = super().__next__()
        while row and not self._row_is_valid(row):
            row = super().__next__()
        return row

    def _row_is_valid(self, row) -> bool:
        # GSNAP outputs bogus alignments (non-positive length /
        # impossible percent identity / NaN e-value) sometimes,
        # and usually they are not the only assignment, so rather than
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
            row["length"] >= self._min_alignment_length,
            -0.25 < row["pident"] < 100.25,
            row["evalue"] == row["evalue"],
            row["evalue"] <= MAX_EVALUE_THRESHOLD,
        ])


class BlastnOutput6Reader(_BlastnOutput6ReaderBase):
    """
    This class is a bit of an oddball due to some compatibility concerns. In addition
    to the normal parsing stuff it also:
    1. Ignores comments (lines starting with '#') in tsv files, rapsearch2 adds them
    2. Supports filtering rows that we consider invalid
    """
    def __init__(self, f: Iterable[Text], filter_invalid: bool = False, min_alignment_length: int = 0):
        super().__init__(f, _BLASTN_OUTPUT_6_SCHEMA, filter_invalid, min_alignment_length)


class BlastnOutput6Writer(_TypedDictTSVWriter):
    def __init__(self, f: Any) -> None:
        super().__init__(f, _BLASTN_OUTPUT_6_SCHEMA)


class BlastnOutput6NTReader(_BlastnOutput6ReaderBase):
    def __init__(self, f: Iterable[Text], filter_invalid: bool = False, min_alignment_length: int = 0):
        super().__init__(f, _BLASTN_OUTPUT_6_NT_SCHEMA, filter_invalid, min_alignment_length)


class BlastnOutput6NTWriter(_TypedDictTSVWriter):
    def __init__(self, f: Any) -> None:
        super().__init__(f, _BLASTN_OUTPUT_6_NT_SCHEMA)


class BlastnOutput6NTRerankedReader(_BlastnOutput6ReaderBase):
    def __init__(self, f: Iterable[Text], filter_invalid: bool = False, min_alignment_length: int = 0):
        super().__init__(f, _BLASTN_OUTPUT_6_NT_RERANKED_SCHEMA, filter_invalid, min_alignment_length)


class BlastnOutput6NTRerankedWriter(_TypedDictTSVWriter):
    def __init__(self, f: Any) -> None:
        super().__init__(f, _BLASTN_OUTPUT_6_NT_RERANKED_SCHEMA)


_HIT_SUMMARY_SCHEMA = [
    ("read_id", str),
    ("level", int),
    ("taxid", int),
    ("accession_id", str),
    ("species_taxid", int),
    ("genus_taxid", int),
    ("family_taxid", int),
]


_HIT_SUMMARY_MERGED_SCHEMA = _HIT_SUMMARY_SCHEMA + [
    ("contig_id", str),
    ("contig_accession_id", str),
    ("contig_species_taxid", int),
    ("contig_genus_taxid", int),
    ("contig_family_taxid", int),
    ("from_assembly", str),
    ("source_count_type", str),
]


class HitSummaryReader(_TypedDictTSVReader):
    def __init__(self, f: Iterable[Text]) -> None:
        super().__init__(f, _HIT_SUMMARY_SCHEMA)


class HitSummaryWriter(_TypedDictTSVWriter):
    def __init__(self, f: Any) -> None:
        super().__init__(f, _HIT_SUMMARY_SCHEMA)


class HitSummaryMergedReader(_TypedDictTSVReader):
    def __init__(self, f: Iterable[Text]) -> None:
        super().__init__(f, _HIT_SUMMARY_MERGED_SCHEMA)


class HitSummaryMergedWriter(_TypedDictTSVWriter):
    def __init__(self, f: Any) -> None:
        super().__init__(f, _HIT_SUMMARY_MERGED_SCHEMA)
