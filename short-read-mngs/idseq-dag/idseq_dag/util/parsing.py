
import csv

from abc import ABC, abstractstaticmethod
from idseq_dag.util.m8 import MAX_EVALUE_THRESHOLD
from typing import Any, Callable, Dict, Iterable, List, TextIO, Tuple

Schema = List[Tuple[str, type]]

class _TSVWithSchmaBase(ABC):
    """
    The _TSVWithSchemaBase class is a base class for reading and writing TSVs without headers with schemas provided by subclasses

    This class has functionality to support several variants on a schema and select one based on the data presented. To use this
    class it must be extended and the `_variants` method must be defined to return a list of named partial schemas. The first
    partial schema becomes the first complete schema and each subsequent schema is appended to the previous complete schema to make
    a complete schema. This guaruntees that there is exactly one schema for each number of fields field count so the schemas can be
    recognized. For reading, the variant is recognized based on how many fields were encountered. For writing the smallest possible
    variant is selected that contains all of the provided fields. This schema recognition is necessary because previously we were using
    optional fields to capture a finite set of different variants.

    Args:
        tsv_stream (TextIO): a text stream to write to/read from
        *schemas (List[Tuple[str, Callable[[str], Any]]]): partial schemas
    """

    @abstractstaticmethod
    def _variants() -> List[Tuple[str, Schema]]:
        pass

    @classmethod
    def _build_variant_map(cls) -> Tuple[Dict[str, Schema], Dict[int, str], Dict[str, str]]:
        schemas = cls._variants()
        assert schemas, "_TSVWithSchemaReader requires at least one schema"
        n = 0
        name_to_variant = num_fields_to_variant = field_to_first_variant = {}
        for i, (variant, schema) in enumerate(schemas):
            assert schema, "_TSVWithSchemaBase does not support empty schemas"
            n += len(schema)
            assert variant not in name_to_variant, f"_TSVWithSchemaBase encountered duplicate variant name: {variant}"
            name_to_variant[variant] = (name_to_variant[num_fields_to_variant[n - len(schema)]] if i > 0 else []) + schema
            num_fields_to_variant[n] = variant
            for field, _ in schema:
                assert field not in field_to_first_variant, f"_TSVWithSchemaBase encountered duplicate field name: {field}"
                field_to_first_variant[field] = variant
        return name_to_variant, num_fields_to_variant, field_to_first_variant

    @classmethod
    def fields(cls, variant) -> List[str]:
        return [k for k, _ in cls._build_variant_map()[0][variant]]

    def __init__(self, tsv_stream: TextIO) -> None:
        self._tsv_stream = tsv_stream
        self._name_to_variant, self._num_fields_to_variant, self._field_to_first_variant = self._build_variant_map()
        self._schema = None

class _TSVWithSchemaReader(_TSVWithSchmaBase, ABC):
    def __init__(self, tsv_stream: TextIO, *schemas: List[Tuple[str, Callable[[str], Any]]]) -> None:
        super().__init__(tsv_stream, *schemas)
        self._generator = ({
            key: _type(value) for ((key, _type), value) in zip(self._get_schema(row), row)
        } for row in csv.reader(self._tsv_stream, delimiter="\t"))

    def __iter__(self):
        return self

    def __next__(self, *args) -> Dict[str, Any]:
        return next(self._generator, *args)

    def _get_schema(self, row):
        if self._schema:
            return self._schema
        if len(row) not in self._num_fields_to_variant:
            raise Exception(f"_TSVWithSchemaReader error. Input: \"{row}\" has {len(row)} fields, no associated schema found in {self._name_to_variant}")
        return self._name_to_variant[self._num_fields_to_variant[len(row)]]

class _TSVWithSchemaWriter(_TSVWithSchmaBase, ABC):
    def __init__(self, tsv_stream: TextIO, *schemas: List[Tuple[str, Callable[[str], Any]]]) -> None:
        super().__init__(tsv_stream, *schemas)
        self._writer = csv.writer(tsv_stream, delimiter="\t")

    def _get_schema(self, row):
        if self._schema:
            return self._schema
        first_inclusive_variant = max(
            (self._field_to_first_variant.get(field, '') for field in row.keys()),
            key=lambda variant: len(self._name_to_variant[variant] if variant else 0)
        )
        if first_inclusive_variant:
            raise Exception(f"_TSVWithSchemaWriter error. Input: \"{row}\" has {len(row)} fields, no associated schema or common fields found in {self._name_to_variant}")
        return self._name_to_variant[first_inclusive_variant]

    def _dict_row_to_list(self, row: Dict[str, Any]) -> List[Any]:
        return [row.get(key) for (key, _) in self._get_schema(row)]

    def write_all(self, rows: Iterable[Dict[str, Any]]) -> None:
        self._writer.writerows(self._dict_row_to_list(row) for row in rows)

    def write(self, row: Dict[str, Any]) -> None:
        self._writer.writerow(self._dict_row_to_list(row))


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

class _BlastnOutput6Base(ABC):
    @staticmethod
    def _variants():
        return [
            ("base", _BLAST_OUTPUT_SCHEMA),
            ("nt", _BLAST_OUTPUT_NT_SCHEMA),
            ("reranked_nt", _RERANKED_BLAST_OUTPUT_NT_SCHEMA),
        ]

class BlastnOutput6Reader(_BlastnOutput6Base, _TSVWithSchemaReader):
    def __init__(self, tsv_stream: TextIO, filter_invalid: bool = False, min_alignment_length: int = 0):
        if filter_invalid:
            # The output of rapsearch2 contains comments that start with '#', these should be skipped
            super().__init__(line for line in tsv_stream if line and line[0] != "#")
            self._generator = (row for row in self._generator if self.row_is_valid(row, min_alignment_length))
        else:
            super().__init__(tsv_stream)

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
            row["length"] >= min_alignment_length,
            -0.25 < row["pident"] < 100.25,
            row["evalue"] == row["evalue"],
            row["evalue"] <= MAX_EVALUE_THRESHOLD,
        ])

class BlastnOutput6Writer(_BlastnOutput6Base, _TSVWithSchemaWriter):
    pass


_HIT_SUMMARY_SCHEMA = [
    ("read_id", str),
    ("level", int),
    ("taxid", int),
    ("accession_id", str),
    ("species_taxid", int),
    ("genus_taxid", int),
    ("family_taxid", int),
]

_HIT_SUMMARY_CONTIG_SCHEMA = [
    ("contig_id", str),
    ("contig_accession_id", str),
    ("contig_species_taxid", int),
    ("contig_genus_taxid", int),
    ("contig_family_taxid", int),
]

_HIT_SUMMARY_ASSEMBLY_SOURCE_SCHEMA = [
    ("from_assembly", str),
]

_HIT_SUMMARY_MERGED_SCHEMA = [
    ("source_count_type", str),
]

class _HitSummaryBase(ABC):
    @staticmethod
    def _variants():
        return [
            ("base", _HIT_SUMMARY_SCHEMA),
            ("contig", _HIT_SUMMARY_CONTIG_SCHEMA),
            ("assembly_source", _HIT_SUMMARY_ASSEMBLY_SOURCE_SCHEMA),
            ("merged", _HIT_SUMMARY_MERGED_SCHEMA),
        ]

class HitSummaryReader(_HitSummaryBase, _TSVWithSchemaReader):
    pass

class HitSummaryWriter(_HitSummaryBase, _TSVWithSchemaWriter):
    pass
