TAB_SCHEMA = {
    "read_id": str,
    "level": int,
    "taxid": int,
    "accession_id": str,
    "species_taxid": int,
    "genus_taxid": int,
    "family_taxid": int,
    "contig_id": str,
    "contig_accession_id": str,
    "contig_species_taxid": int,
    "contig_genus_taxid": int,
    "contig_family_taxid": int,
    "from_assembly" : bool,
}

TAB_SCHEMA_MERGED = dict(TAB_SCHEMA, **{
    "source_count_type": str
})
