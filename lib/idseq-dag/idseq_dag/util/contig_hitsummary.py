import csv
from collections import defaultdict

from lineage import NULL_LINEAGE
from dict import open_file_db_by_extension
from m8 import build_should_keep_filter
from parsing import HitSummaryMergedWriter
from parsing import BlastnOutput6Reader
from parsing import BlastnOutput6NTRerankedReader


def summarize_hits(
    blast6_path: str,
    db_type: str,
    read_to_contig_tsv_path: str,
    deuterostome_path: str,
    taxon_allow_list_path: str,
    taxon_ignore_list_path: str,
    taxid_to_lineage_path: str,
    accession_to_taxid_path: str,
    min_alignment_length: int,
    hitsummary_output_path: str,
):
    with open(read_to_contig_tsv_path) as f:
        contig_to_reads = defaultdict(list)
        for row in csv.reader(f, delimiter='\t'):
            contig_to_reads[row[1]].append(row[0])

    should_keep = build_should_keep_filter(
        deuterostome_path, taxon_allow_list_path, taxon_ignore_list_path)

    with open(blast6_path) as in_f, \
        open_file_db_by_extension(accession_to_taxid_path, "L") as accession_to_taxid, \
        open_file_db_by_extension(taxid_to_lineage_path, "lll") as taxid_to_lineage, \
    open(hitsummary_output_path, 'w') as out_f:

        writer = HitSummaryMergedWriter(out_f)

        if db_type.lower() == "nt":
            reader = BlastnOutput6NTRerankedReader(
                in_f,
                filter_invalid=True,
                min_alignment_length=min_alignment_length,
            )
        else:
            reader = BlastnOutput6Reader(in_f, min_alignment_length=min_alignment_length)

        for hit in reader:
            qseqid = hit["qseqid"]
            accession = hit["sseqid"]
            taxid = accession_to_taxid.get(accession)
            lineage = taxid_to_lineage.get(str(taxid), NULL_LINEAGE)
            species_taxid, genus_taxid, family_taxid = lineage

            if not should_keep(lineage):
                continue

            if qseqid not in contig_to_reads:
                writer.writerow({
                    "read_id": qseqid,
                    "level": 1,  # NOTE: this is always 1 regardless of the actual level
                    "taxid": species_taxid,
                    "accession_id": accession,
                    "species_taxid": species_taxid,
                    "genus_taxid": genus_taxid,
                    "family_taxid": family_taxid,
                    "contig_id": "*",
                    "contig_accession_id": accession,
                    "contig_species_taxid": species_taxid,
                    "contig_genus_taxid": genus_taxid,
                    "contig_family_taxid": family_taxid,
                })
                continue

            for read_id in contig_to_reads[qseqid]:
                writer.writerow({
                    "read_id": read_id,
                    "level": 1,  # NOTE: this is always 1 regardless of the actual level
                    "taxid": species_taxid,
                    "accession_id": accession,
                    "species_taxid": species_taxid,
                    "genus_taxid": genus_taxid,
                    "family_taxid": family_taxid,
                    "contig_id": qseqid,
                    "contig_accession_id": accession,
                    "contig_species_taxid": species_taxid,
                    "contig_genus_taxid": genus_taxid,
                    "contig_family_taxid": family_taxid,
                    "from_assembly": "from_assembly",
                })
