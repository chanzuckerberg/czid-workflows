import csv
from collections import defaultdict

from idseq_dag.util.lineage import NULL_LINEAGE
from idseq_dag.util.dict import open_file_db_by_extension
from idseq_dag.util.m8 import build_should_keep_filter
from idseq_dag.util.parsing import HitSummaryMergedWriter
from idseq_dag.util.parsing import BlastnOutput6Reader, BlastnOutput6Writer
from idseq_dag.util.parsing import BlastnOutput6NTRerankedReader, BlastnOutput6NTRerankedWriter


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
    m8_reassigned_output_path: str,
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
            open(m8_reassigned_output_path, 'w') as m8_out_f, \
            open(hitsummary_output_path, 'w') as hitsummary_out_f:

        if db_type.lower() == "nt":
            m8_writer = BlastnOutput6NTRerankedWriter(m8_out_f)
        else:
            m8_writer = BlastnOutput6Writer(m8_out_f)

        hitsummary_writer = HitSummaryMergedWriter(hitsummary_out_f)

        if db_type.lower() == "nt":
            reader = BlastnOutput6NTRerankedReader(
                in_f,
                filter_invalid=True,
                min_alignment_length=min_alignment_length,
            )
        else:
            reader = BlastnOutput6Reader(in_f, min_alignment_length=min_alignment_length)

        hit_by_qseqid = {}
        for hit in reader:
            qseqid = hit["qseqid"]
            accession = hit["sseqid"]
            taxid = accession_to_taxid.get(accession.split(".")[0], "NA")
            lineage = taxid_to_lineage.get(str(taxid), NULL_LINEAGE)
            species_taxid, genus_taxid, family_taxid = lineage

            if not should_keep(lineage):
                continue

            if prev_hit := hit_by_qseqid.get(qseqid, False):
                # take the weighted mean
                prev_hit["pident"] = (prev_hit["pident"] * prev_hit["length"] + hit["pident"] * hit["length"]) / (prev_hit["length"] + hit["length"])
                prev_hit["evalue"] = (prev_hit["evalue"] * prev_hit["length"] + hit["evalue"] * hit["length"]) / (prev_hit["length"] + hit["length"])

                # take the sum for these values
                prev_hit["mismatch"] += hit["mismatch"]
                prev_hit["gapopen"] += hit["gapopen"]
                prev_hit["bitscore"] += hit["bitscore"]
                prev_hit["length"] += hit["length"]

                # take the min/max for these values
                prev_hit["qstart"] = min(prev_hit["qstart"], hit["qstart"])
                prev_hit["qend"] = max(prev_hit["qend"], hit["qend"])
                prev_hit["sstart"] = min(prev_hit["sstart"], hit["sstart"])
                prev_hit["send"] = max(prev_hit["send"], hit["send"])
            else:
                hit_by_qseqid[qseqid] = hit

        for hit in hit_by_qseqid.values():
            qseqid = hit["qseqid"]
            accession = hit["sseqid"]
            taxid = accession_to_taxid.get(accession.split(".")[0], "NA")
            lineage = taxid_to_lineage.get(str(taxid), NULL_LINEAGE)
            species_taxid, genus_taxid, family_taxid = lineage
            if qseqid not in contig_to_reads:
                m8_writer.writerow(hit)
                hitsummary_writer.writerow({
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
                hit = dict(**hit)
                hit["qseqid"] = read_id
                m8_writer.writerow(hit)
                hitsummary_writer.writerow({
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
