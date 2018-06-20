from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.taxid_lineage as taxid_lineage
import idseq_dag.util.s3 as s3
import shelve


class PipelineStepGenerateTaxidFasta(PipelineStep):
    """Generate taxid FASTA from hit summaries. Intermediate conversion step
    that includes handling of non-specific hits with artificial tax_ids.
    """

    def run(self):
        # Setup
        input_files = self.input_files_local[0]
        hit_summary_files = {'NT': input_files[2], 'NR': input_files[3]}

        # Open lineage db
        lineage_db = s3.fetch_from_s3(
            self.additional_files["lineage_db"],
            self.ref_dir_local,
            allow_s3mi=True)
        lineage_map = shelve.open(lineage_db.replace(".db", ""))

        # Get primary hit mappings
        valid_hits = PipelineStepGenerateTaxidFasta.parse_hits(
            hit_summary_files)

        input_fa = open(input_files[0], 'rb')
        output_fa = open(self.output_files_local()[0], 'wb')
        seq_name = input_fa.readline()
        seq_data = input_fa.readline()
        while len(seq_name) > 0 and len(seq_data) > 0:
            # Example read_id: "NR::NT:CP010376.2:NB501961:14:HM7TLBGX2:1:23109
            # :12720:8743/2"
            # Translate the read information into our custom format with fake
            # taxids at non-specific hit levels.
            annotated_read_id = seq_name.decode("utf-8").rstrip().lstrip('>')
            read_id = annotated_read_id.split(":", 4)[-1]

            nr_taxid_species, nr_taxid_genus, nr_taxid_family = PipelineStepGenerateTaxidFasta.get_valid_lineage(
                valid_hits, lineage_map, read_id, 'NR')
            nt_taxid_species, nt_taxid_genus, nt_taxid_family = PipelineStepGenerateTaxidFasta.get_valid_lineage(
                valid_hits, lineage_map, read_id, 'NT')

            fields = ["family_nr", nr_taxid_family, "family_nt", nt_taxid_family]
            fields += ["genus_nr", nr_taxid_genus, "genus_nt", nt_taxid_genus]
            fields += ["species_nr", nr_taxid_species, "species_nt", nt_taxid_species]
            fields += [annotated_read_id]
            new_read_name = ('>' + ':'.join(fields) + '\n').encode()

            output_fa.write(new_read_name)
            output_fa.write(seq_data)
            seq_name = input_fa.readline()
            seq_data = input_fa.readline()
        input_fa.close()
        output_fa.close()

    def count_reads(self):
        pass

    @staticmethod
    def parse_hits(hit_summary_files):
        """Return map of {NT, NR} => read_id => (hit_taxid_str, hit_level_str)"""
        valid_hits = {}
        for count_type, summary_file in hit_summary_files.items():
            hits = {}
            with open(summary_file, "rb") as sf:
                for hit_line in sf:
                    hit_line_columns = hit_line.decode("utf-8").strip().split(
                        "\t")
                    if len(hit_line_columns) >= 3:
                        hit_read_id = hit_line_columns[0]
                        hit_level_str = hit_line_columns[1]
                        hit_taxid_str = hit_line_columns[2]
                        hits[hit_read_id] = (hit_taxid_str, hit_level_str)
            valid_hits[count_type] = hits
        return valid_hits

    @staticmethod
    def get_valid_lineage(valid_hits, lineage_map, read_id, count_type):
        # If the read aligned to something, then it would be present in the
        # summary file for count type, and correspondingly in valid_hits[
        # count_type], even if the hits disagree so much that the
        # "valid_hits" entry is just ("-1", "-1"). If the read didn't align
        # to anything, we also represent that with ("-1", "-1"). This ("-1",
        # "-1") gets translated to NULL_LINEAGE.
        hit_taxid_str, hit_level_str = valid_hits[count_type].get(
            read_id, ("-1", "-1"))
        hit_lineage = lineage_map.get(hit_taxid_str,
                                      taxid_lineage.NULL_LINEAGE)
        return taxid_lineage.validate_taxid_lineage(hit_lineage, hit_taxid_str,
                                                    hit_level_str)
