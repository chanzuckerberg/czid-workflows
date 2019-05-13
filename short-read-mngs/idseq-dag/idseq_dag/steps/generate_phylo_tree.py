''' Generate Phylogenetic tree '''
import os
import glob
import json
import traceback
import xml.etree.ElementTree as ET

from collections import defaultdict

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.steps.generate_alignment_viz import PipelineStepGenerateAlignmentViz
import idseq_dag.util.command as command
import idseq_dag.util.s3 as s3
import idseq_dag.util.log as log
import idseq_dag.util.count as count
import idseq_dag.util.convert as convert

from idseq_dag.util.dict import IdSeqDictValue, open_file_db_by_extension

class PipelineStepGeneratePhyloTree(PipelineStep):
    '''
    Generate a phylogenetic tree from the input fasta files using kSNP3:
    http://gensoft.pasteur.fr/docs/kSNP3/01/kSNP3.01%20User%20Guide%20.pdf
    Augment the inputs with
      (a) NCBI sequences for the top accession IDs found in hitsummary2_files
    and/or
      (b) Genbank full reference genomes from the same taxid.
    '''
    def run(self):
        output_files = self.output_files_local()
        local_taxon_fasta_files = [f for input_item in self.input_files_local for f in input_item]
        taxid = self.additional_attributes["taxid"]
        reference_taxids = self.additional_attributes.get("reference_taxids", [taxid]) # Note: will only produce a result if species-level or below
        superkingdom_name = self.additional_attributes.get("superkingdom_name")

        # knsp3 has a command (MakeKSNP3infile) for making a ksnp3-compatible input file from a directory of fasta files.
        # Before we can use the command, we symlink all fasta files to a dedicated directory.
        # The command makes certain unreasonable assumptions:
        # - current directory is parent directory of the fasta file directory
        # - file names do not have dots except before extension (also no spaces)
        # - file names cannot be too long (for kSNP3 tree building).
        input_dir_for_ksnp3 = f"{self.output_dir_local}/inputs_for_ksnp3"
        command.execute(f"mkdir {input_dir_for_ksnp3}")
        for local_file in local_taxon_fasta_files:
            command.execute(f"ln -s {local_file} {input_dir_for_ksnp3}/{os.path.basename(local_file)}")

        # Retrieve Genbank references (full assembled genomes).
        genbank_fastas = self.get_genbank_genomes(reference_taxids, input_dir_for_ksnp3, superkingdom_name, 0)

        # Retrieve NCBI NT references for the accessions in the alignment viz files.
        # These are the accessions (not necessarily full genomes) that were actually matched
        # by the sample's reads during GSNAP alignment.
        accession_fastas = self.get_accession_sequences(input_dir_for_ksnp3, taxid, 10)

        # Retrieve NCBI metadata for the accessions
        metadata_by_node = self.get_metadata_by_tree_node({**accession_fastas, **genbank_fastas})
        metadata_output = output_files[1]
        with open(metadata_output, 'w') as f:
            json.dump(metadata_by_node, f)

        # Run MakeKSNP3infile.
        ksnp3_input_file = f"{self.output_dir_local}/inputs.txt"
        command.execute(f"cd {input_dir_for_ksnp3}/..; MakeKSNP3infile {os.path.basename(input_dir_for_ksnp3)} {ksnp3_input_file} A")

        # Specify the names of finished reference genomes.
        # Used for annotation & variant-calling.
        annotated_genome_input = f"{self.output_dir_local}/annotated_genomes"
        reference_fasta_files = list(genbank_fastas.values()) + list(accession_fastas.values())
        if reference_fasta_files:
            grep_options = " ".join([f"-e '{path}'" for path in reference_fasta_files])
            command.execute(f"grep {grep_options} {ksnp3_input_file} | cut -f2 > {annotated_genome_input}")

        # Now build ksnp3 command:
        k_config = {
            # All entries to be revisited and benchmarked.
            # Values for viruses and bacteria come from kSNP3 recommendations (13-15 / 19-21).
            "Viruses": 13,
            "Bacteria": 19,
            "Eukaryota": 19,
            None: 13
        }
        k = k_config[superkingdom_name]
        ksnp_output_dir = f"{self.output_dir_local}/ksnp3_outputs"
        ksnp_cmd = (f"mkdir {ksnp_output_dir}; cd {os.path.dirname(ksnp_output_dir)}; "
                    f"kSNP3 -in inputs.txt -outdir {os.path.basename(ksnp_output_dir)} -k {k}")

        # Annotate SNPs using reference genomes:
        # TODO: fix gi vs accession problem
        if os.path.isfile(annotated_genome_input):
            ksnp_cmd += f" -annotate {os.path.basename(annotated_genome_input)}"
            self.optional_files_to_upload.append(f"{ksnp_output_dir}/SNPs_all_annotated")

        # Produce VCF file with respect to first reference genome in annotated_genome_input:
        if os.path.isfile(annotated_genome_input):
            ksnp_cmd += " -vcf"

        # Run ksnp3 command:
        command.execute(ksnp_cmd)

        # Postprocess output names in preparation for upload:
        command.execute(f"mv {ksnp_output_dir}/tree.parsimony.tre {output_files[0]}")
        ksnp_vcf_file = glob.glob(f"{ksnp_output_dir}/*.vcf")
        if ksnp_vcf_file:
            target_vcf_file = f"{ksnp_output_dir}/variants_reference1.vcf"
            self.name_samples_vcf(ksnp_vcf_file[0], target_vcf_file)
            self.additional_files_to_upload.append(target_vcf_file)

        # Upload all kSNP3 output files for potential future reference
        supplementary_files = [f for f in glob.glob(f"{ksnp_output_dir}/*")
                               if os.path.isfile(f) and
                               f not in self.additional_files_to_upload + self.optional_files_to_upload]
        self.additional_files_to_upload.extend(supplementary_files)

    def count_reads(self):
        pass

    def get_genbank_genomes(self, reference_taxids, destination_dir, superkingdom_name, n=10):
        '''
        Retrieve up to n GenBank reference genomes under the reference_taxids.
        Assumes reference_taxids are species-level or below.
        Also assumes they are all in the same superkingdom, which is the only thing we need in our application.
        Saves the references under file names compatible with MakeKSNP3infile.
        TODO: Retrieve the genomes from S3 rather than ftp.ncbi.nih.gov (JIRA/IDSEQ-334).
        '''
        if n == 0 or not reference_taxids:
            return {}
        n_per_taxid = max(n // len(reference_taxids), 1)
        genbank_categories_by_superkingdom = {
            "Viruses": ["viral"],
            "Bacteria": ["bacteria"],
            "Eukaryota": ["fungi", "protozoa"],
            None: ["bacteria", "viral", "fungi", "protozoa"]
        }
        # additional options in genbank that we probably don't need right now:
        # ["archaea", "plant",
        # "vertebrate_mammalian", "vertebrate_other", "invertebrate",
        # "other", "metagenomes"]
        categories = genbank_categories_by_superkingdom[superkingdom_name]
        for cat in categories:
            genome_list_path_s3 = f"s3://idseq-database/genbank/{cat}/assembly_summary.txt" # source: ftp://ftp.ncbi.nih.gov/genomes/genbank/{cat}/assembly_summary.txt
            genome_list_local = s3.fetch_from_s3(genome_list_path_s3, destination_dir)
            genomes = []
            for taxid in reference_taxids:
                cmd = f"cut -f1,6,7,8,20 {genome_list_local}" # columns: 1 = assembly_accession; 6 = taxid; 7 = species_taxid, 8 = organism_name, 20 = ftp_path
                cmd += f" | awk -F '\t' '$2 == {taxid}'" # try to find taxid in the taxid column (2nd column of the piped input)
                cmd += f" | head -n {n_per_taxid}" # take only top n_per_taxid results
                taxid_genomes = list(filter(None, command.execute_with_output(cmd).split("\n")))
                genomes += [entry for entry in taxid_genomes if entry not in genomes]
            genomes = genomes[:n]
            command.execute_with_output(f"rm {genome_list_local}")
            if genomes:
                genbank_fastas = {}
                for line in genomes:
                    assembly_accession, taxid, species_taxid, organism_name, ftp_path = line.split("\t")
                    ftp_fasta_gz = f"{ftp_path}/{os.path.basename(ftp_path)}_genomic.fna.gz"
                    tree_node_name = f"genbank_{self.clean_name_for_ksnp3(assembly_accession)}"
                    local_fasta = f"{destination_dir}/{tree_node_name}.fasta"
                    if os.path.isfile(local_fasta):
                        local_fasta = f"{local_fasta.split('.')[0]}__I.fasta"
                    command.execute(f"wget -O {local_fasta}.gz {ftp_fasta_gz}")
                    command.execute(f"gunzip {local_fasta}.gz")
                    genbank_fastas[assembly_accession] = local_fasta
                return genbank_fastas
        return {}

    @staticmethod
    def parse_tree(current_dict, results, key = None):
        """
        Produce a dictionary like:
          { "accession 1": { "coverage_summary": ... },
            "accession 2": { "coverage_summary": ... },
            ...
          }
        from a dictionary like:
          { "family taxid 1": {
              "genus taxid 1": {
                "species taxid 1": {
                  "accession 1": { "coverage_summary": ... },
                }
              },
              "genus taxid 2": {
                "species taxid 2": {
                  "accession 2": { "coverage_summary": ... },
                }
              }
            }
          }
        """
        if "coverage_summary" in current_dict:
            results[key] = current_dict
        else:
            for key2, sub_dict in current_dict.items():
                PipelineStepGeneratePhyloTree.parse_tree(sub_dict, results, key2)

    def get_accession_sequences(self, dest_dir, taxid, n=10):
        '''
        Retrieve NCBI NT references for the most-matched accession in each hitsummary2 file, up to a maximum of n references.
        Write each reference to a separate fasta file.
        '''
        if n == 0:
            return {}

        # Retrieve files
        nt_db = self.additional_attributes["nt_db"]
        nt_loc_db = s3.fetch_from_s3(
            self.additional_files["nt_loc_db"],
            self.ref_dir_local,
            allow_s3mi=True)

        # Choose accessions to process.
        s3_hitsummary2_files = self.additional_attributes["hitsummary2_files"].values()
        accessions = defaultdict(lambda: 0)
        for file_list in s3_hitsummary2_files:
            tally = defaultdict(lambda: 0)
            for s3_file in file_list:
                local_basename = s3_file.replace("/", "-").replace(":", "-")
                local_file = s3.fetch_from_s3(
                    s3_file,
                    os.path.join(self.ref_dir_local, local_basename))
                if local_file is None:
                    continue
                with open(local_file, 'r') as f:
                    for line in f:
                        acc, species_taxid, genus_taxid, family_taxid = line.rstrip().split("\t")[3:7]
                        if any(int(hit_taxid) == taxid for hit_taxid in [species_taxid, genus_taxid, family_taxid]):
                            tally[acc] += 1
            if tally:
                best_acc, max_count = max(tally.items(), key=lambda x: x[1])
                accessions[best_acc] += max_count
        if len(accessions) > n:
            accessions = dict(sorted(accessions.items(), key=lambda x: x[1], reverse=True)[:n])
        accessions = set(accessions.keys())

        # Make map of accession to sequence file
        accession2info = dict((acc, {}) for acc in accessions)
        nt_loc_dict = open_file_db_by_extension(nt_loc_db, IdSeqDictValue.VALUE_TYPE_ARRAY)
        PipelineStepGenerateAlignmentViz.get_sequences_by_accession_list_from_s3(
            accession2info, nt_loc_dict, nt_db)

        # Put 1 fasta file per accession into the destination directory
        accession_fastas = {}
        for acc, info in accession2info.items():
            if info['seq_file'] is None:
                log.write(f"WARNING: No sequence retrieved for {acc}")
                continue
            clean_accession = self.clean_name_for_ksnp3(acc)
            local_fasta = f"{dest_dir}/NCBI_NT_accession_{clean_accession}.fasta"
            command.execute(f"ln -s {info['seq_file']} {local_fasta}")
            command.execute(f"echo '>{acc}' | cat - {local_fasta} > temp_file && mv temp_file {local_fasta}")
            accession_fastas[acc] = local_fasta

        # Return kept accessions and paths of their fasta files
        return accession_fastas

    @staticmethod
    def get_accession_metadata(accession):
        '''
        Retrieve metadata of an NCBI accession (e.g. name, country, collection date)
        TODO: Put this data in S3 instead and get it from there.
        '''
        accession_metadata = {}
        efetch_command = ";".join([
            f"QUERY={accession}",
            "BASE=https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
            "SEARCH_URL=${BASE}/esearch.fcgi?db=nuccore\&term=${QUERY}\&usehistory=y",
            "OUTPUT=$(curl $SEARCH_URL)",
            "WEB=$(echo $OUTPUT | sed -e 's/.*<WebEnv>\(.*\)<\/WebEnv>.*/\\1/')",
            "KEY=$(echo $OUTPUT | sed -e 's/.*<QueryKey>\(.*\)<\/QueryKey>.*/\\1/')",
            "FETCH_URL=${BASE}/efetch.fcgi?db=nuccore\&query_key=${KEY}\&WebEnv=${WEB}\&rettype=gb\&retmode=xml",
            f"curl $FETCH_URL"
        ])
        genbank_xml = command.execute_with_output(efetch_command)
        root = ET.fromstring(genbank_xml).find('GBSeq')
        if not root:
            log.write(f"WARNING: {efetch_command} did not give a result")
            return accession_metadata
        accession_metadata['name'] = root.find('GBSeq_definition').text
        qualifiers_needed = {'country', 'collection_date'}
        for entry in root.find('GBSeq_feature-table')[0].find('GBFeature_quals'):
            if all(key in accession_metadata for key in qualifiers_needed):
                break
            for key in qualifiers_needed - accession_metadata.keys():
                if entry.find('GBQualifier_name').text == key:
                    accession_metadata[key] = entry.find('GBQualifier_value').text
        return accession_metadata

    @staticmethod
    def get_metadata_by_tree_node(accession2fasta_map):
        metadata_by_node = {}
        for acc, fasta in accession2fasta_map.items():
            node = os.path.basename(os.path.splitext(fasta)[0]) # that's what kSNP3 chooses as the tree node name
            metadata = PipelineStepGeneratePhyloTree.get_accession_metadata(acc)
            metadata["accession"] = acc
            metadata_by_node[node] = metadata
        return metadata_by_node

    def name_samples_vcf(self, input_file, output_file):
        # The VCF has standard columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
        # followed by 1 column for each of the pipeline_run_ids. This function replaces the pipeline_run)_ids
        # by the corresponding sample names so that users can understand the file.
        sample_names_by_run_ids = self.additional_attributes["sample_names_by_run_ids"]
        vcf_columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        column_description_line = command.execute_with_output(f"awk /'^{vcf_columns}'/ {input_file}")
        column_description_line = column_description_line.strip()
        run_ids_in_order = [id for id in column_description_line.split("FORMAT\t")[1].split("\t") if convert.can_convert_to_int(id)]
        sample_names_in_order = [sample_names_by_run_ids.get(id, f"pipeline_run_{id}") for id in run_ids_in_order]
        new_column_description = '\t'.join([vcf_columns] + sample_names_in_order)
        command.execute(f"sed 's/^#CHROM.*/{new_column_description}/' {input_file} > {output_file}")

    @staticmethod
    def clean_name_for_ksnp3(name):
        return name.replace(' ', '-').replace('.', '')
