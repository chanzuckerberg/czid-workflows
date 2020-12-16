''' Generate Phylogenetic tree '''
import os
import glob
import json
import itertools
import xml.etree.ElementTree as ET

from collections import defaultdict

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.steps.generate_alignment_viz import PipelineStepGenerateAlignmentViz
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.s3 as s3
import idseq_dag.util.log as log
import idseq_dag.util.convert as convert

from idseq_dag.util.dict import open_file_db_by_extension

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
        reference_taxids = self.additional_attributes.get("reference_taxids", [taxid])  # Note: will only produce a result if species-level or below
        # During phylo tree creation, if the taxon is in an unknown superkingdom then the k selected from k_config is supposed to be from the key None.
        superkingdom_name = self.additional_attributes.get("superkingdom_name") if self.additional_attributes.get("superkingdom_name") != '' else None

        # knsp3 has a command (MakeKSNP3infile) for making a ksnp3-compatible input file from a directory of fasta files.
        # Before we can use the command, we symlink all fasta files to a dedicated directory.
        # The command makes certain unreasonable assumptions:
        # - current directory is parent directory of the fasta file directory
        # - file names do not have dots except before extension (also no spaces)
        # - file names cannot be too long (for kSNP3 tree building).
        input_dir_for_ksnp3 = f"{self.output_dir_local}/inputs_for_ksnp3"
        command.make_dirs(input_dir_for_ksnp3)
        for local_file in local_taxon_fasta_files:
            command.execute(
                command_patterns.SingleCommand(
                    cmd="ln",
                    args=[
                        "-s",
                        local_file,
                        os.path.join(input_dir_for_ksnp3, os.path.basename(local_file))
                    ]
                )
            )

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
        command.execute(
            command_patterns.SingleCommand(
                cd=self.output_dir_local,
                cmd='MakeKSNP3infile',
                args=[
                    os.path.basename(input_dir_for_ksnp3),
                    ksnp3_input_file,
                    "A"
                ]
            )
        )

        # Specify the names of finished reference genomes.
        # Used for annotation & variant-calling.
        annotated_genome_input = f"{self.output_dir_local}/annotated_genomes"
        reference_fasta_files = list(genbank_fastas.values()) + list(accession_fastas.values())
        if reference_fasta_files:
            grep_options = (("-e", path) for path in reference_fasta_files)
            grep_options = list(itertools.chain.from_iterable(grep_options))  # flatmap
            command.execute(
                command_patterns.ShellScriptCommand(
                    script=r'''grep "${grep_options[@]}" "${ksnp3_input_file}" | cut -f2 > "${annotated_genome_input}";''',
                    named_args={
                        'ksnp3_input_file': ksnp3_input_file,
                        'annotated_genome_input': annotated_genome_input,
                        'grep_options': grep_options
                    }
                )
            )

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
        command.make_dirs(ksnp_output_dir)
        ksnp_cd = os.path.dirname(ksnp_output_dir)
        ksnp_cmd = "kSNP3"
        ksnap_args = [
            "-in",
            "inputs.txt",
            "-outdir",
            os.path.basename(ksnp_output_dir),
            "-k",
            k
        ]

        # Annotate SNPs using reference genomes:
        # TODO: fix gi vs accession problem
        if os.path.isfile(annotated_genome_input):
            ksnap_args.extend([
                "-annotate",
                os.path.basename(annotated_genome_input)
            ])
            snps_all_annotated = f"{ksnp_output_dir}/SNPs_all_annotated"
            if os.path.isfile(snps_all_annotated):
                self.additional_output_files_hidden.append(snps_all_annotated)

        # Produce VCF file with respect to first reference genome in annotated_genome_input:
        if os.path.isfile(annotated_genome_input):
            ksnap_args.append("-vcf")

        # Run ksnp3 command:
        command.execute(
            command_patterns.SingleCommand(
                cd=ksnp_cd,
                cmd=ksnp_cmd,
                args=ksnap_args
            )
        )
        tree_output = os.path.join(ksnp_output_dir, "tree.parsimony.tre")

        # kSNP3 sometimes outputs an empty file when it fails instead of
        #   returning a non-zero code. This effectively causes the
        #   step to silently fail which causes problems downstream
        #   If this happens we should fail here, and not upload
        #   the output
        assert os.path.getsize(tree_output) > 0, "kSNP3 encountered an issue and produced a tree file of size 0"

        # Postprocess output names in preparation for upload:
        command.move_file(tree_output, output_files[0])
        ksnp_vcf_file = glob.glob(f"{ksnp_output_dir}/*.vcf")
        if ksnp_vcf_file:
            target_vcf_file = f"{ksnp_output_dir}/variants_reference1.vcf"
            self.name_samples_vcf(ksnp_vcf_file[0], target_vcf_file)
            self.additional_output_files_hidden.append(target_vcf_file)

        # Upload all kSNP3 output files for potential future reference
        supplementary_files = [f for f in glob.glob(f"{ksnp_output_dir}/*")
                               if os.path.isfile(f) and
                               f not in self.additional_output_files_hidden]
        self.additional_output_files_hidden.extend(supplementary_files)

    def count_reads(self):
        pass

    @staticmethod
    def get_taxid_genomes(genome_list_local, taxid, n_per_taxid):
        cmd = command_patterns.ShellScriptCommand(
            script=(
                # columns: 1 = assembly_accession; 6 = taxid; 7 = species_taxid, 8 = organism_name, 20 = ftp_path
                r'''cut -f1,6,7,8,20 "${genome_list_local}" '''
                # try to find taxid in the taxid column (2nd column of the piped input)
                r''' | awk -F '\t' "${awk_find_pattern}" '''
                # take only top n_per_taxid results
                r''' | head -n "${n_per_taxid}";'''
            ),
            named_args={
                'genome_list_local': genome_list_local,
                'awk_find_pattern': f'$2 == "{taxid}"',
                'n_per_taxid': n_per_taxid
            }
        )
        taxid_genomes = list(filter(None, command.execute_with_output(cmd).split("\n")))
        return taxid_genomes

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
            genome_list_path_s3 = f"s3://idseq-public-references/genbank/{cat}/assembly_summary.txt"  # source: ftp://ftp.ncbi.nih.gov/genomes/genbank/{cat}/assembly_summary.txt
            genome_list_local = s3.fetch_from_s3(genome_list_path_s3, destination_dir)
            genomes = []
            for taxid in reference_taxids:
                taxid_genomes = PipelineStepGeneratePhyloTree.get_taxid_genomes(genome_list_local, taxid, n_per_taxid)
                genomes += [entry for entry in taxid_genomes if entry not in genomes]
            genomes = genomes[:n]
            command.remove_file(genome_list_local)
            if genomes:
                genbank_fastas = {}
                for line in genomes:
                    assembly_accession, taxid, _species_taxid, _organism_name, ftp_path = line.split("\t")
                    ftp_fasta_gz = f"{ftp_path}/{os.path.basename(ftp_path)}_genomic.fna.gz"
                    tree_node_name = f"genbank_{self.clean_name_for_ksnp3(assembly_accession)}"
                    local_fasta = f"{destination_dir}/{tree_node_name}.fasta"
                    if os.path.isfile(local_fasta):
                        local_fasta = f"{local_fasta.split('.')[0]}__I.fasta"
                    command.execute(
                        command_patterns.SingleCommand(
                            cmd='wget',
                            args=[
                                "-O",
                                f"{local_fasta}.gz",
                                ftp_fasta_gz
                            ]
                        )
                    )
                    command.execute(
                        command_patterns.SingleCommand(
                            cmd='gunzip',
                            args=[
                                f"{local_fasta}.gz"
                            ]
                        )
                    )
                    genbank_fastas[assembly_accession] = local_fasta
                return genbank_fastas
        return {}

    @staticmethod
    def parse_tree(current_dict, results, key=None):
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
        nt_loc_db = s3.fetch_reference(
            self.additional_files["nt_loc_db"],
            self.ref_dir_local,
            allow_s3mi=True)

        # Choose accessions to process.
        s3_hitsummary2_files = self.additional_attributes["hitsummary2_files"].values()
        accessions = defaultdict(lambda: 0)
        # TODO: Address issue where accessions in nr can be chosen in the following code.
        # These accessions will not be found in nt_loc and will be subsequently omitted.
        for file_list in s3_hitsummary2_files:
            tally = defaultdict(lambda: 0)
            for s3_file in file_list:
                local_basename = s3_file.replace("/", "-").replace(":", "-")
                local_file = s3.fetch_from_s3(
                    s3_file,
                    os.path.join(self.output_dir_local, local_basename))
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
        with open_file_db_by_extension(nt_loc_db) as nt_loc_dict:
            PipelineStepGenerateAlignmentViz.get_sequences_by_accession_list_from_s3(
                accession2info, nt_loc_dict, nt_db)

        # Put 1 fasta file per accession into the destination directory
        accession_fastas = {}
        for acc, info in accession2info.items():
            if 'seq_file' not in info or info['seq_file'] is None:
                log.write(f"WARNING: No sequence retrieved for {acc}")
                continue
            clean_accession = self.clean_name_for_ksnp3(acc)
            local_fasta = f"{dest_dir}/NCBI_NT_accession_{clean_accession}.fasta"
            command.execute(
                command_patterns.SingleCommand(
                    cmd="ln",
                    args=[
                        "-s",
                        info['seq_file'],
                        local_fasta
                    ]
                )
            )
            command.execute_with_output(
                command_patterns.ShellScriptCommand(
                    script=r'''echo ">${acc}" | cat - "${local_fasta}" > temp_file;''',
                    named_args={
                        'acc': acc,
                        'local_fasta': local_fasta
                    }
                )
            )
            command.move_file('temp_file', local_fasta)

            accession_fastas[acc] = local_fasta

        # Return kept accessions and paths of their fasta files
        return accession_fastas

    @staticmethod
    def fetch_ncbi(accession):
        query = accession
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        search_url = f"{base}/esearch.fcgi?db=nuccore&term={query}&usehistory=y"
        output = command.execute_with_output(
            command_patterns.SingleCommand(
                cmd="curl",
                args=[search_url]
            )
        )
        root = ET.fromstring(output)
        web = root.find('WebEnv').text
        key = root.find('QueryKey').text
        fetch_url = f"{base}/efetch.fcgi?db=nuccore&query_key={key}&WebEnv={web}&rettype=gb&retmode=xml"
        genbank_xml = command.execute_with_output(
            command_patterns.SingleCommand(
                cmd="curl",
                args=[fetch_url]
            )
        )
        return {
            'search_url': search_url,
            'fetch_url': fetch_url,
            'genbank_xml': genbank_xml
        }

    @staticmethod
    def get_accession_metadata(accession):
        '''
        Retrieve metadata of an NCBI accession (e.g. name, country, collection date)
        TODO: Put this data in S3 instead and get it from there.
        '''
        accession_metadata = {}
        fetch_ncbi_result = PipelineStepGeneratePhyloTree.fetch_ncbi(accession)
        genbank_xml = fetch_ncbi_result['genbank_xml']

        root = ET.fromstring(genbank_xml).find('GBSeq')
        if not root:
            log.write(f"WARNING: {fetch_ncbi_result} did not return a result")
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
            node = os.path.basename(os.path.splitext(fasta)[0])  # that's what kSNP3 chooses as the tree node name
            metadata = PipelineStepGeneratePhyloTree.get_accession_metadata(acc)
            metadata["accession"] = acc
            metadata_by_node[node] = metadata
        return metadata_by_node

    @staticmethod
    def _vcf_new_column_description(input_file, sample_names_by_run_ids):
        vcf_columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        column_description_line = command.execute_with_output(
            command_patterns.SingleCommand(
                cmd='awk',
                args=[
                    f"/^{vcf_columns}/",
                    input_file
                ]
            )
        )
        column_description_line = column_description_line.strip()
        additional_columns = column_description_line.split("FORMAT\t")[1].split("\t")
        sample_names_in_order = [
            sample_names_by_run_ids.get(id, f"pipeline_run_{id}") if convert.can_convert_to_int(id) else id
            for id in additional_columns
        ]
        new_column_description = '\t'.join([vcf_columns] + sample_names_in_order)
        return new_column_description

    @staticmethod
    def _vcf_replace_column_description(input_file, output_file, new_column_description):
        escaped_new_column_description = new_column_description.replace("\\", "\\\\").replace("&", "\\&").replace("/", r"\/")
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''sed "${sed_pattern}" "${input_file}" > "${output_file}"''',
                named_args={
                    'sed_pattern': f"s/^#CHROM.*/{escaped_new_column_description}/",
                    'input_file': input_file,
                    'output_file': output_file
                }
            )
        )

    def name_samples_vcf(self, input_file, output_file):
        # The VCF has standard columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
        # followed by 1 column for each of the pipeline_run_ids. This function replaces the pipeline_run)_ids
        # by the corresponding sample names so that users can understand the file.
        sample_names_by_run_ids = self.additional_attributes["sample_names_by_run_ids"]
        new_column_description = self._vcf_new_column_description(input_file, sample_names_by_run_ids)
        self._vcf_replace_column_description(input_file, output_file, new_column_description)

    @staticmethod
    def clean_name_for_ksnp3(name):
        return name.replace(' ', '-').replace('.', '')
