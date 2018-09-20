''' Generate Phylogenetic tree '''
import os
import json
import shelve
import traceback
import xml.etree.ElementTree as ET

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.steps.generate_alignment_viz import PipelineStepGenerateAlignmentViz
import idseq_dag.util.command as command
import idseq_dag.util.s3 as s3
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepGeneratePhyloTree(PipelineStep):
    ''' 
    Generate a phylogenetic tree from the input fasta files using kSNP3:
    http://gensoft.pasteur.fr/docs/kSNP3/01/kSNP3.01%20User%20Guide%20.pdf
    Augment the inputs with
      (a) NCBI sequences for the accession IDs specified in align_viz_json
    and/or
      (b) Genbank full reference genomes from the same taxid.
    '''
    def run(self):
        output_files = self.output_files_local()
        taxid = self.additional_attributes["taxid"]
        reference_taxids = self.additional_attributes.get("reference_taxids", [taxid]) # Note: will only produce a result if species-level or below
        superkingdom_name = self.additional_attributes.get("superkingdom_name")

        # Retrieve IDseq taxon fasta files
        local_taxon_fasta_files = []
        for pipeline_run_id, byterange_dict in self.additional_attributes["taxon_byteranges"].items():
            partial_fasta_files = []
            for hit_type, byterange in byterange_dict.items():
                first_byte = byterange[0]
                last_byte = byterange[1]
                s3_file = byterange[2]
                local_basename = f"{pipeline_run_id}_{hit_type}.fasta"
                bucket, key = s3.split_identifiers(s3_file)
                local_file = os.path.join(self.output_dir_local, local_basename)
                s3.fetch_byterange(first_byte, last_byte, bucket, key, local_file)
                partial_fasta_files.append(local_file)
            full_taxon_fasta = f"{self.output_dir_local}/{pipeline_run_id}.fasta"
            self.fasta_union(partial_fasta_files, full_taxon_fasta)
            local_taxon_fasta_files.append(full_taxon_fasta)
            for fasta in partial_fasta_files + [full_taxon_fasta]:
                print(f"{count.reads(fasta)} reads in {fasta}")

        # Trim Illumina adapters
        # TODO: consider moving this to the beginning of the main pipeline
        self.trim_adapters_in_place(local_taxon_fasta_files)

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
        genbank_fastas = self.get_genbank_genomes(reference_taxids, input_dir_for_ksnp3, superkingdom_name, 3)

        # Retrieve NCBI NT references for the accessions in the alignment viz files.
        # These are the accessions (not necessarily full genomes) that were actually matched
        # by the sample's reads during GSNAP alignment.
        accession_fastas = self.get_accession_sequences(input_dir_for_ksnp3, 10)

        # Retrieve NCBI metadata for the accessions
        metadata_by_node = self.get_metadata_by_tree_node({**accession_fastas, **genbank_fastas})
        metadata_output = output_files[1]
        with open(metadata_output, 'w') as f:
            json.dump(metadata_by_node, f)

        # Run MakeKSNP3infile.
        ksnp3_input_file = f"{self.output_dir_local}/inputs.txt"
        command.execute(f"cd {input_dir_for_ksnp3}/..; MakeKSNP3infile {os.path.basename(input_dir_for_ksnp3)} {ksnp3_input_file} A")

        # Specify which genomes should be used for annotation.
        # Specify the names of the genomes that should be used for annotation.
        # Here, we use the full genomes from genbank.
        annotated_genome_input = f"{self.output_dir_local}/annotated_genomes"
        genbanbk_fasta_files = list(genbank_fastas.values())
        if genbanbk_fasta_files:
            grep_options = " ".join([f"-e '{path}'" for path in genbanbk_fasta_files])
            command.execute(f"grep {grep_options} {ksnp3_input_file} | cut -f2 > {annotated_genome_input}")

        # Now run ksnp3.
        # We can choose among 4 different output files, see http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081760#s2:
        # (1) tree.parsimony.tre: basic, includes no node labels
        # (2) tree_AlleleCounts.parsimony.tre: labels the internal nodes with the number of SNPs that are shared exclusively by the descendants of that node
        # (3) tree_tipAlleleCounts.parsimony.tre: same as (2), but also labels the strain names at the tips with the number of SNPs that are exclusive to that strain.
        # (4) tree_AlleleCounts.parsimony.NodeLabel.tre: labels the internal nodes with the node number separated by an underscore from the number of SNPs that are
        #     shared exclusively by the descendants of that node.
        # Note: for integration with idseq-web, the node names need to be the pipeline_run_ids. So if we wanted to use outputs (2)/(3)/(4),
        # we would need to parse the appended information out from the newick node names and put it in a separate data structure.
        k_config = {
            # All entries to be revisited and benchmarked.
            # Values for viruses and bacteria come from kSNP3 recommendations (13-15 / 19-21).
            "Viruses": 14,
            "Bacteria": 20,
            "Eukaryota": 20,
            None: 14
        }
        k = k_config[superkingdom_name]
        ksnp_cmd = (f"cd {self.output_dir_local}; mkdir ksnp3_outputs; "
                    f"kSNP3 -in inputs.txt -outdir ksnp3_outputs -k {k}")
        if os.path.isfile(annotated_genome_input):
            ksnp_cmd += f" -annotate {os.path.basename(annotated_genome_input)}"
            # Note: produces SNP annotation file in a human-readable format that's very space inefficient (~100 MB).
            # May want to do some postprocessing once we know better what exactly we need for the web app.
        command.execute(ksnp_cmd)
        command.execute(f"mv {self.output_dir_local}/ksnp3_outputs/tree.parsimony.tre {output_files[0]}")
        snp_annotation_output = f"{self.output_dir_local}/ksnp3_outputs/SNPs_all_annotated"
        if os.path.isfile(snp_annotation_output):
            self.additional_files_to_upload.append(snp_annotation_output)
        else:
            log.write(f"Warning: {snp_annotation_output} was not generated!")


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

    def get_accession_sequences(self, dest_dir, n=10):
        '''
        Retrieve NCBI NT references for the most-matched accession in each alignment viz file, up to a maximum of n references.
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
        s3_align_viz_files = self.additional_attributes["align_viz_files"].values()
        local_align_viz_files = []
        for s3_file in s3_align_viz_files:
            local_basename = s3_file.replace("/", "-").replace(":", "-") # needs to be unique locally
            local_file = s3.fetch_from_s3(
                s3_file,
                os.path.join(self.ref_dir_local, local_basename))
            if local_file != None:
                local_align_viz_files.append(local_file)

        # Choose accessions to process.
        # align_viz files are a bit brittle, so we just log exceptions rather than failing the job.
        accessions = set()
        for local_file in local_align_viz_files:
            try:
                with open(local_file, 'rb') as f:
                    align_viz_dict = json.load(f)
                most_matched_accession = None
                max_num_reads = 0
                for acc, info in align_viz_dict.items():
                    num_reads = info["coverage_summary"]["num_reads"]
                    if num_reads > max_num_reads:
                        max_num_reads = num_reads
                        most_matched_accession = acc
                accessions.add(most_matched_accession)
                if len(accessions) >= n:
                    break
            except:
                log.write(f"Warning: couldn't get accession from {local_file}!")
                traceback.print_exc()
        if len(accessions) > n:
            accessions = set(list(accessions)[0:n])

        # Make map of accession to sequence file
        accession2info = dict((acc, {}) for acc in accessions)
        nt_loc_dict = shelve.open(nt_loc_db.replace(".db", ""))
        PipelineStepGenerateAlignmentViz.get_sequences_by_accession_list_from_s3(
            accession2info, nt_loc_dict, nt_db)

        # Put 1 fasta file per accession into the destination directory
        accession_fastas = {}
        for acc, info in accession2info.items():
            clean_accession = self.clean_name_for_ksnp3(acc)
            local_fasta = f"{dest_dir}/NCBI_NT_accession_{clean_accession}"
            command.execute(f"ln -s {info['seq_file']} {local_fasta}")
            command.execute(f"echo '>{clean_accession}' | cat - {local_fasta} > temp_file && mv temp_file {local_fasta}")
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

    @staticmethod
    def trim_adapters_in_place(local_input_files):
        for local_file in local_input_files:
            local_file_trimmed = os.path.join(os.path.dirname(local_file), "trimmed_" + os.path.basename(local_file))
            command.execute(f"cutadapt -a AGATCGGAAGAGCACACGTCT -o {local_file_trimmed} {local_file}")
            command.execute(f"mv {local_file_trimmed} {local_file}")

    @staticmethod
    def clean_name_for_ksnp3(name):
        return name.replace(' ', '-').replace('.', '')

    @staticmethod
    def fasta_union(partial_fasta_files, full_fasta_file):
        ''' Takes a list of fasta file paths and writes the union of the fasta records to full_fasta_file. '''
        if len(partial_fasta_files) == 1:
            command.execute(f"ln -s {partial_fasta_files[0]} {full_fasta_file}")
            return
        # For now, just load the files into memory. They are relatively small and
        # the same approach is used in the web app to serve taxon fasta files.
        # TODO: consider changing it to a solution that requires less memory.
        full_fasta = set()
        for fasta_file in partial_fasta_files:
            with open(fasta_file, 'r') as f:
                fasta = set(f.read().split(">"))
            full_fasta.update(fasta)
        full_fasta.discard('')
        with open(full_fasta_file, 'w') as f:
            f.write('>' + '>'.join(full_fasta))
