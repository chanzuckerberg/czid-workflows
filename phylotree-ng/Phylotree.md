# PhyloTree Workflow: Technical Description

The following documentation applies to the legacy kSNP3 PhyloTree workflow. The new SKA-based PhyloTree-NG workflow is a
work in progress. The documentation below will be adapted to the new workflow over time.

## Inputs

- **Samples**: An array with information about each sample/pipeline run used to make the phylo tree

  - **Sample Name**: The name the user has given the sample. Used to add user-friendly names to the VCF file that is parsed to construct the phylo tree for the user.

  - **Run ID**: The run ID of the pipeline run.

  - **Hit Summary NR**: The refined hit summary produced by rapsearch2 + our refinement. Used to compute the top N matched accessions.

  - **Hit Summary NT**: The refined hit summary produced by gsnap + our refinement. Used to compute the top N matched accessions.

  - **Taxon Byterange NR**: Contains an s3 path to the refined, sorted, annotated nr fasta as well as the byte range for the reads corresponding to the taxon ID.

  - **Taxon Byterange NT**: Contains an s3 path to the refined, sorted, annotated nt fasta as well as the byte range for the reads corresponding to the taxon ID.

- **Taxon ID**: Supplied by the user, the taxon ID of the organism of interest

- **Superkingdom Name**: The superkingdom that the organism of interest belongs to

- **Reference Taxon IDs**: A unique list of the species-level taxon IDs with the most hits from each pipeline run/sample. To compute these, for each pipeline run we take the select the species-level taxon with the highest count with the user-supplied taxon somewhere in it's lineage. Then the resulting list is de-duplicated. Note: if the original taxon ID is species-level or below the only reference taxon ID will be the original taxon ID.

- **NT Location DB**: A shelve key value store where the keys are accession IDs and the values are the byte ranges in the NT database for that accession
 
- **NT DB S3 Path**: An s3 path to the NT database. Note: this is an s3 path not a file because we only download parts of this based on byte ranges.

## Outputs

- **Phylo Tree Newick**: File Representing the phylo tree itself. This is saved as a string to the database and sent directly to the front end to build the phylo tree. The nodes on the tree have unique IDs corresponding to the name of the tree node. The IDs either correspond to pipeline run IDs for nodes that come from samples/pipeline runs, `genbank_{assembly_accession_id}.fasta` for nodes that come from genomes downloaded from genbank for reference taxon IDs, or `NCBI_NT_accession_{accession_id}.fasta` for nodes that come from reference sequences downloaded from the NCBI NT reference database.

- **NCBI Metadata JSON**: JSON containing metadata for sequences downloaded for reference taxon IDs from genbank the NCBI NT reference database. The keys are node IDs as described in the **Phylo Tree Newick** output and the values are the metadata from that ID. The metadata is: accession name, the accession ID, the country the accession was collected in, and the date the accession was collected on. This data is added to the user-facing phylo tree visualization.

- **SNPs Annotated**: (Optional) This output is only produced if there are any reference genomes. This is not used for anything by the web application but it is made available for the user to download if it is produced. This is the `SNPs_all_annotated` output from `kSNP3`. *Note*: in practice this is often absent or empty.

- **Variants Reference**: (Optional) This output is only produced if there are any reference genomes. This is not used for anything by the web application but it is made available for the user to download if it is produced. This is the output produced by `kSNP3` when the `-vcf` flag is provided. It is a "vcf file using the first genome specified in the -positions file as the reference genome".

## Algorithm

### Prepare Taxon Fastas

For each sample/pipeline run a fasta is constructed by downloading the reads from the **refined, sorted, annotated fasta** using the **taxon byteranges for both NR and NT**. The resulting NR and NT fastas are then concatenated. Finally, `cutadapt` is run to trim illumina adapters from each fasta. Each fasta is named based on it's pipeline run ID (ex. `123.fasta`). This naming is relevant to the later phylo tree construction as the file names are used directly by `MakeKSNP3infile`.

### Prepare kSNP3 Input

#### Place Taxon Fastas

Symlink the **taxon fastas** (output from the previous step) into a directory, which will be referred to as the kSNP3 inputs directory.

#### Download Genbank Genomes

Currently this step is skipped. See note at the end for explaination.

`n` genbank genomes per genbank category are then downloaded and added to the kSNP3 inputs directory in the following mannor.

The **superkingdom name** is used to get a list of category names. The mapping of superkingdom to category is: `Viruses` -> [`viral`], `Bacteria` -> [`bacteria`], `Eukaryota` -> [`fungi`, `protozoa`]. If there is no **superkingdom name** provided all possible category names are used, though this never happens in practice. For each of these category names the assembly summary for that category is downloaded from a copy kept in S3, though the original comes from `ftp://ftp.ncbi.nih.gov/genomes/genbank/{CATEGORY}/assembly_summary.txt`.

For each taxon ID in the **reference taxon IDs** the top `m` genomes are selected from the assembly summary where `m = max(n / len(reference taxon IDs), 1)`. These genomes are slected by using `cut` to parse the tsv, `awk` to filter for rows containing the taxon ID in question, and head to select the first `m` rows in the order they appear in the files.

```bash
cut -f1,6,7,8,20 $ASSEMBLY_SUMMARY | awk -F '\t' "$2 == $TAXON_ID" | head -n $M
```

The selected genomes for each taxon ID in the **reference taxon IDs** are combined into a single, deduplicated, list so in practice, fewer than `n` genomes may be downloaded. Fewer than `n` genomes may also be downloaded if there aren't enough hits for the taxon ID. There is a guard that ensures no more than `n` genomes from the resulting genome list are selected in case there are more than `n` reference taxon IDs. In that case each taxon ID would select at least one genome so there could be more than `n` provided the selections are not duplicates. If there are more than `n` selected genomes the ones later in the list are discarded. There is no particular order to the **reference taxon IDs** so the discard is non-deterministic with respect to the whole process of phylo tree generation, though the pipeline requires them in a specific order so the pipeline itself is deterministic. Once the genomes have been selected, the assembly summary file is deleted and each selected genome is downloaded directly from genbank, unzipped, and placed in the kSNP3 inputs directory. The ftp download URI is constructed by appending `/{basename(ftp_path)}_genomic.fna.gz` to the `ftp_path` from the assembly summary file. The downloaded files are named `genbank_{assembly_accession_id}.fasta` where `assembly_accession_id` corresponds to the field from the assembly summary file with spaces replaced by `-`s and `.`s removed. This file naming is relevant because these file names are used directly by `MakeKSNP3infile`. A map of that `assembly_accession_id` to the full path of the fasta is also saved from this step for use in generating the NCBI Metadata JSON.

*Note*: Currently `n` is hard coded to be `0` in which case this entire step is skipped. The default parameter for `n` the function that governs this step is `10`. There is no way for this step to be run by the caller and it has been disabled in this way for two years. It is unclear if this was ever used.


#### Download Accession Sequences

The NCBI NT references are downloaded for the most matched accession for each sample up to `n` (hardcoded to `10`) references. If there are more than `n` references to be downloaded, the top `n` by number of matches are selected. The selection and download are conducted in the following mannor.

For each sample, a count of the hits for each accession with the **taxon ID** anywhere in it's lineage are computed from the **refined hitsummaries from gsnap and rapsearch2**. The accession with the most hits for each sample has it's hits added to a counter. This means that if two samples have the same accession with the most hits their counts are combined. This also serves to deduplicate the accession sequences. This deduplication means that there may be fewer references downloaded than samples though there can never be more. If there are more than `n` accessions in the best accessions counter the `n` accessions with the highest combined number of hits are selected.

For each selected accession, the accession's byteranges are looked up in the **NT location database** and then these byteranges are used to download the reference from the **NT reference database** in S3. If no hit is found a warning is printed and the job continues, this is another way this step can result in fewer than `n` sequences per category. The references files are added to the kSNP3 input directory by symlinking them, then `echo` and `cat` are used to turn the file into a fasta file by prepending the line `>{accession_id}` to it and piping that output to a temporary file. This temporary file is then moved to where the symlink used to be. The names of the resulting files are `NCBI_NT_accession_{accession_id}.fasta` where `accession_id` is the accession ID with spaces replaced by `-`s and `.`s removed. This file naming is relevant because these file names are used directly by `MakeKSNP3infile`. A map of that `accession_id` to the full path of the fasta is also saved from this step for use in generating the NCBI Metadata JSON.

#### Run MakeKSNP3infile

`MakeKSNP3infile` is run on the kSNP3 input directory. It creates the input file required to run `kSNP3`. `MakeKSNP3infile` uses the basenames of the files in that directory as node names for the resulting phylo tree. This is why the file names are so important. Automatic mode is used so the `MakeKSNP3infile` is called with the following positional parameters: the kSNP3 input directory, the kSNP3 input file path, and `A`.

```bash
MakeKSNP3infile $KSNP3_INPUT_DIRECTORY $KSNP3_INPUT_FILE A
```

### Run kSNP3

The main stage of the phylo generation is running `kSNP3`. There are a few preparation steps before it can be run.

If the **reference taxon IDs** produced any reference genomes an annotation file is created (for `kSNP3`s `-annotate` input). This file is created by using `grep` to search for every line in the `kSNP3` input file containing one of the fasta outputs produced by the reference taxon IDs, then using `cut` to grab the second column from that output, which corresponds to the name for that genome. This output is then piped to a file to be used for the `-annotate` input to `kSNP3`.

```bash
grep -e path/to/name-2.fasta -e path/to/name-2.fasta $KSNP3_INPUT_FILE | cut -f2 > $ANNOTATE_INPUT_FILE
```

The kmer length parameter (`-k`) parameter is decided based on the **superkingdom name**. The kmer length for each superkingdom is: `Viruses` -> `13`, `Bacteria` -> `19`, `Eukaryota` -> `19`. If none is provided `13` is used though this never happens in practice.

An output directory, which will be referred to as the kSNP3 output directory, is then created. The final `kSNP3` parameters are assembled. The are:

- `-in`: the kSNP3 input file generated by `MakeKSNP3infile`
- `-outdir`: the kSNP3 output directory
- `-k`: the kmer length as computed above

If an annotation file was created earlier the following parameters are also added:

- `-annotate`: the annotation file prepared above
- `-vcf`: this is a boolean flag that will cause `kSNP3` to produce a VCF file with respect to first reference genome in the annotation file

The resulting tree output is checked to ensure it is not empty as kSNP3 sometimes produces an empty tree instead of crashing when encountering certain errors. The tree output is then renamed with `.newick`. The VCF output has it's column headers re-mapped from **pipeline run IDs** to **sample names** to help users manually inspecting this file. The VCF output is also renamed such that the output name is the same each time.

### Create NCBI Metadata JSON

The metadata JSON file is computed based maps of accession IDs to fasta files generated when fetching reference genomes for the **reference taxon IDs**. A JSON file is created with a map of phylo tree node IDs for every node that is from a reference genome to metadata for that node. The metadata is: accession name, the accession ID, the country the accession was collected in, and the date the accession was collected on.

To make this map the key value pairs of accession IDs and fasta file paths are iterated through and each pair is used to create an entry in the node ID to metadata map. The node ID is the basename of the fasta file without the file extension. The accession ID is already in the key value pair as well.

The rest of the metadata is fetched from NCBI. First the accession ID is searched for by calling `curl` on `https://eutils.ncbi.nlm.nih.gov/entrez/eutil/esearch.fcgi?db=nuccore&term=$ACCESSION_ID&usehistory=y`. The XML from this response is then parsed. The first `WebEnv` and `QueryKey` tag's text values are then used to construct another URL for `curl`: `https://eutils.ncbi.nlm.nih.gov/entrez/eutil/efetch.fcgi?db=nuccore&query_key=$QUERY_KEY&WebEnv=$WEB_ENV&rettype=gb&retmode=xml`. This XML is also parsed. The first `GBSeq_definition` is used for the accession name. Then the first `GBSeq_feature-table` tag is selected. This is used to select first `GBFeature_quals` tag that is a child of that tag's first child. All of the children of that tag are checked for the first instance of a `GBQualifier_name` tag. If that tag has text equal to `country` or `collection_date` the first instance of the tag `GBQualifier_value` is searched for from the same starting node, the text is used as the value for either country or collection date, and the key is no longer searched for in subsequent children. These values are used to construct a map of metadatum names to metadatum values. The overall map of node IDs to this metadata is then encoded as JSON.
