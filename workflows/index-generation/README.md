# Index Generation

### Directory overview:

#### **index-generation.wdl**:
workflow code to create the following assets:
* NT and NR indexes from NCBI with redundant sequences removed - we use this to develop additional indexes (notebaly minimap and diamond indexes).
  Subsequently, we compile a "short DB" comprising solely the accessions that were identified from non-host alignment. This database is then used
  to blast (blastx for NR blastn for NT) assembled contigs during the mNGS postprocessing phase.

    * By using blastx we translates nucleotide hits from NT into protein sequences using the three different reading frames, we then search this sequence against the protein database (NR in this case). We also get out alignment statistics which we display in the mNGS sample report.
    * By using blastn we can find regions of similarity between nucleotide sequences (in this case our contigs and NT). We also get out alignment statistics which we display in the mNGS sample report.

* nt_loc.marisa and nr_loc.marisa:
  * index to quickly access an accession and it's sequence from either NT or NR.
  * This is used in consensus genome and phylotree workflows to pull reference accessions (when CG run is kicked off from the mNGS report)
  * Generating the JSON file for read alignment visualizations
  * In postprocess step of mNGS (where we generate taxon_counts) this file is used to download accessions based on hit summary from diamond and minimap to create the "short DB" mentioned above.
* nt_info.marisa and nr_info.marisa:
  * index to quickly go from accession ID to accession name, sequence length, and offset information from either NT and NR
  * mainly used to generate JSON files for coverage viz to be consumed by the web app.
* Minimap indexes - chunked minimap index for non host alignment for NT
* diamond indexes - chunked diamond index for non host alignment for NR
* accession2taxid.marisa - index to quickly go from accession ID to taxon ID - used for determining the optimal taxon assignment for each read from the alignment
   results
* taxid-lineages.marisa:
  * index used to go from tax ID to taxonomy IDs (taxid for species, genus, family)
  * used for determining the optimal taxon assignment for each read from the alignment
   results (calling hits), for example if a read aligns to multiple distinct references, we need to assess at which level in the taxonomic hierarchy the multiple alignments reach consensus.
  * We also use this file for generating taxon counts.
* deuterostome_taxids.txt - used to filter out eukaryotic sequences which helps narrow down taxon_counts to microbial DNA (bacteria, viruses, fungi, and parasites).
* taxon_ignore_list.txt - taxa that we would like to ignore (synthetic, constructs, plasmids, vectors, etc) in taxon_counts from non host alignment


#### **ncbi-compress**:
compression code written in rust to remove redundant sequences from NT and NR

#### **helpful_for_analysis**
jupyter notebooks used for helpful common analysis steps including:
* querying NCBI for an accession and lineage (used to investigate reads in the "all taxa with neither family nor genus classification" report for mNGS)
* querying marisa trie files - notebook to easily query all marisa trie files generated from index generation workflow
* compare non host alignment times between two projects - this was used to benchmark how long it took to do non host alignment on old and new indexes.
* generate taxon lineage changelog for a sample report - get a readout on which reads from a sample report have a taxon / lineage change between new and old index runs. Used for comp bio validation purposes mainly.
* checking sequence retention for different index compression runs - this notebook was handy for running multiple compression runs and summarizing which reads were dropped, helpful for early analysis and benchmarking of compression workflow.

#### **ncbitax2lin.py**
used to generate taxid-lineages.marisa

#### **generate_accession2taxid.py**
used to generate accession2taxis.marisa

#### **generate_ncbi_db_index.py**
used to generate nt_loc.marisa and nr_loc.marisa

#### **generate_lineage_csvs.py**
used to generate versioned-taxid-lineages.csv file used to populate taxon_lineage database table, also generates changelogs for deleted taxa, changed taxa, and new taxa.

### Debugging Notes:
* usually you need to launch an EC2 instance to test this workflow out at scale: `aegea launch --instance-type i3.16xlarge --no-provision-user --storage /=128GB --iam-role idseq-on-call <instance-name>`
* install rust: `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
* install python:
  ```
  curl -O https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
  chmod u+x Anaconda3-2023.09-0-Linux-x86_64.sh
  ./Anaconda3-2023.09-0-Linux-x86_64.sh
  ```

* if running miniwdl:
     make sure the version is 1.5.1 (`pip3 install miniwdl==1.5.1`)
    * downgrade importlib-metadata: `pip3 install importlib-metadata==4.13.0`

* if building the images on an EC2 machine:
    Follow [these instructions](https://docs.docker.com/engine/install/ubuntu/) to update docker before building

    * change docker build directory to be on /mnt:
        ```
        sudo systemctl stop docker
        sudo mv /var/lib/docker /mnt/docker
        sudo ln -s /mnt/docker /var/lib/docker
        sudo systemctl start docker
        ```

    * add the current user to the docker group:
        * `sudo usermod -aG docker $USER`
        * logout and login to reflect new group status

* build index-generation docker image: `make rebuild WORKFLOW=index-generation`
* run and exec into container: `docker run -v /mnt:/mnt --rm -it index-generation:latest bash`
* to run a certain task in index-geneation.wdl with miniwdl:
* `miniwdl run index-generation.wdl --task GenerateLocDB --input generate_loc_db_input.json`



