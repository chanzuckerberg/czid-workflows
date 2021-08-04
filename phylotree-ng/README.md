# IDseq PhyloTree-NG workflow

This is the IDseq PhyloTree-NG (next generation) workflow. Given a set of 
IDseq samples and a taxon ID of the organism of interest, the workflow 
computes the relative genomic distance between samples and generates phylogenies of the samples. Each sample is represented by contigs (assembled 
reads) from the original input sample mapped to reference sequences for the 
organism of interest.

Depending on coverage and phylogenetic diversity, the resulting phylogenies 
range in precision from a heatmap (pairwise distance matrix) to a phylogram.

The workflow uses [SKA](https://github.com/simonrharris/SKA) split k-mers to 
create clusters and phylogenetic trees. This works for complete genomes as 
well as raw sequences (currently supporting only `.fasta` inputs).

## Implementation details

### Current workflow overview
The workflow takes as input the samples (each represented by sample name, mngs workflow run ID, contig fasta file, contig summary file) and reference accessions (idetified via a heuristic). The workflow then prepares FASTA files for each sample and reference:

- For samples, we scan their
  [merged hit summaries](https://github.com/chanzuckerberg/idseq-workflows/blob/main/short-read-mngs/idseq-dag/idseq_dag/util/parsing.py)
  to determine which contigs map to taxa under the given reference taxon ID, then subset the sample's contigs FASTA
  file.

- For references, we download their FASTA files from the AWS NCBI BLAST database S3 bucket using taxoniq.

We then feed the FASTA files to SKA to create kmer profiles and generate phylogenies.

### Detailed workflow description
A detailed outline of the workflow steps is available in [Phylotree-ng.md](Phylotree-ng.md).

### Old workflow
An in-depth description of the first generation PhyloTree pipeline from idseq-dag (in IDseq until August 2021) is available in
[Phylotree.md](Phylotree.md). The authoritative copies of the associated idseq-dag steps are in
[prepare_taxon_fasta.py](../short-read-mngs/idseq-dag/idseq_dag/steps/prepare_taxon_fasta.py) and
[enerate_phylo_tree.py](../short-read-mngs/idseq-dag/idseq_dag/steps/enerate_phylo_tree.py).

# Running the workflow locally

First, clone the repo

```
git clone https://github.com/chanzuckerberg/idseq-workflows.git
```

Then, build the docker image:

```
docker build -t phylotree-ng phylotree-ng
```

We then use the docker image to run the pipeline as follows:

```
miniwdl run --verbose phylotree-ng/run.wdl --input phylotree-ng/test/test_inputs.json
```

The inputs.json file will have the following structure:

```
{"samples":[

], "reference_taxon_id":64320,
"superkingdom_name":"viruses",
"additional_reference_accession_ids":["NC_012532.1", "NC_035889.1"],
"docker_image_id":"phylo-ng"}

# where the "samples" array contains instances of the following format
{"workflow_run_id":1, "sample_name":"X", "combined_contig_summary":".json", "contig_fasta":".fasta"}
```



The `/analysis_support` directory contains scripts that have been used to support the experimentation and validation of
the phylotree pipeline.

## Reference files
Filename                | Provenance
------------------------|------------
`test/full_zika.tar.gz` | https://github.com/katrinakalantar/clustered-phylotree/blob/main/test/full_zika.tar.gz?raw=true
