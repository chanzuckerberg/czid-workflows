# IDseq PhyloTree-NG workflow

This is the IDseq PhyloTree-NG (next generation) workflow. Given a set of IDseq samples and a taxon ID of the organism
of interest, the workflow computes a hierarchically clustered set of phylogenies of the samples (where each sample is
represented by contigs (assembled reads) from that sample mapped to reference sequences for the organism of interest).
Depending on coverage and phylogenetic diversity, the resulting phylogenies range in precision from a heatmap (pairwise
distance matrix) to a phylogram.

The workflow uses [SKA](https://github.com/simonrharris/SKA) split k-mers to create clusters and phylogenetic trees. This
works for complete genomes as well as raw sequences (currently supporting only `.fasta` inputs).

## Implementation details

### Current workflow
The workflow takes as input the main reference (represented by its taxon ID), additional references (represented by
taxon IDs or accession IDs), and samples (each represented by sample name, mngs workflow run ID, contig fasta file,
NT hit summary, and NR hit summary). The workflow then prepares FASTA files for each sample and reference:

- For samples, we scan the hit summaries to determine which contigs map to taxa under the given reference taxon ID,
  then subset the sample's contigs FASTA file.

- For references, we download their FASTA files from the AWS NCBI BLAST database S3 bucket using taxoniq.

We then feed the FASTA files to SKA to create kmer profiles and generate phylogenies.

### Old workflow
An in-depth description of the first generation PhyloTree pipeline from idseq-dag is available in
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
miniwdl run --verbose phylotree-ng/run.wdl docker_image_id=phylotree-ng data_directory=full_entero_data.tar.gz cut_height=.14 ska_align_p=.9
```

Note: the pipeline requires as input the `data_directory`, which is a directory containing .fasta files (one file per
sample) which is then tar zipped using the following command:

```
tar -czf full_entero_data.tar.gz full_entero_data
```

The optional parameters `cut_height` and `ska_align_p` allow you to specify the dendrogram cut height for pre-clustering
and the ska alignment proportion parameters, respectively. The default options are shown in the command above.

The `/analysis_support` directory contains scripts that have been used to support the experimentation and validation of
the phylotree pipeline.

## Reference files
Filename                | Provenance
------------------------|------------
`test/full_zika.tar.gz` | https://github.com/katrinakalantar/clustered-phylotree/blob/main/test/full_zika.tar.gz?raw=true
