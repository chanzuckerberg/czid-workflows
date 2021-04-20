# IDseq PhyloTree-NG workflow

This is the IDseq PhyloTree-NG (next generation) workflow. Given a set of IDseq samples and a taxon ID of the organism
of interest, the workflow computes a split-kmer phylogeny with pre-filtering by hierarchical clustering, using
[SKA](https://github.com/simonrharris/SKA) split k-mers to create phylogenetic trees for an input set of samples. This
works for complete genomes as well as raw sequences (currently supporting only `.fasta` inputs).

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
Filename    | Provenance
------------|------------
To be added | To be added
