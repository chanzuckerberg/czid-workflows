# CZ ID host filter indexing

This folder contains the WDL workflow used to generate the precomputed index files for short-read-mngs host filtering, as well as a JSON file linking to the raw inputs for each supported host species.

## Adding a new host genome

First determine whether the host's RNA/cDNA transcripts are available for indexing (GTF and transcripts FASTA), in addition to the host's DNA genome (FASTA). The transcripts are needed for RNA quantification (abundance.tsv) and to maximize host filtering sensitivity. However, if the transcripts aren't available/convenient, then it's acceptable to index non-human host genomes without them.

To locate the genome & transcripts for a new host species, we usually refer to [Ensembl](https://uswest.ensembl.org/index.html) (including [Ensembl Metazoa](https://metazoa.ensembl.org/index.html), [Ensembl Fungi](http://fungi.ensembl.org/index.html), etc.), which curates and serves them in a uniform way. For example, consider Ensembl's [*C. elegans* landing page](https://uswest.ensembl.org/Caenorhabditis_elegans/Info/Index):

1. Genomic DNA: follow "Download DNA sequence (FASTA)" and find the file URL ending in `.dna_sm.toplevel.fa.gz`
2. Transcripts GTF: follow "Download GTF" and find the file URL ending in `.gtf.gz` (*not* with `abinitio`)
3. Transcripts FASTA: follow "Download FASTA" and find the `.fa.gz` files under `cdna/` and `ncrna/`. (Array of two URLs)

Add the necessary entry with these data file URLs to the JSON raw data file, then run the driver WDL on the new host.

If Ensembl doesn't cover the desired species, then find the genome FASTA in [NCBI Genomes](https://www.ncbi.nlm.nih.gov/genome/), and set `"transcripts_gtf_gz": null` and `"transcripts_fasta_gz": []`. For some hosts we use custom NCBI queries to get concatenated FASTA genomes for several species in a clade. We've also seen one case (white_shrimp) in which we could not use transcripts because they make the HISAT2 index build need an unreasonable amount of memory. (HISAT2's index building needs large memory when given transcripts, but the index file size and aligner memory usage are generally reasonable.) 

Lastly, upload the generated index files from `out/indexes/{species}` in the miniwdl run folder into `s3://czid-public-references/host_filtering_v2/{datestamp}/{species}/`.

## Example WDL invocations

TODO
