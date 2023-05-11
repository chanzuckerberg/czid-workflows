# CZ ID Long-Read mNGS workflow
​
The CZ ID mNGS Nanopore pipeline was developed in collaboration with the bioinformatics team from Oxford Nanopore Technologies (ONT). 
​
The documentation presented here reflects the current CZ ID long-read-mngs pipeline status and processes, with descriptions of the parameters used. 
​
Further documentation on how to use the CZ ID long-read-mngs workflow can be found on the CZ ID help center -- including a [pipeline overview](https://chanzuckerberg.zendesk.com/hc/en-us/articles/13756558532884-CZ-ID-Pipeline-Overview) outlining the pipeline steps and details of the [initial validation](https://chanzuckerberg.zendesk.com/hc/en-us/articles/13895641006100-mNGS-Nanopore-Initial-Pipeline-Validation).
​
### Changelog
​

**v0.7.2** -- April 25, 2023 -- Algorithmic updates to improve coverage of taxa identified in samples
​
 - Stitches together multiple hits from the same read or contig to the same accession into longer hits. These longer hits are not used for the coverage visualization.
​

**v0.7.0** -- April 25, 2023 -- Algorithmic updates to improve coverage of taxa identified in samples
​
 - Updates the Flye version used from `v2.9.0` to `v2.9.2`
 - Increases the stringency of alignment when mapping reads to their respective contigs (in the `RunReadsToContigs` step) by requiring that <20% of the read be clipped. This increases the number of "loose reads" but reduces spurious assignment of reads to contigs.
 - Modifies the minimap2 command for NCBI NT alignment (in the `RunNTAlignment` step) to use `map-ont` instead of previously-used `asm20`. This adjustment compliments the new stringency by improving the mapping of reads that don't belong to contigs.
​

**v0.6.0** -- April 4, 2023 -- Initial pipeline release
