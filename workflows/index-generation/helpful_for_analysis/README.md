### jupyter notebooks used for helpful common analysis steps including:
* querying NCBI for an accession and lineage (used to investigate reads in the "all taxa with neither family nor genus classification" report for mNGS)
* querying marisa trie files - notebook to easily query all marisa trie files generated from index generation workflow above.
* compare non host alignment times between two projects - this was used to benchmark how long it took to do non host alignment on old and new indexes.
* generate taxon lineage changelog for a sample report - get a readout on which reads from a sample report have a taxon / lineage change between new and old index runs. Used for comp bio validation purposes mainly.
* checking sequence retention for different index compression runs - this notebook was handy for running multiple compression runs and summarizing which reads were dropped, helpful for early analysis and benchmarking of compression workflow.
