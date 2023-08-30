# CZ ID AMR Workflow

CZ ID's AMR pipeline implements the [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) tool for AMR sequence detection. The RGI tool is used to compare quality controlled reads and assembled contigs against AMR references sequences from the [Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/) (CARD).  Further documentation on how to use the CZ ID AMR workflow, including a [pipeline workflow](https://chanzuckerberg.zendesk.com/hc/en-us/articles/15091031482644-AMR-Pipeline-Workflow) can be found in the [CZ ID help center](https://chanzuckerberg.zendesk.com/hc/en-us/categories/15001531592980-Antimicrobial-Resistance-Analysis).

# Changelog

## [1.2.3] -- 2023-05-24 -- Initial pipeline release

## [0.1.0] -- 2022-08-25

### Added

- Python requirements.txt file to run workflow without conda
- Dockerfile based on Ubuntu 20.04 image
- CARD RGI dependencies installed in Docker build
- CARD databases added as part of Docker image build
- WDL workflow support for samples with non-host reads and contigs
- Add task RunRgiMain to analyze input contigs using `rgi main`
- Add task RunRgiBwtKma to analyze input non-host reads using `rgi bwt`
- Add task RunRgiKmerMain to run k-mer taxonomic classification on input contigs
- Add task RunRgiKmerBwt torun k-mer taxonomic classification on non-host reads
- Add task RunResultsPerSample that collects RGI output into TSV files; primary output in final_summary.tsv
- Add Python tests for RgiMain and RgiBwtKma WDL tasks
- Add test contigs and non-host reads data

### Fixed


### Changed


### Removed


## Reference Files

Filename | Provenance 
---------|-----------
s3://czid-public-references/card/2023-05-22/* | All files downloaded from https://card.mcmaster.ca/download on May 23, 2023. CARD version: 3.2.6 Wildcard version: 4.0.0
