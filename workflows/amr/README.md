# CZ ID AMR Workflow

CZ ID's AMR workflow implements the [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) tool for AMR sequence detection. The RGI tool is used to compare quality controlled reads and assembled contigs against AMR references sequences from the [Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/) (CARD).  Further documentation on how to use the CZ ID AMR workflow, including a [workflow visualization](https://chanzuckerberg.zendesk.com/hc/en-us/articles/15091031482644-AMR-Pipeline-Workflow) can be found in the [CZ ID help center](https://chanzuckerberg.zendesk.com/hc/en-us/categories/15001531592980-Antimicrobial-Resistance-Analysis).

# Changelog

## 1.2.10 - 2023-07-10

NOTE: Due to a bug in our release script, there were no releases with version numbers 1.2.6 - 1.2.9.

## 1.2.5 - 2023-06-13

### Fixed

- Fixed bug in calculating gene coverage. Contigs are now grouped by ARO accession instead of model ID before calculating coverage.
- Pinned `dask` version to `2023.5.0` to prevent installation of newer versions incompatible with Python 3.8.

### Changed
- Updated `requests` to version 2.31.0.

## 1.2.4 - 2023-05-23

### Added

- Added this README file.

### Changed

- Updated CARD database files.

## 1.2.3 - 2023-05-18

### Added

- Nudged results from running RGI on contigs are now included in the primary AMR report. These rows have the value "Nudged" in the "Cut_Off_contig_amr" column.

### Changed

- RGI is now called with the `--include_nudge` flag on contigs in task `RunRgiMain`.

## 1.2.2 - 2023-05-16

### Fixed

- Fixed a bug that caused creation of the indexed contig BAM file to error out if SPAdes assembly failed. The task now outputs an empty BAM file.

## 1.2.1 - 2023-05-10

Version 1.2.1 is unchanged from 1.2.0.

## 1.2.0 - 2023-05-10

### Fixed

- Fixed a bug that caused all gene names in the AMR reports to be entirely lowercase.

## 1.1.2 - 2023-05-10

### Added

- Added Python test for the WDL task which generates reports as TSV files, `RunResultsPerSample`.

### Fixed

- Fixed a bug that occurred when determining ARO accessions for rows in the primary AMR report.
- Removed carriage return characters from workflow test data.

## 1.1.1 - 2023-04-25

### Fixed

- Fixed the filename specified for the non-host reads file produced by `ZipOutputs`.

## 1.1.0 - 2023-04-25

### Changed

- The non-host reads file output by `ZipOutputs` is no longer compressed with gzip.

## 1.0.0 - 2023-04-25

### Added

- Contigs indexed by gene are now output as a BAM/BAI file pair from new task `tsvToSam`.
- A new column to the primary AMR report `read_gene_id` that can be used to look up sequences in the above-mentioned contigs bam file.

### Changed

- The non-host reads used for workflow runs on unfiltered samples now use the subsampled output from host filtering. Previously, the workflow used reads that were output from the host filtering process before subsampling.

### Fixed

- The workflow will now only interleave non-host reads as part of the `ZipOutputs` task if two files are present. Interleaving is skipped for workflows where all non-host reads are contained in a single file.

## 0.2.3 - 2023-04-13

### Added

- Dockerfile: Add SeqFu for interleaving input files.
- Non-host reads interleaved in a single file is now an output of `ZipOutputs`.
- Contigs are also now an output of `ZipOutputs`.

### Changed

- Renamed `ZipOutputs` input parameter `output_files` to `contigs_in`.

### Removed

- Removed input parameter `Int? total_reads`
- Removed task `GetTotalReads`.

## 0.2.2 - 2023-02-07

### Added

- New optional input parameter `Int? total_reads`.
- New task `GetTotalReads`.

### Fixed

- Load json with all CARD models when running `rgi load`.`

## 0.2.1 - 2022-12-15

### Added

- More detailed % coverage, now divided into contig coverage and read coverage. New columns have been added to the primary AMR report and some columns have been renamed in order to reflect this.
- New workflow input parameter `File wildcard_index`.

### Fixed

- Fixed typo in `MakeGeneCoverage` task.
- Fixed a bug in `MakeGeneCoverage` that did not account for zero-length output when writing results to `gene_coverage.tsv`.
- Output RunSpades.scaffolds is now optional type `File?` to reflect the fact that the `scaffolds.fasta` file may be missing in some SPAdes runs.`

### Changed

- AMR workflow files have been moved from `amr/` to `workflows/amr/`.`
- Dockerfile: `bedtools` installed via apt instead of source to eliminate the long compile time when building the image.
- Dockerfile: install RGI from commit 20b22dab on the rzlim08 fork.
- CARD databases are now loaded at runtime, and have been removed from the Dockerfile.
- For contigs, RGI is no longer called with the `--include_nudge` flag.
- Updated `RunRgiMain` Python test.

## 0.2.0 - 2022-10-03

### Added

- Sample name now supported as input and appears in workflow output.
- Synthesized report `primary_AMR_report.tsv`, a table containing the most important workflow output.
- Intermediate file with allele analysis output from k-mer taxonomic classification on non-host reads is now a workflow output.
- Mapped read stats from RGI main analysis of non-host reads available in several output files as workflow outputs.
- Mapped reads from RGI main analysis of non-host reads available as a BAM file with BAI index in workflow outputs.
- Intermediate outputs added to outputs.zip, organized into subfolders.
- Added test for RGI main analysis of failed assembly output.

### Fixed

- Workflow now continues properly with a sample that failed assembly (no contigs).
- Instead of erroring out, workflow now writes empty AMR report if there is no output from RGI.

### Changed

- Use main branch of rzlim08 RGI fork.
- Updated lxml dependency to version `4.9.1`.
- Rename final summary file to `comprehensive_AMR_metrics.tsv`.
- Primary AMR workflow output now in `primary_AMR_report.tsv`, from task `RunResultsPerSample`.
- Several intermediate output files have their names change to be consistent with snake_case.

## 0.1.2 - 2022-09-15

### Fixed

- Run SPAdes with a fixed number of threads (36) to prevent issues with using `$(nproc --all)`.

### Changed

- SPAdes stdout stream redirected to stderr for logging purposes.

## 0.1.1 - 2022-09-09

### Added

- Support for raw reads as input in paired fasta/fastq files.
- Host filtering via calling the short read mNGS host filter workflow.
- WDL task `RunSpades` to assemble contigs from host filter output using SPAdes genome assembler using short read mNGS Docker image.
- WDL task `MakeGeneCoverage` to calculate % coverage of reference genome for identified sequences.

### Changed

- Uses different rzlim08 RGI fork branch that adds gene hit coordinates.
- RGI main tasks now take in sequences from either top-level workflow inputs or host filtering/SPAdes output using `select_first()`.

## 0.1.0 - 2022-08-25

### Added

- Python requirements.txt file to run workflow without using conda.
- Dockerfile based on Ubuntu 20.04 image, with Python 3.8 installed.
- rzlim08 RGI fork that fixes header parsing.
- CARD RGI dependencies installed in Docker build.
- CARD databases added as part of Docker image build.
- WDL workflow support for samples with non-host reads and contigs.
- Workflow input parameter `Array[File] non_host_reads`.
- Workflow input parameter `File contigs`.
- Workflow input parameter `File card_json`.
- Workflow input parameter `File kmer_db`.
- Workflow input parameter `File amr_kmer_db`.
- Workflow input parameter `File wildcard_data`.
- WDL task `RunRgiMain` to identify input contigs using `rgi main`.
- WDL task `RunRgiBwtKma` to identify input non-host reads using `rgi bwt`.
- WDL task `RunRgiKmerMain` to run k-mer taxonomic classification on input contigs.
- WDL task `RunRgiKmerBwt` to run k-mer taxonomic classification on non-host reads.
- WDL task `RunResultsPerSample` that collects RGI output into TSV files; primary output in final_summary.tsv.
- WDL task `ZipOutputs` to collect primary analysis in one file location.
- Python tests for `RgiMain` and `RgiBwtKma` WDL tasks.
- Test contigs and non-host reads data.


# Reference Files

Filename | Provenance 
---------|-----------
s3://czid-public-references/card/2023-05-22/* | All files downloaded from https://card.mcmaster.ca/download on May 23, 2023. CARD version: 3.2.6 Wildcard version: 4.0.0
