# [IDseq](https://idseq.net/) &middot; [![GitHub license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://github.com/chanzuckerberg/idseq-web/blob/master/LICENSE) ![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)

![logo](https://assets.idseq.net/assets/Logo_Black.png)

#### Infectious Disease Sequencing Platform
IDseq is a hypothesis-free global software platform that helps scientists identify pathogens in metagenomic sequencing data.

- **Discover** - Identify the pathogen landscape
- **Detect** - Monitor and review potential outbreaks
- **Decipher** - Find potential infecting organisms in large datasets

A collaborative open project of [Chan Zuckerberg Initiative](https://www.chanzuckerberg.com/) and [Chan Zuckerberg Biohub](https://czbiohub.org).

Check out our repositories:
- [idseq-web](https://github.com/chanzuckerberg/idseq-web) - Frontend portal
- [idseq-dag](https://github.com/chanzuckerberg/idseq-dag) - Bioinformatics pipeline and workflow engine (here)
- [idseq-cli](https://github.com/chanzuckerberg/idseq-cli) - Command line upload interface
- [idseq-bench](https://github.com/chanzuckerberg/idseq-bench) - Pipeline benchmarking tools

## IDSEQ-DAG

idseq_dag is the pipeline execution engine for idseq (see idseq.net). It is a pipelining system that implements a directed acyclic graph (DAG) where the nodes (steps) correspond to individual python classes. The graph is defined using JSON.

The pipeline would be executed locally with local machine resources. idseq-dag could be installed inside a docker container and run inside the container. See the [Dockerfile](Dockerfile) for our setup.  More details could be found below.

## Requirements

Python 3.6 and up

## Installation

```
cd idseq-dag; pip3 install -e .
```

## Example

```
idseq_dag examples/host_filter_dag.json

```

` idseq_dag --help ` for more options

## Test

```
cd idseq-dag; python3 -m unittest

```

or `python3 -m unittest tests/<module_file> ` for testing individual modules.

## DAG Execution Details

### Composing an example dag json file
There are five basic elements of an IdSeq Dag
 - **name**: the name of the IdSeq Dag.
 - **output_dir_s3**: the s3 directory where the output files will be copied to.
 - **targets**: the outputs that are to be generated through dag execution. Each target consists of a list of files that will be copied to *output_dir_s3*
 - **steps**: the steps that will be executed in order to generate the targets. For each step, the following attributes can be specified:
   * *in*: the input targets
   * *out*: the output target
   * *module*: name of python module
   * *class*: name of python class that inherits from [PipelineStep](idseq_dag/engine/pipeline_step.py)
   * *additional_files*: additional S3 files required for dag execution, i.e. reference files.
   * *additional_attributes*: additional input parameters for the pipeline class
 - **given_targets**: the list of targets that are given. Given targets will be downloaded from S3 before the pipeline execution starts.

The following is an example dag for generating alignment output for idseq. The *host_filter_out* is given, and once downloaded, *gsnap_out* and *rapsearch2_out* steps will run in parallel. When *gsnap_out* and *rapsearch2_out* are both completed, *taxon_count_out* and *annotated_out* will be run simultaneously and the pipeline will be complete once everything is uploaded to S3.

```
{
  "name": "alignment",
  "output_dir_s3": "s3://idseq-samples-prod/samples/12/5815/results_test",
  "targets": {
    "host_filter_out": [
        "gsnap_filter_1.fa"
          , "gsnap_filter_2.fa"
          , "gsnap_filter_merged.fa"
    ],
    "gsnap_out": [
      "gsnap.m8",
      "gsnap.deduped.m8",
      "gsnap.hitsummary.tab",
      "gsnap_counts.json"
    ],
    "rapsearch2_out": [
      "rapsearch2.m8",
      "rapsearch2.deduped.m8",
      "rapsearch2.hitsummary.tab",
      "rapsearch2_counts.json"
    ],
    "taxon_count_out": ["taxon_counts.json"],
    "annotated_out": ["annotated_merged.fa", "unidentified.fa"]
  },
  "steps": [
    {
      "in": ["host_filter_out"],
      "out": "gsnap_out",
      "class": "PipelineStepRunAlignment",
      "module": "idseq_dag.steps.run_alignment_remotely",
      "additional_files": {
        "lineage_db": "s3://idseq-public-references/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/taxid-lineages.db",
        "accession2taxid_db": "s3://idseq-public-references/alignment_data/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/accession2taxid.db"

          ,"deuterostome_db": "s3://idseq-public-references/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/deuterostome_taxids.txt"

      },
      "additional_attributes": {
        "alignment_algorithm": "gsnap",
      }
    },
    {
      "in": ["host_filter_out"],
      "out": "rapsearch2_out",
      "class": "PipelineStepRunAlignment",
      "module": "idseq_dag.steps.run_alignment_remotely",
      "additional_files": {
        "lineage_db": "s3://idseq-public-references/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/taxid-lineages.db",
        "accession2taxid_db": "s3://idseq-public-references/alignment_data/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/accession2taxid.db"

          ,"deuterostome_db": "s3://idseq-public-references/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/deuterostome_taxids.txt"

      },
      "additional_attributes": {
        "alignment_algorithm": "rapsearch2",
      }
    },
    {
      "in": ["gsnap_out", "rapsearch2_out"],
      "out": "taxon_count_out",
      "class": "PipelineStepCombineTaxonCounts",
      "module": "idseq_dag.steps.combine_taxon_counts",
      "additional_files": {},
      "additional_attributes": {}
    },
    {
      "in": ["host_filter_out", "gsnap_out", "rapsearch2_out"],
      "out": "annotated_out",
      "class": "PipelineStepGenerateAnnotatedFasta",
      "module": "idseq_dag.steps.generate_annotated_fasta",
      "additional_files": {},
      "additional_attributes": {}
    }
  ],
  "given_targets": {
    "host_filter_out": {
      "s3_dir": "s3://idseq-samples-prod/samples/12/5815/results_test/2.4"
    }
  }
}

```

### Design

idseq-dag follows the KISS design principle. There are only two major components for dag execution:

 -  *[PipelineFlow](idseq_dag/engine/pipeline_flow.py)* validates the DAG, downloads the given targets, starts prefetching the additional files in a different thread, starts the pipeline steps in parallel and coordinates the execution.
 -  *[PipelineStep](idseq_dag/engine/pipeline_step.py)* waits for the input targets to be available, executes the run method, validates the output and uploads the files to S3.

Here is a quick example of a PipelineStep implementation for generating [taxon_count_out](master/idseq_dag/steps/combine_taxon_counts.py). Anyone can implement their own step by subclassing PipelineStep and implementing the *run* and *count_reads* method. *count_reads* is a method we use to count the output files. If you are not sure how to implement the count_reads function, just put a dummy function there as shown below.

In the *run* function, you need to make sure your implementation generates all the files specified in the `self.output_files_local()` list. Otherwise, the step will fail, which will trigger the whole pipeline to also fail.

By default, idseq-dag will only execute a step when the output target is not generated yet. You can turn off this caching mechanism with `--no-lazy-run` option with `idseq_dag` command.

```

import json
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count

class PipelineStepCombineTaxonCounts(PipelineStep):
    '''
    Combine counts from gsnap and rapsearch
    '''
    def run(self):
        input_files = []
        for target in self.input_files_local:
            input_files.append(target[3])
        output_file = self.output_files_local()[0]
        self.combine_counts(input_files, output_file)


    def count_reads(self):
        pass
    @staticmethod
    def combine_counts(input_json_files, output_json_path):
        taxon_counts = []
        for input_file_name in input_json_files:
            with open(input_file_name, 'r') as input_file:
                data = json.load(input_file)
                taxon_counts += data["pipeline_output"]["taxon_counts_attributes"]
        output_dict = {"pipeline_output": { "taxon_counts_attributes": taxon_counts}}
        with open(output_json_path, 'w') as output_file:
            json.dump(output_dict, output_file)

```

## Developers
When merging a commit to master, you need to increase the version number in `idseq_dag/__init__.py`:
  - if results are expected to change, increase the 2nd number
  - if results are not expected to change, increase the 3rd number.

## Generating reference indices
IDseq DAGs require the use of several indices prepared from files in NCBI. If you need to generate a new version of these indices, please refer to https://github.com/chanzuckerberg/idseq-pipeline, specifically the following commands:
  - host_indexing
  - gsnap_indexing
  - rapsearch_indexing
  - blacklist
  - lineages
  - curate_accession2taxid
  - curate_accessionid2seq
  - push_reference_update.
TODO: Move this code over to the idseq-dag repo.

## Release notes
Version numbers for this repo take the form X.Y.Z.
- We increase Z for a change that does not add or change results in any way. Example: adding a log statement.
- We increase Y for a change that adds results, changes results, or has the potential to change results. Example: changing a parameter in the GSNAP command.
- We increase X for a paradigm shift in how the pipeline is conceived. Example: adding a de-novo assembly step and then reassigning hits based on the assembled contigs.
Changes to X or Y force recomputation of all results when a sample is rerun using idseq-web. Changes to Z do not force recomputation when the sample is rerun - the pipeline will lazily reuse existing outputs in AWS S3.

When releasing a new version, please add a Git tag of the form `vX.Y.Z`.

- 4.11.9
  - Switch from assertion that cd-hit-clusters only emit one read to a warning if they emit more than one read

- 4.11.8
  - Support local non-host alignment

- 4.11.7
  - Replace references to idseq-database with idseq-public-references

- 4.11.5
  - Unify remote and local non host alignment commands

- 4.11.4
  - Correctly count reads in empty file

- 4.11.3
  - Make STAR outputs deterministic by sorting

- 4.11.2
  - Disable lazy loading remote alignment chunks

- 4.11.1
  - Add support for running alignment locally.
  - Add fetch_from_s3 compatibility hack for transition to miniwdl handled downloads

- 4.11.0
  - Normalize error handling by raising a specific class of exceptions when expected errors occur
  - Simplify idseq_dag.util.count.reads() and remove invalid assumptions
  - Add unit tests for changed components
  - Fix unit test autodiscovery
  - Improve portability of existing unit tests and run them in CI

- 4.10.0
  - Add e-value threshold to require all internal alignments (short read and assembly-based) have e-value below threshold.

- 4.9.0
  - Update NCBI index databases to those downloaded on 2020-04-20.

- 4.8.0
  - Fix counting of reads that have a tax ID but do not have an accession to ensure all reads mapped to taxa only by assembly are counted in r value.

- 4.7.1
  - Fix occasional error in unidentified.fa counting.

- 4.7.0
  - Add a step-level entry point CLI.

- 4.6.0
  - Include non-unique reads in unidentified fastas for download.

- 4.5.0
  - Include non-unique reads in non-host fastqs for download.

- 4.4.2
  - Make the pipeline deterministic with hard coded seeded pseudo random number generation
  - Re-enable reusing previous chunks when re-running alignment.

- 4.4.1
  - Disable reusing previous chunks when re-running alignment. This was casing an error because the input to this step is non-determinsitic and cdhitdup requires all reads from the input to be present in the `duplicate_cluster_sizes` file.

- 4.4
  - When assigning contigs to their best-matching accessions, prioritize the accession that has the most best matches across all contigs when resolving ties.

- 4.3.8
  - Increase gsnap and rapsearch2 chunk size by 4x to reduce the number of batch jobs generated
  - Decreased alignment max chunks in flight from 32 to 16 to better balance chunk execution between large and small jobs

- 4.3.7
  - Remove old alignment code
  - Remove sqlite

- 4.3.6
  - Run alignment chunks first on spot then on on-demand after two failures

- 4.3.5
  - Check for filename collisions on input and output

- 4.3.4
  - Refactor call_hits_m8 to improve memory usage

- 4.3.3
  - Update pipeline stage status directly on s3 file for compatibility with SFN pipeline.

- 4.3.2
  - Clean up log statement for AMR bug.

- 4.3.1
  - Add compatibility for idseq-web environment name 'production'

- 4.3.0
  - Generate betacoronavirus fastq files for user download if use_taxon_whitelist is specified in the DAG.

- 4.2.4
  - Update RunAlignmentRemotely to name batch jobs with the chunk id, project id, and sample id
  - Update RunAlignmentRemotely to download results using boto3 rather than fetch_from_s3

- 4.2.3
  - Validate input step now properly rejects invalid gzip files.

- 4.2.2
  - Fix bug in phylo tree creation for organisms with an unknown superkingdom.

- 4.2.1
  - Switch RunAlignmentRemotely to distribute alignment chunks use AWS Batch instead of custom Autoscaling Group Logic logic

- 4.2.0
  - Apply deuterostome, blacklist, whitelist and human filters to contigs.

- 4.1.1
  - Removed `DAG_SURGERY_HACKS_FOR_READ_COUNTING` and related code.

- 4.1.0
  - Fix the behavior of blacklist and whitelist so that they are applied during compilation of taxon counts and not during identification of candidate accessions. This means that the pipeline now finds the global best taxon match for each read/contig first and only then excludes blacklisted/non-whitelisted taxa from the resulting counts. Previously, it was artifically restricting the search space by applying the blacklist and whitelist to candidate taxa upfront, thus reporting non-optimal matches for affected reads/contigs.

- 4.0.0
  - Switched to computing relative abundance values (r, rPM, contig r) without duplicate removal.

- 3.21.2
  - fix v4 counts for nonhuman hosts human filter steps; no effect in v3

- 3.21.1
  - bugfix affecting samples where SPADES crashed and we only have read alignment data

- 3.21.0
  - work around cdhitdup bug affecting unpaired reads that sometimes discards half the unique reads in an unpaired sample
  - set stage for 4.0 release by changing cdhit identity threshold to 100%
  - emit new files taxon_count_with_dcr.json, duplicate_cluster_sizes.tsv, dedup1.clstr
  - compute ReadCountingMode.COUNT_ALL but still emit COUNT_UNIQUE

- 3.20.1
  - Update s3parcp
  - Switch uploads to s3parcp
  - Support uploading directories with checksum metadata
  - Standardize gsnap and rapsearch index location
  - Update gsnap index generation to upload the raw index directory in addition to a tarball

- 3.20.0
  - Add a custom taxon whitelist mode. Fix taxon blacklist reference downloads.

- 3.19.6
  - Finished removal of optional_files_to_upload

- 3.19.5
  - Switch title of STAR description to be above first line of description

- 3.19.4
  - Copy change for STAR step
  - Don't break STAR if picard can't generate metrics

- 3.19.3
  - Handle case of null nucleotide type for collecting insert metrics

- 3.19.2
  - Use additional_output_files_visible

- 3.19.1
  - Upload additional cdhitdup output.

- 3.19.0
  - Compute insert size metrics for humans only

- 3.18.1
  - Update log statement for AMR bug for alerting purposes.

- 3.18.0
  - Version marker: Update NCBI index databases to those downloaded on 2020-02-10.

- 3.17.0
  - Version marker: Update NCBI index databases to those downloaded on 2020-02-03.

- 3.16.6
  - Isolate directories on alignment instances to chunks rather than whole samples
  - Clean up intermediate files from alignment instances after running alignment on a chunk
  - Ensure alignment instance is clean before running alignment on a chunk

- 3.16.5
  - Fix dag validation of Rapsearch index generation template.

- 3.16.4
  - Fail gracefully with INSUFFICIENT_READS error if all reads drop out during LZW filtering.

- 3.16.3
  - Use s3parcp with checksum for rapsearch index uploads

- 3.16.2
  - fix botocore import issue for util.s3
  - switch to subprocess command for `util.s3.list_s3_keys` for thread safety

- 3.16.1
  - implement LRU policy for reference downloads cache

- 3.16.0
  - Only compute insert size metrics for RNA reads if we have an organism-specific gtf file

- 3.15.1-4
  - change PySAM concurrency pattern to improve performance and eliminate deadlock
  - reduce logging from run_in_subprocess decorator
  - avoid using corrupt reference downloads
  - increase default Rapsearch timeout

- 3.15.0
  - Compute insert size metrics for all hosts for paired end DNA reads
  - Compute insert size metrics for hosts with gtf files for paired end RNA reads

- 3.14.1-4
  - aws credential caching and other stability improvements
  - fix bug that made reverse-strand alignments appear very short in the coverage viz
  - limit blastn BATCH_SIZE to avoid out-of-memory errors (results are unchanged)
  - reduce RAM footprint of blast rerank python to avoid out-of-memory errors (results are unchanged)

- 3.14
  - add average insert size computation

- 3.13.1 - 3.13.3
  - Make phylo tree and alignment viz steps more robust to missing accessions in index.
  - Ensure reference caching respects version.
  - Reduce frequency of s3 requests, other stability fixes.

- 3.13
  - Rerank NT blast results by taking into account all fragments for a given contig and reference sequence, not just the highest scoring one.

- 3.12.1
  - Fix bug where reads unassigned during alignment that were assembled into contigs were also being counted as loose reads.

- 3.12.0
  - Update NCBI databases to those downloaded on 2019-09-17.

- 3.11.0
  - Modify the LZW filter to apply a more stringent cutoff at higher read lengths.

- 3.10.2
  - Better logging for a rare AMR bug.

- 3.10.1
  - Increase GSNAP threads to 48 for better utilization of r5d.metal instances.

- 3.10.0
  - Apply a length filter, requiring all NT alignments (GSNAP and BLAST) be >= 36 nucleotides long.

- 3.9.4
  - Additional performance improvements in run_srst2 step, so that the step uses less RAM.

- 3.9.3
  - Fix typo in phylo tree generation step.

- 3.9.2
  - Fixed error in run_srst2 that failed to take into account different naming patterns from srst2 for the sorted bam file that it outputs.

- 3.9.1
   - Refactoring of command execution patterns and logs.
   - Removed some false error log messages related to lz4 file download support.

- 3.9.0
   - Add number of reads, reads per million, and depth per million to the output of PipelineStepRunSRST2.

- 3.8.0
   - Creates a [status name]_status.json file for each dag_json received, which each step updates with information
     about itself and its status.

- 3.7.6
   - Fail with an informative user error if the input contains broken read pairs.

- 3.7.0 .. 3.7.5
   - Validate whether input files to a pipeline step contain a sufficient number of reads.
     Output invalid_step_input.json file if validation fails.
   - Log output in JSON format. Change TraceLock log level to DEBUG.
   - Upgrade to python 3.7.3
   - Remove db_hack. Standardize db_open/db_assert_table/db_close log entries.
   - Fix division by zero error in coverage viz step.
   - Modify trimmomatic command to reduce MINLEN parameter to 35 and allow reads from fragments with small
     insert sizes (where R1 and R2 are reverse complements of each other) through the QC steps.

- 3.6.6
   - Another fix related to sqlite3 concurrency

- 3.6.0 .. 3.6.5
   - Fix an issue with the log event function when trying to log non json serializable fields.
   - A possible fix to some hanging issues in the pipeline that seem to be related to sqlite3 concurrency.
   - Address array index rounding error in coverage viz.
   - Extra logs to help detecting potential deadlocks in the pipeline
   - Add pipeline step to generate data for coverage visualization for IDseq report page. Data includes an index
     file that maps taxons to accessions with available coverage data, as well as data files for each accession
     that list various metrics including the coverage of the accession.

- 3.5.0 ... 3.5.4
   - New log methods to write log events. Added and replaced a few log entries.
   - Add ability to run STAR further downstream from input validation. This can be used to filter human reads
     after the host has been filtered out (if host is non-human).
   - Handle absence of m8 hits in PipelineStepBlastContigs.
   - Choose most-represented accessions of assembly/gsnap.hitsummary2.tab and assembly/rapsearch2.hitsummary2.tab
     as the NCBI references to include on phylogenetic trees, as opposed to making the choice based on pre-assembly
     align_viz files.
   - Improve the efficiency of S3 downloads and uploads in PipelineStepGenerateAlignmentViz.

- 3.4.0
   - switch from shelve to sqlite3 for all the lookup tables
   - add lineage generation step

- 3.3.2
   - Add input validation step.

- 3.3.1
   - Add step to generate nonhost fastq files by filtering the original fastq files for nonhost reads.

- 3.3.0
   - Upgrade GSNAP executable to version 2018-10-26.  Index remains unchanged at 2018-12-01.
     In comprehensive testing on a diverse set of samples, this has shown just a few minor
     effects on overall results, mostly for reads that align at the limit of detection.
     The benefit of the change is 3x-8x faster performance.  A/B test data is archived
     in slack channel #idseq-benchmarking.

- 3.2.5-3.2.1 only affect staging environment
   - 3.2.5 GSNAP Pre-release 2018-10-26, this time for real.
   - 3.2.4 Revert 3.2.3.
   - 3.2.3 GSNAP Pre-release 2018-10-20 (briefly thought to be 2018-10-20 by mistake).
   - 3.2.2 Revert 3.2.1.
   - 3.2.1 GSNAP Pre-release 2018-10-20.

- 3.2.0
   - Assembly with paired ends if available
   - Coverage Stats Step
- 3.1.0
   - Assembly based pipeline. Add assembly and blast the contigs to the aligned accessions

- 2.11.0
   - Add adapter trimming step.

- 2.10.0
   - Relax LZW filter for reads longer than 150 bp, linearly with read length.

- 2.9.0
   - Change how blacklist filter works so that if a read maps both to
     blacklisted and non-blacklisted entries, it isn't dropped, and
     only the non-blacklisted entries are used.  This improves recall
     for organisms whose DNA is used in cloning vectors, such
     as Chikunguya virus.

- 2.8.0
   - Add taxon blacklist filtering at hit calling

- 2.7.1 ... 2.7.4
   - Reduce LZW runtime from 2h 35m to 24 min on the largest samples
   - Increase GSNAP threads to 36 for i3.metal instances.
   - Addded Antimicrobial resistance step. Results for other steps won't change; only new results for AMR are expected.
   - Acquire lock before fork in run_in_subprocess decorator

- 2.7
   - ?

- 2.6
   - ?

- 2.5
   - ?

- 2.4.0
   - New directed acyclic graph-based execution model for the pipeline. Changes integration with the web app as well.

Below is copied from https://github.com/chanzuckerberg/idseq-pipeline :

- 1.8.7
   - Bug fix for count_reads and non-host read counts.

- 1.8.4 ... 1.8.6
   - Minor code quality, documentation, and logging improvements.

- 1.8.0 ... 1.8.3
   - Upload a status file that indicates when a job has completed.
   - Add a dedicated semaphore for S3 uploads.
   - Code quality and documentation improvements.
   - Restore capability to run non-host alignment from the development environment.
   - Try a more relaxed LZW fraction if the initial filter leaves 0 reads

- 1.7.2 ... 1.7.5
   - General code style changes and code cleanup.
   - Convert string exceptions and generic exceptions to RuntimeErrors.
   - Change some print statements for python3.
   - Add more documentation.

- 1.7.1
   - Truncate enormous inputs to 75 mil paired end / 150 mil unpaired reads.
   - Support input fasta with pre-filtered host, e.g. project NID.
   - Many operational improvements.

- 1.7.0
    - Add capability to  further filter out host reads by filtering all the hits
      from gsnapping host genomes. (i.e. gsnap hg38/patron5 for humans).

- 1.6.3 ... 1.6.1
    - Handle bogus 0-length alignments output by gsnap without crashing.
    - Fix crash for reruns which reuse compatible results from a previous run.
    - Fix crash for samples with unpaired reads.
    - Improve hit calling performance.

- 1.6.0
    - Fix fasta downloads broken by release 1.5.0, making sure only
      hits at the correct level are output in the deduped m8.
    - Fix fasta download for samples with unpaired reads by eliminating
      merged fasta for those samples.
    - Extend the partial fix in release 1.5.1 to repair more of the
      broken reports.  Full fix requires rerun with updated webapp.
    - Correctly aggregate counts for species with unclassified genera,
      such as e.g. genus-less species 1768803 from family 80864.
    - Fix total count in samples with unpaired reads (no longer doubled).
    - Fix crash when zero reads remain after host filtering.
    - Fix bug in enforcing command timeouts that could lead to hangs.
    - Fix performance regression in stage 2 (non-host alignment)
      introduced with 1.5.0.
    - Deduplicate and simplify much of stage 2, and improve performance
      by parallelizing uploads and downloads.

- 1.5.1
    - Fix bug introduced in 1.5.0 breaking samples with non-species-specific
      deuterostome hits.

- 1.5.0
    - Identify hits that match multiple species within the same genus as
      "non species specific" hits to the genus.

- 1.4.0
    - Version result folder.

- 1.3.0
    - Fix bug causing alignment to run before host subtraction in samples
      with unpaired reads.
    - Include ERCC gene counts from STAR.

- 1.2.0
    - Synchronize pair order after STAR to improve sensitivity in 10% of
      samples with paired-end reads.
