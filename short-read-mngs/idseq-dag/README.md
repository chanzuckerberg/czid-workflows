# [IDseq](https://idseq.net/) &middot; [![GitHub license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://github.com/chanzuckerberg/idseq-web/blob/master/LICENSE) ![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)

![logo](https://assets.idseq.net/Logo_Black.png)

#### Infectious Disease Sequencing Platform
IDseq is an unbiased global software platform that helps scientists identify pathogens in metagenomic sequencing data.

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
There are four basic elements of an IdSeq Dag

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
      "class": "PipelineStepRunAlignmentRemotely",
      "module": "idseq_dag.steps.run_alignment_remotely",
      "additional_files": {
        "lineage_db": "s3://idseq-database/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/taxid-lineages.db",
        "accession2taxid_db": "s3://idseq-database/alignment_data/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/accession2taxid.db"

          ,"deuterostome_db": "s3://idseq-database/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/deuterostome_taxids.txt"

      },
      "additional_attributes": {
        "service": "gsnap",
        "chunks_in_flight": 32,
        "chunk_size": 15000,
        "max_concurrent": 3,
        "environment": "prod"
      }
    },
    {
      "in": ["host_filter_out"],
      "out": "rapsearch2_out",
      "class": "PipelineStepRunAlignmentRemotely",
      "module": "idseq_dag.steps.run_alignment_remotely",
      "additional_files": {
        "lineage_db": "s3://idseq-database/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/taxid-lineages.db",
        "accession2taxid_db": "s3://idseq-database/alignment_data/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/accession2taxid.db"

          ,"deuterostome_db": "s3://idseq-database/taxonomy/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/deuterostome_taxids.txt"

      },
      "additional_attributes": {
        "service": "rapsearch2",
        "chunks_in_flight": 32,
        "chunk_size": 10000,
        "max_concurrent": 6,
        "environment": "prod"
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

- 3.2.3-3.2.1 only on staging environment
   - 3.2.3 GSNAP Pre-release 2018-10-26.
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
