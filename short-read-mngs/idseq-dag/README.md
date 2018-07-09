# IDSEQ-DAG

## Installation

```
cd idseq-dag; pip3 install -e .
```

## Example

```
idseq_dag examples/host_filter_dag.json

```

## Test

```
cd idseq-dag; python3 -m unittest

```

## Developers
When merging a commit to master, you need to increase the version number in `idseq_dag/__init__.py`:
  - if results are expected to change, increase the 2nd number
  - if results are not expected to change, increase the 3rd number.


## Release notes

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
