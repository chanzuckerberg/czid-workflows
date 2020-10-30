# Running Consensus-Genome locally

For Consensus Genome workflow we can follow a similar workflow to the `short-read-mngs` presented in wiki: [Running-WDL-workflows-locally](https://github.com/chanzuckerberg/idseq-workflows/wiki/Running-WDL-workflows-locally).

We first build a local image with the consensus genome workflow:

```bash
docker build -t idseq-consensus-genome idseq-workflows/consensus-genome
```

We then use our local sample configuration file that points to IDseq's public references for smaller runs:

```bash
miniwdl run --verbose consensus-genome/run.wdl \
    docker_image_id=idseq-consensus-genome \
    fastqs_0=idseq-workflows/tests/consensus-genome/sample_sars-cov-2_paired_r1.fastq.gz \
    fastqs_1=idseq-workflows/tests/consensus-genome/sample_sars-cov-2_paired_r2.fastq.gz \
    sample=sample_sars-cov-2_paired \
    -i idseq-workflows/tests/consensus-genome/local_test.yml
```

Where:

* `docker_image_id=` should be set to the docker image tag you used when building the image (in our example, `idseq-consensus-genone`)
* `consensus-genome/run.wdl` is the WDL for the consensus genome sequencing workflow.
* `fastqs_0` and `fastqs_1` are the pair of FASTQ files. The ones referred to are small files to run locally.
* `sample` is the name to use where referencing the sample in the output files.
* `local_test.yml` supplies boilerplate workflow inputs, such as the S3 paths for the reference databases. For local run purposes, we use lighter references:
  * The human database for host removal only contains chromosome 1.
  * The kraken db used locally only has coronavirus sequences.
