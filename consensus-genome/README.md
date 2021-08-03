# IDseq Consensus Genome workflow

This workflow performs reference-based consensus genome mapping from metagenomic sequencing assays with spiked primer
enrichment or from amplicon sample sequencing assays.

The workflow has several modes:

- Oxford Nanopore SARS-CoV-2 samples
- Illumina SARS-CoV-2 samples
- Illumina samples of other viruses
- Illumina reads isolated from metagenomic samples via the IDseq mngs workflow

Based on original work at:
- CZ Biohub SARS-CoV-2 pipeline, https://github.com/czbiohub/sc2-illumina-pipeline
- ARTIC Oxford Nanopore MinION SARS-CoV-2 SOP, https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html

With enhancements and additional modules by the CZI Infectious Disease team.

## Running Consensus-Genome locally

For Consensus Genome workflow we can follow a similar workflow to the `short-read-mngs` presented in wiki: [Running-WDL-workflows-locally](https://github.com/chanzuckerberg/idseq-workflows/wiki/Running-WDL-workflows-locally).

We first build a local Docker container image with the consensus genome workflow:

```bash
docker build -t idseq-consensus-genome idseq-workflows/consensus-genome
```

TIPS: For more detailed setup information
 - for miniwdl - https://github.com/openwdl/learn-wdl/blob/master/6_miniwdl_course/0_setup.md
 - for this example - https://github.com/openwdl/learn-wdl/blob/master/6_miniwdl_course/5c_cloud_spec_consensus-genome.md 

## Run 

We then use our local sample configuration file that points to IDseq's public references for smaller runs:

```bash
miniwdl run --verbose idseq-workflows/consensus-genome/run.wdl \
    docker_image_id=idseq-consensus-genome \
    fastqs_0=idseq-workflows/consensus-genome/test/sample_sars-cov-2_paired_r1.fastq.gz \
    fastqs_1=idseq-workflows/consensus-genome/test/sample_sars-cov-2_paired_r2.fastq.gz \
    sample=sample_sars-cov-2_paired \
    technology=Illumina \
    ref_fasta=s3://idseq-public-references/consensus-genome/MN908947.3.fa \
    -i idseq-workflows/consensus-genome/test/local_test.yml
```

Where:

* `docker_image_id=` should be set to the docker image tag you used when building the image (in our example, `idseq-consensus-genome`)
* `idseq-workflows/consensus-genome/run.wdl` is the WDL for the consensus genome sequencing workflow.
* `fastqs_0` and `fastqs_1` are the pair of FASTQ files. The ones referred to are small files to run locally.
* `sample` is the name to use where referencing the sample in the output files.
* `technology` is the sequencing technology (options = Illumina or ONT)
* `local_test.yml` supplies boilerplate workflow inputs, such as the S3 paths for the reference databases. For local run purposes, we use lighter references:
  * The human database for host removal only contains chromosome 1.
  * The kraken db used locally only has coronavirus sequences.

## Reference files
Filename | Provenance
---------|-----------
s3://idseq-public-references/consensus-genome/MN908947.3.fa | Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/MN908947 in July 2020
s3://idseq-public-references/consensus-genome/ampliseq_primers.bed | The .bed file was obtained from the Illumina Ampliseq protocol documentation https://www.illumina.com/products/by-brand/ampliseq/community-panels/sars-cov-2.html on 2021-01-26
s3://idseq-public-references/consensus-genome/artic_v3_primers.bed | The .bed file was obtained from the CZ Biohub sc2 pipeline repository: https://raw.githubusercontent.com/czbiohub/sc2-illumina-pipeline/master/data/nCoV-2019.bed in July 2020. The master file can be downloaded from ARTIC network here: https://github.com/artic-network/fieldbioinformatics/blob/master/test-data/primer-schemes/nCoV-2019/V3/nCoV-2019.bed
s3://idseq-public-references/consensus-genome/artic_v3_short_275_primers.bed | The .bed file was received from scientists at UCSF on 2021-03-11 and links to this protocol https://www.protocols.io/view/covid-19-artic-v3-illumina-library-construction-an-bh4zj8x6
s3://idseq-public-references/consensus-genome/combined_msspe_artic_primers.bed | The .bed file was obtained from scientists at the CZ Biohub on 2021-01-26.
s3://idseq-public-references/consensus-genome/ercc_sequences.fasta | ERCC sequence file was obtained from the CZ Biohub sc2 pipeline repository: https://github.com/czbiohub/sc2-illumina-pipeline/blob/cd37a25cdf3c0260d082bd0146daa5e192704893/data/ercc_sequences.fasta in July 2020. The initial sequences can be downloaded from ThermoFisher here: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
s3://idseq-public-references/consensus-genome/hg38.fa.gz | The human genome was downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz in July 2020.
s3://idseq-public-references/consensus-genome/human_chr1.fa | Test file was obtained from the CZ Biohub sc2 pipeline repository: https://github.com/czbiohub/sc2-illumina-pipeline/blob/master/data/human_chr1.fa in July 2020.
s3://idseq-public-references/consensus-genome/kraken2_h+v_20200319.tar.gz | Kraken2 database of sars-cov-2 + human that was downloaded from https://genexa.ch/sars2-bioinformatics-resources/ in July 2020, however this resource is no longer kept up-to-date.
s3://idseq-public-references/consensus-genome/kraken_coronavirus_db_only.tar.gz | The smaller kraken2 database was taken from the CZ Biohub sc2 pipeline configuration in July 2020.
s3://idseq-public-references/consensus-genome/msspe_primers.bed | The primer .bed file was initially generated by scientists at CZ Biohub and was taken from the Biohub sc2 pipeline https://raw.githubusercontent.com/czbiohub/sc2-illumina-pipeline/master/data/SARS-COV-2_spikePrimers.bed in July 2020. 
s3://idseq-public-references/consensus-genome/msspe_primers-v2.bed | The primer .bed file was updated by scientists at CZ Biohub to flip the orientation of the primer sequences.
s3://idseq-public-references/consensus-genome/snap_primers.bed | The .bed file was obtained from the swift representatives in December 2020. It can also be found here https://swiftbiosci.com/wp-content/uploads/2020/09/sarscov2_v1_masterfile.txt.zip
s3://idseq-public-references/consensus-genome/vadr-models-corona-1.1.3-1.tar.gz | Downloaded from https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/coronaviridae/1.2-1/vadr-models-corona-1.2-1.tar.gz on 2021-03-05
s3://idseq-public-references/consensus-genome/artic-primer-schemes.tar.gz | `primer_schemes` directory of https://github.com/artic-network/artic-ncov2019/commit/7e359dae37d894b40ae7e35c3582f14244ef4d36
`test/MT007544.fastq.gz` | Copied from https://github.com/artic-network/fieldbioinformatics/blob/master/test-data/MT007544/MT007544.fastq on 2021-03-06

## More Information

For more information, including a screencast of this example, see the `learn-miniwdl` open source course
- Screencast at - https://www.youtube.com/watch?v=bnXOoPm_F2I
- Miniwdl Course at - https://github.com/openwdl/learn-wdl/tree/master/6_miniwdl_course
- WDL course at - https://github.com/openwdl/learn-wdl
