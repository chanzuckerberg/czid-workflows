# CZ ID Short-Read MNGS workflow 

CZ ID is a hypothesis-free global software platform that helps scientists identify pathogens in metagenomic sequencing data.

* Discover - Identify the pathogen landscape
* Detect - Monitor and review potential outbreaks
* Decipher - Find potential infecting organisms in large datasets
A collaborative open project of Chan Zuckerberg Initiative and Chan Zuckerberg Biohub.

The documentation presented here reflects the current CZ ID pipeline status and processes, with descriptions of the parameters used. It is divided into three main analysis sections.

1. Host Filtering and QC
2. Alignment and Taxonomic Aggregation
3. Assembly and Realignment
4. Reporting and Visualization
The pipeline relies on a variety of different software, whose versions are listed below. The most up-to-date version information can be found in the Dockerfile.

## Running Short-Read MNGS locally

### Short-Read MNGS viral pipeline dependencies 
* The short-read mngs viral pipeline uses a database that has been restricted to only viral reads. This is only used for end-to-end testing without having to download the entire databse. 

#### System Requirements
* \> 10GB of harddrive space
* 16GB of RAM

#### Dependencies
* Docker 
* miniwdl 
* python3


### Short-Read MNGS full pipeline dependencies

#### System Requirements
* \>1TB of harddrive space 
    * Mostly used to store full NT and NR databases 
* 64GB of RAM
* \>=4 CPUS

#### Dependencies
* Docker
* miniwdl
* python3 


### Runing the Pipeline Locally Step-by-Step
Clone the `czid-workflows` repo with: 

```bash 
git clone https://github.com/chanzuckerberg/czid-workflows.git
cd czid-workflows
```

Either pull the short-read-mngs docker container with: 

```bash 
docker pull ghcr.io/chanzuckerberg/czid-workflows/czid-short-read-mngs-public:latest
docker tag ghcr.io/chanzuckerberg/czid-workflows/czid-short-read-mngs-public:latest czid-short-read-mngs
```

Or build the container from scratch with: 

```bash
./scripts/docker-build.sh ./workflows/short-read-mngs/ -t czid-short-read-mngs
```

Install the miniwdl dependencies with: 

```bash
pip3 install -r requirements-dev.txt
```
If running the viral pipeline, run with: 
```bash 
miniwdl run workflows/short-read-mngs/local_driver.wdl \
    docker_image_id=czid-short-read-mngs \
    fastqs_0=workflows/short-read-mngs/test/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R1.fastq.gz \
    fastqs_1=workflows/short-read-mngs/test/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R2.fastq.gz \
    -i workflows/short-read-mngs/test/local_test_viral.yml --verbose
```

If running the full pipeline, run:
```bash 
miniwdl run workflows/short-read-mngs/local_driver.wdl \
    docker_image_id=czid-short-read-mngs \
    fastqs_0=workflows/short-read-mngs/test/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R1.fastq.gz \
    fastqs_1=workflows/short-read-mngs/test/norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v6__R2.fastq.gz \
    -i workflows/short-read-mngs/test/local_test.yml --verbose
```
Where `fastqs_0` and `fastqs_1` are paired-end fastq files. For single-end reads, just use `fastqs_0`.

The `.yml` file contains default databases for host filtering and alignment as well as supplemental helper files. You can change the values in these files to suit your analysis. 

#### Files and NCBI Indexes
In the local_test.yml file and within the `wdl` files, there are references to some default files and databases. 

For host read filtering, two versions of the human genomes are available: [GRCh38 or HG38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/) and the latest assembly from the Telomere-to-Telomore (T2T) Consortium completed in 2022 [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/). 

GRCh38 | T2T-CHM13
-------|----------
s3://public-test-bucket-idseq/host_filter/human/2022/bowtie2_index_tar/GRCh38_ERCC.bowtie2.tar | s3://czid-public-references/host_filter/human_telomere/2023-07-05/host-genome-generation-1/human_telomere.bowtie2.tar
s3://public-test-bucket-idseq/host_filter/human/2022/hisat2_index_tar/GRCh38_ERCC.hisat2.tar | s3://czid-public-references/host_filter/human_telomere/2023-07-05/host-genome-generation-1/human_telomere.hisat2.tar
s3://public-test-bucket-idseq/host_filter/human/2022/kallisto_idx/GRCh38_ERCC.kallisto.idx | s3://czid-public-references/host_filter/human_telomere/2023-07-05/host-genome-generation-1/human_telomere.kallisto.idx

<br>

We generate our databases regularly from NCBI. The files are available at `s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/`. The `{date}` used on the web application is either `2021-01-22` or `2024-02-06`.

Filename | Description
---------|------------
s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/nt | The NT database downloaded from NCBI
s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/nr | The NR database downloaded from NCBI
s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/accession2taxid.marisa | A mapping from accession to tax id generated from the NCBI databases using the `marisa-trie` package
s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/nt_loc.marisa | A mapping from accession to the location of the accession in NT. Generated using the `marisa-trie` package
s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/nr_loc.marisa | A mapping from accession to the location of the accession in NR. Generated using the `marisa-trie` package
s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/nt_info.marisa | A mapping from accession to the name and length of the entry generated using the `marisa-trie` package



