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

#### Files and Databases
In the local_test.yml file and within the `wdl` files, there are references to some default files and databases. A description of some of these files is below:  

Filename | Description
---------|------------
s3://czid-public-references/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/human_STAR_genome.tar | The database used to filter human reads using STAR. Other host genomes can be found at s3://czid-public-references/host_filter/
s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nt | The NT database downloaded from NCBI
s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nr | The NR database downloaded from NCBI
s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/accession2taxid.marisa | A mapping from accession to tax id generated from the NCBI databases using the `marisa-trie` package
s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nt_loc.marisa | A mapping from accession to the location of the accession in NT. Generated using the `marisa-trie` package
s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nr_loc.marisa | A mapping from accession to the location of the accession in NR. Generated using the `marisa-trie` package
s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nt_info.marisa | A mapping from accession to the name and length of the entry generated using the `marisa-trie` package

#### NCBI Indexes
We generate our databases regularly from NCBI. The files are available at `s3://czid-public-references/ncbi-indexes-prod/{date}/index-generation-2/` 

