# czid-workflows - portable CZ ID production pipeline logic

#### Infectious Disease Sequencing Platform
CZ ID is a hypothesis-free global software platform that helps scientists identify pathogens in metagenomic sequencing
data.

- **Discover** - Identify the pathogen landscape
- **Detect** - Monitor and review potential outbreaks
- **Decipher** - Find potential infecting organisms in large datasets

CZ ID is a collaborative open project of [Chan Zuckerberg Initiative](https://www.chanzuckerberg.com/) and
[Chan Zuckerberg Biohub](https://czbiohub.org).

## Workflows
Currently we have 5 main workflows. The details of each pipeline are in a README in each of the workflow folders. 

* [short-read-mngs](workflows/short-read-mngs/README.md) 
* [consensus-genome](workflows/consensus-genome/README.md)
* [phylotree-ng](workflows/phylotree-ng/README.md)
* long-read-mngs (Beta)
* amr (Beta)

## Running these workflows
This repository contains [WDL](https://openwdl.org/) workflows that the [CZ ID](https://czid.org) platform uses in
production. See [Running WDL workflows locally](https://github.com/chanzuckerberg/czid-workflows/wiki/Running-WDL-workflows-locally)
to get started with them.

### System Requirements 
The system requirements differs for each workflow and depending on the database being used. For example running the short-read-mngs workflow with the full NT and NR databases would require an instance with >1TB of disk space and >100GB of memory. Running other workflows (e.g. consensus-genome, amr) requires much less space. 

### Software requirements
* docker with buildx support (version >= 19.03)
* python3 
* virtualenv
* requirements-dev.txt - to automatically install this run `make python-dependencies`

### Quick Setup
To get setup, first set the workflow you want to run with 

```export WORKFLOW=<workflow-name>``` e.g.

```export WORKFLOW=amr```

You can see available workflows with `make ls`

Either `build` or `pull` the workflow docker image with 

```make pull  ## The faster option``` or 

```make build ## The slower option, but necessary if you're modifying the docker container```

### Running a workflow
Run a workflow with 

```make run```

Which simply runs the ```miniwdl run path_to_wdl.wdl``` command with some defaults

Each workflow has a number of required and optional inputs, and all require at least an input file (usually a fastq). Default inputs are set from the `workflows/<workflow-name>/test/local_test.yml` file. These may or may not be accurate for every analysis. You can override these defaults and add your own with:

```make run INPUT='-i your_file_here.yml'```

If you're happy with the defaults, you can add arguments to the `miniwdl` command using 

```make run EXTRA_INPUTS='input_fastq=/path/to/input.fastq' ```

### An example
Lets say I want to run a consensus-genome workflow. I would run the following:

```
export WORKFLOW=consensus-genome
make pull # pull the latest docker container from github packages
make python-dependencies # create a .venv and install the requirements-dev.txt dependencies 

make run EXTRA_INPUTS='fastqs_0=workflows/consensus-genome/test/sample_sars-cov-2_paired_r1.fastq.gz \
        fastqs_1=workflows/consensus-genome/test/sample_sars-cov-2_paired_r2.fastq.gz \
        technology="Illumina" \
        sample="my_sample_name"'
```
## CI/CD

We use GitHub Actions for CI/CD. Lint and unit tests run on GitHub from jobs in `.github/workflows/wdl-ci.yml`
(triggered on every commit).

## Contributing

This project adheres to the Contributor Covenant code of conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to opensource@chanzuckerberg.com.

## Security

Please disclose security issues responsibly by contacting security@chanzuckerberg.com.
