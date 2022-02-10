# czid-workflows - portable CZ ID production pipeline logic

#### Infectious Disease Sequencing Platform
CZ ID is a hypothesis-free global software platform that helps scientists identify pathogens in metagenomic sequencing
data.

- **Discover** - Identify the pathogen landscape
- **Detect** - Monitor and review potential outbreaks
- **Decipher** - Find potential infecting organisms in large datasets

CZ ID is a collaborative open project of [Chan Zuckerberg Initiative](https://www.chanzuckerberg.com/) and
[Chan Zuckerberg Biohub](https://czbiohub.org).

## Running these workflows
This repository contains [WDL](https://openwdl.org/) workflows that the [CZ ID](https://czid.org) platform uses in
production. See [Running WDL workflows locally](https://github.com/chanzuckerberg/czid-workflows/wiki/Running-WDL-workflows-locally)
to get started with them.

## CI/CD

We use GitHub Actions for CI/CD. Lint and unit tests run on GitHub from jobs in `.github/workflows/wdl-ci.yml`
(triggered on every commit).

## Contributing

This project adheres to the Contributor Covenant code of conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to opensource@chanzuckerberg.com.

## Security

Please disclose security issues responsibly by contacting security@chanzuckerberg.com.
