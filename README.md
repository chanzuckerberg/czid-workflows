# idseq-workflows - portable IDseq production pipeline logic

#### Infectious Disease Sequencing Platform
IDseq is a hypothesis-free global software platform that helps scientists identify pathogens in metagenomic sequencing
data.

- **Discover** - Identify the pathogen landscape
- **Detect** - Monitor and review potential outbreaks
- **Decipher** - Find potential infecting organisms in large datasets

IDseq is a collaborative open project of [Chan Zuckerberg Initiative](https://www.chanzuckerberg.com/) and
[Chan Zuckerberg Biohub](https://czbiohub.org).

This repository contains [WDL](https://openwdl.org/) workflows that the [IDseq](https://idseq.net/) platform uses in
production.

# CI/CD

We use GitHub Actions for CI/CD. Lint and unit tests run on GitHub from jobs in `.github/workflows/wdl-ci.yml`
(triggered on every commit).
