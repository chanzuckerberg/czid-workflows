# short-read-mngs auto benchmarks

This subdirectory has scripts and reference materials for "drift detection" in short-read-mngs pipeline results.

1. **benchmarks.yml**: catalogue of benchmark samples (FASTQs), reference databases, and pipeline settings, defining the test scenarios
2. **run_local.py**: run one or more of the scenarios by locally invoking `miniwdl run` using the WDL code in the current repo checkout
3. **run_dev.py**: submit one or more of the scenarios to the idseq-dev SFN-WDL backend, with a specified point version of the pipeline (requires the invoking session preconfigured with an appropriate AWS profile)
4. **harvest.py**: consolidate results of either runner script into a JSON file
5. **ref_libs/**: library of "reference" (expected) results
6. **short-read-mngs-benchmarks.ipynb**: Jupyter notebook template which generates comparison tables of harvested & reference results

A GitHub Actions workflow runs limited tests on every code push (two small samples & viral reference databases), automatically generating the notebook as a build artifact.

## Manual steps to run

### (1A) run_local

### (1B) run_dev

### (2) harvest

### (3) jupyter

## Updating reference library
