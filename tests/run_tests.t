#!/bin/bash
# bash-tap harness to trigger test WDL workflows. The test workflows (also in this directory) run
# the idseq workflows on various inputs using files supplied from command-lines here, and include
# validity assertions. Error message assertions are checked here (as WDL doesn't have an exception
# handling mechanism).
set -o pipefail

# bash-tap boilerplate; see: https://github.com/illusori/bash-tap
BASH_TAP_ROOT="tests/util/bash-tap"
source tests/util/bash-tap/bash-tap-bootstrap

# Build idseq docker images. Currently requires checkout of idseq monorepo alongside (FIXME)
cd "$(dirname $0)/.."
REPO="$(pwd)"
docker build -t idseq_workflows_tests:latest "${REPO}/../idseq/workflows"

# provision temp directory for workflow runs
RUNDIR=$(mktemp -d --tmpdir idseq_tests_XXXXXX)
miniwdl_run="tests/util/miniwdl_run_with_custom_error.py --verbose --dir $RUNDIR"

# declare expected total number of bash-tap "is" assertions to follow
plan tests 3

###################################################################################################
# host_filter
###################################################################################################

# positive cases: run through a testing workflow
run_host_filter_test="$miniwdl_run ${REPO}/tests/host_filter/test.wdl \
    docker_image_id=idseq_workflows_tests:latest \
    test_fastq_gz=${REPO}/tests/host_filter/test.fastq.gz \
    run_validate_input_multiline_fastq=${REPO}/tests/host_filter/run_validate_input_multiline.fastq"

$run_host_filter_test
is "$?" "0"

# error cases: capture stdout JSON and check error messages with jq & bash-tap assertion
$run_host_filter_test invalid_fastq=true | tee ${RUNDIR}/stdout
is "$?" "1"
is "$(jq -r .cause.info.cause ${RUNDIR}/stdout)" "File does not follow fastq format"
