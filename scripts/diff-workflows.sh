#!/bin/bash

set -euo pipefail

WORKFLOWS=$(git diff --name-status HEAD^ workflows/ | cut -d/ -f2)

if [ -n "$(git diff --name-status HEAD^ lib)" ]; then 
    # if the lib/ folder changes, run tests for anything that depends on it
    WORKFLOWS+='\n'"short-read-mngs"
    WORKFLOWS+='\n'"long-read-mngs"
fi

echo -e "$WORKFLOWS" | sort | uniq 