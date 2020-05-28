#!/bin/bash

set -euo pipefail

for tag in $(git tag); do
    for file in $(git ls-tree -r --name-only "$tag" | grep '.wdl$'); do
        echo "Uploading: $tag $file"
        git show "${tag}:${file}" | aws s3 cp --acl public-read - "s3://idseq-workflows/${tag}/${file}"
    done
done
