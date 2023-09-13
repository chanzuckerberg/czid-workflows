#!/bin/bash

set -euo pipefail

if [[ $# != 3 ]]; then
    echo "This script creates and pushes a git tag representing a new workflow release."
    echo "Usage: $(basename $0) workflow_name release_type release_notes"
    echo "Example: $(basename $0) short-read-mngs patch 'Fixing a bug'"
    exit 1
fi

export WORKFLOW_NAME=$1 RELEASE_TYPE=$2 RELEASE_NOTES="$3"

if ! [[ -d "$(dirname $0)/../workflows/$WORKFLOW_NAME" ]]; then
    echo "Unable to locate workflow under a top level directory workflows/$WORKFLOW_NAME" 1>&2
    exit 1
fi

OLD_TAG=$(git describe --tags --match "${WORKFLOW_NAME}-v*" || echo "${WORKFLOW_NAME}-v0.0.0")
# if [[ $RELEASE_TYPE == major ]]; then
#     TAG=$(echo "$OLD_TAG" | perl -ne '/(.+)-v(\d+)\.(\d+)\.(\d+)/; print "$1-v@{[$2+1]}.0.0"')
# elif [[ $RELEASE_TYPE == minor ]]; then
#     TAG=$(echo "$OLD_TAG" | perl -ne '/(.+)-v(\d+)\.(\d+)\.(\d+)/; print "$1-v$2.@{[$3+1]}.0"')
# elif [[ $RELEASE_TYPE == patch ]]; then
#     TAG=$(echo "$OLD_TAG" | perl -ne '/(.+)-v(\d+)\.(\d+)\.(\d+)/; print "$1-v$2.$3.@{[$4+1]}"')
# else
#     echo "RELEASE_TYPE should be one of major, minor, patch" 1>&2
#     exit 1
# fi

TAG="amr-v1.3.0-beta.5"

if [[ $( git branch --show-current) != "main" ]]; then 
    COMMIT=$(git rev-parse --short HEAD)
    TAG=$TAG"-$COMMIT"
fi

TAG_MSG=$(mktemp)
echo "# Changes for ${TAG} ($(date +%Y-%m-%d))" > $TAG_MSG
echo "$RELEASE_NOTES" >> $TAG_MSG
git log --pretty=format:%s ${OLD_TAG}..HEAD >> $TAG_MSG || true
git config --get user.name || git config user.name "CZ ID release action triggered by ${GITHUB_ACTOR:-$(whoami)}"
git tag --annotate --file $TAG_MSG "$TAG"
git push origin "$TAG"

export TAG DEPLOYMENT_ENV=wdl DEPLOY_REF=main GITHUB_TOKEN=${GH_DEPLOY_TOKEN:-GITHUB_TOKEN}
deployment_args=$(jq -n ".auto_merge=false | .ref=env.DEPLOY_REF | .environment=env.DEPLOYMENT_ENV | .required_contexts=[] | .payload={workflow_tag: env.TAG}")
gh api repos/:owner/idseq/deployments --input - <<< "$deployment_args"
