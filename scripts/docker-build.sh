#!/bin/bash

set -euxo pipefail

docker buildx build --build-context lib=lib "$@"
