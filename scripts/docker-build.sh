#!/bin/bash

set -euxo pipefail

docker buildx build --platform linux/amd64 --build-context lib=lib "$@"
