#!/bin/bash

docker buildx build --build-context lib=lib "$@"
