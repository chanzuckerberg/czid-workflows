# Makefile for czid-workflows

## OPTIONAL VARIABLES
WORKFLOW?=short-read-mngs # default needed to build dag-test
VERSION?=latest

.PHONY: help
help: 
	@grep -E '^[a-zA-Z_-]+%?:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: build
build: ## Build docker images for a given workflow. eg. make build WORKFLOW=short-read-mngs
	./scripts/docker-build.sh workflows/$(WORKFLOW) -t czid-$(WORKFLOW)

.PHONY: pull
pull: ## Pull docker image from public github repository. Faster than build. Defaults to latest eg. make pull WORKFLOW=long-read-mngs
	docker pull ghcr.io/chanzuckerberg/czid-workflows/czid-$(WORKFLOW)-public:$(VERSION)
	docker tag ghcr.io/chanzuckerberg/czid-workflows/czid-$(WORKFLOW)-public:$(VERSION) czid-$(WORKFLOW)

.PHONY: lint
lint: ## lint files
	pre-commit run --all-files
	flake8 .
	flake8 --ignore E302,F401,E501,W504,E711,E712,E722,E741 lib/idseq-dag/idseq_dag

.PHONY: release
release:
	scripts/release.sh

test-%: ## run miniwdl step tests eg. make test-short-read-mngs
	DOCKER_IMAGE_ID=czid-$* pytest -v -n `python3 -c 'import multiprocessing as mp; print(max(1,mp.cpu_count()-1))'` --durations=0 --tb=short --log-cli-level=11 workflows/$*/test

integration-test-%: ## run miniwdl integration tests eg. make integration-test-short-read-mngs
	DOCKER_IMAGE_ID=czid-$* pytest -v -n `python3 -c 'import multiprocessing as mp; print(max(1,mp.cpu_count()-1))'` --durations=0 --tb=short --log-cli-level=11 workflows/$*/integration_test
	
.PHONY: test
test: ## run miniwdl step tests for all workflows eg. make test
	for i in $$(ls workflows); do $(MAKE) test-$$i; done

.PHONY: dag-test
dag-test: build ## run tests for idseq-dag
	docker run -it -v ${PWD}/lib/idseq-dag:/work -w /work czid-$(WORKFLOWS) pytest -s
