# Makefile for czid-workflows

## OPTIONAL VARIABLES
WORKFLOW?=short-read-mngs# default needed to build dag-test
VERSION?=latest
EXTRA_INPUTS?=
ifeq ($(WORKFLOW),short-read-mngs)
    MINIWDL_INPUT?=-i workflows/$(WORKFLOW)/test/local_test_viral.yml
else
    MINIWDL_INPUT?=-i workflows/$(WORKFLOW)/test/local_test.yml
endif

.PHONY: help
help: 
	@grep -E '^[a-zA-Z_-]+%?:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

.PHONY: build
build: .build-$(WORKFLOW) ## Build docker images for a given workflow. eg. make build WORKFLOW=short-read-mngs

.build-%:
	./scripts/docker-build.sh workflows/$* -t czid-$*
	touch .build-$*
	
.PHONY: rebuild
rebuild: 
	rm .build-$(WORKFLOW) | true
	$(MAKE) build WORKFLOW=$(WORKFLOW)

.PHONY: pull
pull: .pull-$(WORKFLOW) ## Pull docker image from public github repository. Faster than build. Possibly less accurate. Defaults to latest eg. make pull WORKFLOW=long-read-mngs

.pull-%: 
	docker pull ghcr.io/chanzuckerberg/czid-workflows/czid-$(WORKFLOW)-public:$(VERSION)
	docker tag ghcr.io/chanzuckerberg/czid-workflows/czid-$(WORKFLOW)-public:$(VERSION) czid-$(WORKFLOW)
	touch .build-$*

.PHONY: lint
lint: ## lint files
	pre-commit run --all-files
	flake8 .
	flake8 --ignore E302,F401,E501,W504,E711,E712,E722,E741 lib/idseq-dag/idseq_dag

.PHONY: release
release:
	scripts/release.sh

test-%: ## run miniwdl step tests eg. make test-short-read-mngs
	pytest -v -n `python3 -c 'import multiprocessing as mp; print(max(1,mp.cpu_count()-1))'` --durations=0 --tb=short --log-cli-level=11 workflows/$*/test

integration-test-%: ## run miniwdl integration tests eg. make integration-test-short-read-mngs
	pytest -v -n `python3 -c 'import multiprocessing as mp; print(max(1,mp.cpu_count()-1))'` --durations=0 --tb=short --log-cli-level=11 workflows/$*/integration_test
	
.PHONY: test
test: ## run miniwdl step tests for all workflows eg. make test
	for i in $$(ls workflows); do $(MAKE) test-$$i; done

.PHONY: dag-test
dag-test: build ## run tests for idseq-dag
	docker run -it -v ${PWD}/lib/idseq-dag:/work -w /work czid-$(WORKFLOW) pytest -s

.PHONY: python-dependencies
python-dependencies: .python_dependencies_installed

.python_dependencies_installed: 
	virtualenv -p python3 .venv
	.venv/bin/pip install -r requirements-dev.txt
	echo "Run: source .venv/bin/activate"
	touch .python_dependencies_installed

.PHONY: run
run: build python-dependencies 
	if [ "$(WORKFLOW)" = "short-read-mngs" ]; then \
		RUNFILE="local_driver.wdl"; \
	else \
		RUNFILE="run.wdl"; \
	fi; \
	.venv/bin/miniwdl run workflows/$(WORKFLOW)/$$RUNFILE docker_image_id=czid-$(WORKFLOW) $(EXTRA_INPUTS) $(MINIWDL_INPUT)
