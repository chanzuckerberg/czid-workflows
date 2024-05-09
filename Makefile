# Makefile for czid-workflows

## OPTIONAL VARIABLES
WORKFLOW?=short-read-mngs# default needed to build dag-test
VERSION?=latest
EXTRA_INPUTS?=
TASK?=

ifeq ($(WORKFLOW),short-read-mngs)
    TEST_FILE?=workflows/$(WORKFLOW)/test/local_test_viral.yml
else
    TEST_FILE?=workflows/$(WORKFLOW)/test/local_test.yml
endif

ifeq ($(TASK),)
ifneq ($(wildcard $(TEST_FILE)), )
	INPUT?=-i $(TEST_FILE)
endif
endif

TASK_CMD := $(if $(TASK), --task $(TASK),)


.PHONY: help
help: 
	@grep -E '^[a-zA-Z_-]+%?:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

.PHONY: build
build: .build-$(WORKFLOW) ## Build docker images for a given workflow. eg. make build WORKFLOW=short-read-mngs

.build-%:
	./scripts/docker-build.sh workflows/$* -t czid-$*
	touch .build-$*
	
.PHONY: rebuild
rebuild: ## Rebuild docker images for a given workflow. If you change anything in the lib/ directory you will likely need to rebuild. 
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
	pytest -v -n $$(( $$(nproc) - 1)) --durations=0 --tb=short --log-cli-level=11 workflows/$*/test

integration-test-%: ## run miniwdl integration tests eg. make integration-test-short-read-mngs
	pytest -v -n $$(( $$(nproc) - 1)) --durations=0 --tb=short --log-cli-level=11 workflows/$*/integration_test

.PHONY: cargo-test-index-generation
cargo-test-index-generation: ## run cargo test for index-generation
	(cd workflows/index-generation/ncbi-compress && cargo test)

.PHONY: test
test: ## run miniwdl step tests for all workflows eg. make test
	for i in $$(ls workflows); do $(MAKE) test-$$i; done

.PHONY: dag-test
dag-test: build ## run tests for idseq-dag
	docker run -it -v ${PWD}/lib/idseq-dag:/work -w /work czid-$(WORKFLOW) pytest -s

.PHONY: python-dependencies
python-dependencies: .python_dependencies_installed # install python dependencies

.python_dependencies_installed: 
	virtualenv -p python3 .venv
	.venv/bin/pip install -r requirements-dev.txt
	echo "Run: source .venv/bin/activate"
	touch .python_dependencies_installed

# TODO: Break this up into 2 functions, one to resolve the RUNFILE and one to run the workflow
.PHONY: run
run: build python-dependencies ## run a miniwdl workflow. eg. make run WORKFLOW=consensus-genome. args: WORKFLOW,EXTRA_INPUT,INPUT,TASK_CMD
	if [ "$(WORKFLOW)" = "short-read-mngs" ]; then \
		RUNFILE="local_driver.wdl"; \
	elif [ $$(ls workflows/$(WORKFLOW)/*.wdl | wc -l) -eq 1 ]; then \
		RUNFILE=$$(ls workflows/$(WORKFLOW)/*.wdl); \
		RUNFILE=$$(basename $$RUNFILE); \
	elif [ -f "workflows/$(WORKFLOW)/$(WORKFLOW).wdl" ]; then \
		RUNFILE="$(WORKFLOW).wdl"; \
	else \
		RUNFILE="run.wdl"; \
	fi; \
	.venv/bin/miniwdl run workflows/$(WORKFLOW)/$$RUNFILE docker_image_id=czid-$(WORKFLOW) $(INPUT) $(EXTRA_INPUTS) $(TASK_CMD)

.PHONY: miniwdl-explore
miniwdl-explore: ## !ADVANCED! explore a previous miniwdl workflow run in the cli. eg. make miniwdl-explore OUTPATH='/mnt/path/to/output/'
	cat $(OUTPATH)/inputs.json | jq ' [values[]] | flatten | .[] | tostring | select(startswith("s3"))' | xargs -I {} aws s3 cp "{}" $(OUTPATH)/work/_miniwdl_inputs/0/
	cat $(OUTPATH)/inputs.json | jq ' [values[]] | flatten | .[] | tostring | select(startswith("/"))' | xargs -I {} cp "{}" $(OUTPATH)/work/_miniwdl_inputs/0/
	docker run -it --entrypoint bash -w /mnt/miniwdl_task_container/work -v$(OUTPATH):/mnt/miniwdl_task_container czid-$(WORKFLOW)

.PHONY: ls
ls: ## list workflows
	@ls -1 workflows/

.PHONY: check
check: python-dependencies ## run miniwdl check on the given workflow
	if [ "$(WORKFLOW)" = "short-read-mngs" ]; then \
                RUNFILE="local_driver.wdl"; \
	elif [ -f "workflows/$(WORKFLOW)/$(WORKFLOW).wdl" ]; then \
		RUNFILE="$(WORKFLOW).wdl"; \
	else \
			RUNFILE="run.wdl"; \
	fi; \
	.venv/bin/miniwdl check workflows/$(WORKFLOW)/$$RUNFILE
