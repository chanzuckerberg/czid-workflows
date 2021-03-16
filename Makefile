lint:
	pre-commit run --all-files

publish:
	scripts/publish.sh

test-%:
	TEST_DIR=tests/$*; if ! [ -d $$TEST_DIR ]; then TEST_DIR=$*/test; fi; pytest -v -n `python3 -c 'import multiprocessing as mp; print(max(1,mp.cpu_count()-1))'` --tb=short --log-cli-level=11 $$TEST_DIR

test:
	for i in $$(dirname */*.wdl | uniq); do $(MAKE) test-$$i; done

.PHONY: lint publish test
