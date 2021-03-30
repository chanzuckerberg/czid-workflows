lint:
	pre-commit run --all-files

publish:
	scripts/publish.sh

test-%:
	pytest -v -n `python3 -c 'import multiprocessing as mp; print(max(1,mp.cpu_count()-1))'` --tb=short --log-cli-level=11 $*/test

test:
	for i in $$(dirname */*.wdl | uniq); do $(MAKE) test-$$i; done

.PHONY: lint publish test
