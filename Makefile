lint:
	pre-commit run --all-files
	flake8 .
	flake8 --ignore E302,F401,E501,W504,E711,E712,E722,E741 short-read-mngs/idseq-dag/idseq_dag

publish:
	scripts/publish.sh

test-%:
	pytest -v -n `python3 -c 'import multiprocessing as mp; print(max(1,mp.cpu_count()-1))'` --durations=0 --tb=short --log-cli-level=11 $*/test

test:
	for i in $$(dirname */*.wdl | uniq); do $(MAKE) test-$$i; done

.PHONY: lint publish test
