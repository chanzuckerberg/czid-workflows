lint:
	pre-commit run --all-files

publish:
	scripts/publish.sh

test:
	pytest -v -n 4 --tb=short --log-cli-level=11 tests/
