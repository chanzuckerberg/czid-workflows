lint:
	pre-commit run --all-files

publish:
	scripts/publish.sh

test:
	prove -v tests/*.t

.PHONY: lint publish test
