lint:
	find . -name '*.wdl' | xargs miniwdl check

publish:
	scripts/publish.sh

test:
	prove -v tests/*.t
