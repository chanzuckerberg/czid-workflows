lint:
	find . -name '*.wdl' | xargs miniwdl check

publish:
	scripts/publish.sh
