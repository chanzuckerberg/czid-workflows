SHELL=/bin/bash

lint:
	python setup.py flake8
	flake8 scripts/*

test: lint
	python -m unittest discover -s tests/unit -v

install:
	-rm -rf dist
	python setup.py bdist_wheel
	pip install --upgrade dist/*.whl

.PHONY: lint test install
