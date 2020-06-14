SHELL=/bin/bash

lint:
	python3 setup.py flake8

test: lint
	python3 -m unittest discover -s tests/unit -v

install:
	-rm -rf dist
	python setup.py bdist_wheel
	pip install --upgrade dist/*.whl

.PHONY: lint test install
