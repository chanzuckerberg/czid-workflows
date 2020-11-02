import os
from ruamel.yaml import YAML
from pathlib import Path

def load_benchmarks_yml():
    with open(Path(__file__).parent.parent / "benchmarks.yml") as benchmarks_yml:
        return YAML(typ="safe").load(benchmarks_yml.read())
