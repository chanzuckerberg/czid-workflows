import unittest
import logging
import os
import sys

DISABLE_LOG_OUTPUT = True

project_root_dir_name = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + "/../../idseq_dag")
print(f"Appending project dir to sys.path: {project_root_dir_name}")
sys.stdout.flush()
sys.path.append(project_root_dir_name)


if DISABLE_LOG_OUTPUT:
    # we don't want any output during unit tests
    logger = logging.getLogger()
    for h in list(logger.handlers):
        logger.removeHandler(h)
