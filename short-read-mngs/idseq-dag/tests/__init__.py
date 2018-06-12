import unittest

from .basic import MyFirstTest
from .advanced import MySecondTest
from .run_bowtie2 import RunBowtie2Test
import idseq_dag.util.log as log

if __name__ == '__main__':
  log.configure_logger()
  unittest.main(verbosity=2)


