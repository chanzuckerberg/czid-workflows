import unittest

from .run_star import RunStarTest
from .run_priceseq import RunPriceSeqTest
from .run_cdhitdup import RunCDHitDupTest
from .run_lzw import RunLZWTest
from .run_bowtie2 import RunBowtie2Test
from .run_subsample import RunSubsampleTest
from .run_gsnap_filter import RunGsnapFilterTest
from .generate_taxid_fasta import GenerateTaxidFastaTest
from .generate_taxid_locator import GenerateTaxidLocatorTest

import idseq_dag.util.log as log

log.configure_logger()

if __name__ == '__main__':
    unittest.main(verbosity=2)
