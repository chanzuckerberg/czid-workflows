from setuptools import setup, find_packages
from idseq_utils import __version__

setup(name='idseq_utils',
      version=__version__,
      description='standalone helper functions for the idseq pipeline',
      url='',
      author='IdSeq Team @ Chan Zuckerberg Initiative',
      author_email='idseqhelp@chanzuckerberg.com',
      license='MIT',
      packages=find_packages(exclude=["tests.*", "tests"]),
      install_requires=["pytz", "boto3", "biopython"],
      tests_require=["coverage", "flake8", "wheel"],
      dependency_links=[],
      entry_points={
          'console_scripts': [
              'sync-pairs = idseq_utils.sync_pairs:main',
              'count-reads = idseq_utils.count_reads:main',
              'save-description = idseq_utils.save_descriptions:main'
          ]
      },
      zip_safe=False)
