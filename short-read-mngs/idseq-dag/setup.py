from setuptools import setup, find_packages
from idseq_dag import __version__

setup(name='idseq_dag',
      version=__version__,
      description='executing a DAG for idseq pipeline',
      url='http://github.com/chanzuckerberg/idseq-dag',
      author='IdSeq Team @ Chan Zuckerberg Initiative',
      author_email='idseqhelp@chanzuckerberg.com',
      license='MIT',
      packages=find_packages(exclude=["tests.*", "tests"]),
      package_data={'idseq_dag': ['scripts/fastq-fasta-line-validation.awk']},
      install_requires=["pytz", "boto3"],
      tests_require=["coverage", "flake8", "wheel"],
      dependency_links=[],
      entry_points={
          'console_scripts': [
              'idseq_dag = idseq_dag.__main__:main',
              'idseq-dag-run-step = idseq_dag.__main__:run_step'
          ]
      },
      zip_safe=False)
