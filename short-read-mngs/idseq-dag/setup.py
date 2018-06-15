from setuptools import setup
from idseq_dag import __version__

setup(name='idseq_dag',
      version=__version__,
      description='executing a DAG for idseq pipeline',
      url='http://github.com/chanzuckerberg/idseq-dag',
      author='IdSeq Team @ Chan Zuckerberg Initiative',
      author_email='idseqhelp@chanzuckerberg.com',
      license='MIT',
      packages=['idseq_dag'],
      install_requires=[],
      dependency_links=[],
      entry_points={
        'console_scripts': [
          'idseq_dag = idseq_dag.main:main'
        ]
      },
      zip_safe=False)
