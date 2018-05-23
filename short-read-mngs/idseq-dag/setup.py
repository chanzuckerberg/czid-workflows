from setuptools import setup

setup(name='idseq_dag',
      version='0.1',
      description='executing a DAG for idseq pipeline',
      url='http://github.com/chanzuckerberg/idseq-dag',
      author='IdSeq Team @ Chan Zuckerberg Initiative',
      author_email='regger@chanzuckerberg.com',
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
