from setuptools import setup, find_packages
from benchmark_helpers import __version__

setup(name='benchmark_helpers',
      version=__version__,
      description='standalone helper functions for the idseq pipeline',
      url='',
      author='CZID Team @ Chan Zuckerberg Initiative',
      author_email='idseqhelp@chanzuckerberg.com',
      license='MIT',
      packages=find_packages(exclude=["tests.*", "tests"]),
      install_requires=["boto3~=1.23.0"],
      tests_require=["coverage", "flake8", "wheel"],
      dependency_links=[],
      entry_points={
          'console_scripts': []
      },
      zip_safe=False)

