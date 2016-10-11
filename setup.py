from setuptools import setup, find_packages
from os import system

system('make')

setup(name='abruijn',
      version='0.4b',
      description='Long read assembly via A-Bruijn graph',
      url='https://github.com/fenderglass/ABruijn',
      author='Mikhail Kolmogorov',
      author_email = '',
      license='BSD-3-Clause',
      packages=['abruijn'],
      package_data={'abruijn': ['resource/*.mat']},
      scripts = ['bin/abruijn-assemble', 'bin/abruijn-polish', 'scripts/abruijn'],
      )

