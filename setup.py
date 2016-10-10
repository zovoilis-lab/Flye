from setuptools import setup

setup(name='abruijn',
      version='0.4b',
      description='Long read assembly via A-Bruijn graph',
      url='https://github.com/fenderglass/ABruijn',
      author='Mikhail Kolmogorov',
      author_email = '',
      license='BSD-3-Clause',
      scripts = ['abruijn.py', 'bin/abruijn-assemble', 'bin/abruijn-polish'],
      data_files=['resource/nano_homopolymers.mat', 'resource/nano_substitutions.mat', 'resource/p6c4_homopolymers.mat', 'resource/p6c4_substitutions.mat', 'resource/pacbio_homopolymers.mat', 'resource/pacbio_substitutions.mat']
      )

