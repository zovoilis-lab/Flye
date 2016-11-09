#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

from distutils.core import setup
from distutils.command.build import build as DistutilsBuild
import subprocess

class MakeBuild(DistutilsBuild):
    def run(self):
        try:
            subprocess.check_call(['make'])
        except subprocess.CalledProcessError as e:
            print "Compilation error: ", e
            return
        DistutilsBuild.run(self)

setup(name='abruijn',
      version='0.4b',
      description='Long read assembly via A-Bruijn graph',
      url='https://github.com/fenderglass/ABruijn',
      author='Mikhail Kolmogorov',
      author_email = '',
      license='BSD-3-Clause',
      packages=['abruijn'],
      package_data={'abruijn': ['resource/*.mat']},
      scripts = ['bin/abruijn-assemble', 'bin/abruijn-polish', 'bin/abruijn'],
      cmdclass={'build': MakeBuild}
      )
