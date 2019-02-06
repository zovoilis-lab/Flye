#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

from __future__ import print_function
import sys

#Check Python version
if sys.version_info[:2] != (2, 7):
    print("Error: Flye requires Python version 2.7 ({0}.{1} detected)."
          .format(sys.version_info[0], sys.version_info[1]))
    sys.exit(-1)

from distutils.core import setup
from distutils.command.build import build as DistutilsBuild
from distutils.spawn import find_executable
import subprocess

from flye.__version__ import __version__


class MakeBuild(DistutilsBuild):
    def run(self):
        if not find_executable("make"):
            print ("ERROR: 'make' command is unavailable")
            sys.exit(1)
        try:
            subprocess.check_call(["make"])
        except subprocess.CalledProcessError as e:
            print ("Compilation error: ", e)
            sys.exit(1)
        DistutilsBuild.run(self)

setup(name='flye',
      version=__version__,
      description='Fast and accurate de novo assembler for single molecule sequencing reads',
      url='https://github.com/fenderglass/Flye',
      author='Mikhail Kolmogorov',
      author_email = 'fenderglass@gmail.com',
      license='BSD-3-Clause',
      packages=['flye', 'flye/assembly', 'flye/config', 'flye/polishing',
                'flye/utils', 'flye/repeat_graph', 'flye/short_plasmids',
                'flye/trestle'],
      package_data={'flye': ['config/bin_cfg/*']},
      scripts = ['bin/flye-assemble', 'bin/flye-polish', 'bin/flye-contigger',
                 'bin/flye-repeat', 'bin/flye', 'bin/flye-minimap2'],
      cmdclass={'build': MakeBuild}
      )
