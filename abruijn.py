#!/usr/bin/env python

#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This script does all the necessary preparations
and invokes ABruijn
"""

import os
import sys

BIN_DIR = "bin"
RESOURCE_DIR = "resource"

#Check Python version
if sys.version_info[:2] != (2, 7):
    print("Error: ABruijn requires Python version 2.7 ({0}.{1} detected)."
          .format(sys.version_info[0], sys.version_info[1]))
    sys.exit(-1)

#Setting executable paths
abruijn_root = os.path.dirname(os.path.realpath(__file__))
bin_absolute = os.path.join(abruijn_root, BIN_DIR)
resource_absolute = os.path.join(abruijn_root, RESOURCE_DIR)
#sys.path.insert(0, bin_absolute)
sys.path.insert(0, abruijn_root)
os.environ["PATH"] = bin_absolute + os.pathsep + os.environ["PATH"]
os.environ["ABRUIJN_RES"] = resource_absolute

#ABruijn entry point
from abruijn.main import main
sys.exit(main())
