#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs repeat analyser binary
"""

import subprocess
import logging
import os

from abruijn.utils import which

REPEAT_BIN = "abruijn-repeat"
logger = logging.getLogger()


class RepeatException(Exception):
    pass


def check_binaries():
    if not which(REPEAT_BIN):
        raise RepeatException("Repeat binary was not found. "
                              "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([REPEAT_BIN, "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        raise RepeatException("Some error inside native {0} module: {1}"
                                .format(REPEAT_BIN, e))


def analyse_repeats(args, input_assembly, out_folder, log_file):
    logger.debug("-----Begin repeat analyser log------")
    cmdline = [REPEAT_BIN, "-k", str(args.kmer_size), "-l", log_file,
               "-t", str(args.threads), "-v", str(args.min_overlap)]
    if args.debug:
        cmdline.append("-d")
    cmdline.extend([input_assembly, args.reads, out_folder])

    try:
        subprocess.check_call(cmdline)
    except (subprocess.CalledProcessError, OSError) as e:
        raise RepeatException("Error in repeat binary: " + str(e))
