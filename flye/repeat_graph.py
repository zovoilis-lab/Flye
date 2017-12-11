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
    except (subprocess.CalledProcessError, OSError) as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise RepeatException(str(e))


def analyse_repeats(args, input_assembly, out_folder, log_file, config_file):
    logger.debug("-----Begin repeat analyser log------")
    cmdline = [REPEAT_BIN, "-k", str(args.kmer_size), "-l", log_file,
               "-t", str(args.threads), "-v", str(args.min_overlap)]
    if args.debug:
        cmdline.append("-d")
    if args.read_type != "subassemblies":
        cmdline.append("-g")
    cmdline.extend([input_assembly, ",".join(args.reads),
                    out_folder, config_file])

    try:
        logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise RepeatException(str(e))
    except OSError as e:
        raise AssembleException(str(e))
