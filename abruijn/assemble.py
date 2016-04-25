#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs assemble binary
"""

import subprocess
import logging
import os

from abruijn.utils import which

ASSEMBLE_BIN = "abruijn-assemble"
logger = logging.getLogger()


class AssembleException(Exception):
    pass


def check_binaries():
    if not which(ASSEMBLE_BIN):
        raise AssembleException("Assemble binary not found. "
                                "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([ASSEMBLE_BIN, "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        raise AssembleException("Some error inside native {0} module: {1}"
                                .format(ASSEMBLE_BIN, e))


def assemble(reads_file, out_file, kmer_size, min_cov):
    logger.info("Assembling reads")
    cmdline = [ASSEMBLE_BIN, reads_file, out_file,
               "-k", str(kmer_size), "-m", str(min_cov)]
    proc = subprocess.Popen(cmdline, stderr=subprocess.PIPE)
    for line in iter(proc.stderr.readline, ""):
        logger.debug(line.strip())
    ret_code = proc.wait()
    if ret_code:
        logger.error("Non-zero return code when calling {0} module: {1}"
                     .format(ASSEMBLE_BIN, ret_code))
        raise AssembleException("Error in assemble binary: " + e)
