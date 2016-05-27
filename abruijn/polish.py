#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs polishing binary in parallel and concatentes output
"""

import logging
import random
import subprocess
import os
from collections import defaultdict
from threading import Thread

import abruijn.bubbles as bbl
import abruijn.fasta_parser as fp
from abruijn.utils import which
import abruijn.config as config


POLISH_BIN = "abruijn-polish"

logger = logging.getLogger()


class PolishException(Exception):
    pass


def check_binaries():
    if not which(POLISH_BIN):
        raise PolishException("polishing binary was not found. "
                              "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([POLISH_BIN, "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        raise PolishException("Some error inside native {0} module: {1}"
                              .format(POLISH_BIN, e))


def polish(bubbles, num_proc, work_dir, profile):
    logger.info("Polishing draft assembly")

    subs_matrix = os.path.join(os.environ["ABRUIJN_RES"],
                               config.vals["profiles"][profile]["subs_matrix"])
    hopo_matrix = os.path.join(os.environ["ABRUIJN_RES"],
                               config.vals["profiles"][profile]["hopo_matrix"])

    bubbles_file = os.path.join(work_dir, "bubbles.fasta")
    bbl.output_bubbles(bubbles, bubbles_file)
    consensus_out = os.path.join(work_dir, "consensus.fasta")

    _run_polish_bin(bubbles_file, subs_matrix, hopo_matrix,
                    consensus_out, num_proc)
    return _compose_sequence([consensus_out])


def _run_polish_bin(bubbles_in, subs_matrix, hopo_matrix,
                    consensus_out, num_threads):
    """
    Invokes polishing binary
    """
    cmdline = [POLISH_BIN, bubbles_in, subs_matrix,
               hopo_matrix, consensus_out, "-t", str(num_threads)]
    try:
        subprocess.check_call(cmdline, stderr=open(os.devnull, "w"))
    except (subprocess.CalledProcessError, OSError) as e:
        raise PolishException("Error while running polish binary: " + str(e))


def _compose_sequence(consensus_files):
    """
    Concatenates bubbles consensuses into genome
    """
    consensuses = defaultdict(list)
    for file_name in consensus_files:
        with open(file_name, "r") as f:
            header = True
            bubble_id = None
            for line in f:
                if header:
                    tokens = line.strip().split(" ")
                    ctg_id = tokens[0][1:]
                    ctg_pos = int(tokens[1])
                else:
                    consensuses[ctg_id].append((ctg_pos, line.strip()))
                header = not header

    polished_fasta = {}
    for ctg_id, seqs in consensuses.iteritems():
        sorted_seqs = map(lambda p: p[1], sorted(seqs, key=lambda p: p[0]))
        polished_fasta[ctg_id] = "".join(sorted_seqs)

    return polished_fasta
