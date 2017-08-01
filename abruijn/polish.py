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

def get_root():
    return os.path.dirname(__file__)


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


def polish(bubbles, num_proc, err_mode, work_dir, iter_id):
    _ROOT = get_root()

    subs_matrix = os.path.join(_ROOT, 'resource',
                               config.vals["err_modes"][err_mode]["subs_matrix"])
    hopo_matrix = os.path.join(_ROOT, 'resource',
                               config.vals["err_modes"][err_mode]["hopo_matrix"])

    bubbles_file = os.path.join(work_dir, "bubbles_{0}.fasta".format(iter_id))
    bbl.output_bubbles(bubbles, bubbles_file)
    consensus_out = os.path.join(work_dir, "consensus_{0}.fasta"
                                                             .format(iter_id))

    _run_polish_bin(bubbles_file, subs_matrix, hopo_matrix,
                    consensus_out, num_proc)
    return _compose_sequence([consensus_out])


def _run_polish_bin(bubbles_in, subs_matrix, hopo_matrix,
                    consensus_out, num_threads):
    """
    Invokes polishing binary
    """
    cmdline = [POLISH_BIN, "-t", str(num_threads), bubbles_in, subs_matrix,
               hopo_matrix, consensus_out]
    try:
        subprocess.check_call(cmdline)
    except (subprocess.CalledProcessError, OSError) as e:
        raise PolishException("Error while running polish binary: " + str(e))


def _compose_sequence(consensus_files):
    """
    Concatenates bubbles consensuses into genome
    """
    consensuses = defaultdict(list)
    coverage = defaultdict(list)
    for file_name in consensus_files:
        with open(file_name, "r") as f:
            header = True
            for line in f:
                if header:
                    tokens = line.strip().split(" ")
                    ctg_id = tokens[0][1:]
                    ctg_pos = int(tokens[1])
                    coverage[ctg_id].append(int(tokens[2]))
                else:
                    consensuses[ctg_id].append((ctg_pos, line.strip()))
                header = not header

    polished_fasta = {}
    for ctg_id, seqs in consensuses.iteritems():
        sorted_seqs = map(lambda p: p[1], sorted(seqs, key=lambda p: p[0]))
        concat_seq = "".join(sorted_seqs)
        mean_coverage = sum(coverage[ctg_id]) / len(coverage[ctg_id])
        extended_id = "{0}_len:{1}_cov:{2}".format(ctg_id, len(concat_seq),
                                                   mean_coverage)
        polished_fasta[extended_id] = concat_seq

    return polished_fasta
