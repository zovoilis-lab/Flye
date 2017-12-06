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
    except (subprocess.CalledProcessError, OSError) as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))


def polish(bubbles_file, num_proc, err_mode, work_dir, iter_id, out_polished):
    _ROOT = os.path.dirname(__file__)

    subs_matrix = os.path.join(_ROOT, "resource",
                               config.vals["err_modes"][err_mode]["subs_matrix"])
    hopo_matrix = os.path.join(_ROOT, "resource",
                               config.vals["err_modes"][err_mode]["hopo_matrix"])

    consensus_out = os.path.join(work_dir, "consensus_{0}.fasta"
                                                .format(iter_id))
    _run_polish_bin(bubbles_file, subs_matrix, hopo_matrix,
                    consensus_out, num_proc)
    polished_fasta, polished_lengths = _compose_sequence([consensus_out])
    fp.write_fasta_dict(polished_fasta, out_polished)

    return polished_lengths


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
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))


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
    polished_stats = {}
    for ctg_id, seqs in consensuses.iteritems():
        sorted_seqs = map(lambda p: p[1], sorted(seqs, key=lambda p: p[0]))
        concat_seq = "".join(sorted_seqs)
        mean_coverage = sum(coverage[ctg_id]) / len(coverage[ctg_id])
        polished_fasta[ctg_id] = concat_seq
        polished_stats[ctg_id] = len(concat_seq)

    return polished_fasta, polished_stats
