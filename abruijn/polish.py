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

POLISH_BIN = "abruijn-polish"
SUBS_MATRIX = "pacbio_error.txt"
logger = logging.getLogger()


class PolishException(Exception):
    pass


def check_binaries():
    if not which(POLISH_BIN):
        raise PolishException("polish binary is not found. Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([POLISH_BIN, "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        raise PolishException("Some error inside native {0} module: {1}"
                              .format(POLISH_BIN, e))


def polish(bubbles, num_proc, work_dir):
    logger.info("Polishing draft assembly")
    buckets = [[] for _ in xrange(num_proc)]
    for i in xrange(len(bubbles)):
        buckets[random.randint(0, num_proc - 1)].append(bubbles[i])

    out_files = _run_parallel(buckets, work_dir)
    return _compose_sequence(out_files)


def _run_polish_bin(bubbles_in, subs_matrix, consensus_out, thread_error):
    """
    Invokes polishing binary
    """
    cmdline = [POLISH_BIN, bubbles_in, subs_matrix, consensus_out]
    try:
        subprocess.check_call(cmdline, stderr=open(os.devnull, "w"))
    except (subprocess.CalledProcessError, OSError) as e:
        print("Error while running polish binary: " + str(e))
        thread_error[0] = 1


def _compose_sequence(consensus_files):
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


def _run_parallel(buckets, work_dir):
    thread_error = [0]
    threads = []
    subs_matrix = os.path.join(os.environ["ABRUIJN_RES"], SUBS_MATRIX)
    output_files = []
    for i, bucket in enumerate(buckets):
        instance_in = os.path.join(work_dir, "bubbles_part_{0}.fasta".format(i))
        bbl.output_bubbles(bucket, instance_in)

        instance_out = os.path.join(work_dir, "consensus_{0}.fasta".format(i))
        output_files.append(instance_out)

        thread = Thread(target=_run_polish_bin, args=(instance_in, subs_matrix,
                                                      instance_out, thread_error))
        thread.start()
        threads.append(thread)

    for t in threads:
        t.join()

    if thread_error[0] != 0:
        raise PolishException("There were errors in polishing threads, exiting")

    return output_files
