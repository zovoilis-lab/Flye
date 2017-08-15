#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Quick and dirty alignment consensus
"""

import logging
from collections import defaultdict
from itertools import izip
import multiprocessing
import signal

from abruijn.alignment import shift_gaps, parse_alignment
import abruijn.config as config

logger = logging.getLogger()

class Profile:
    def __init__(self):
        self.insertions = defaultdict(str)
        #self.deletions = 0
        self.matches = {"A" : 0, "C": 0, "G" : 0, "T" : 0, "N" : 0, "-" : 0}
        self.nucl = None


def _thread_worker(args_tuple):
    ctg_id, contigs_info, ctg_aln, results_queue = args_tuple
    #ctg_aln = parse_alignment(aln_path, ctg_id)

    profile = _contig_profile(ctg_aln, contigs_info[ctg_id].length)
    sequence = _flattern_profile(profile)
    results_queue.put((ctg_id, sequence))


def get_consensus(alignment_path, contigs_info, num_proc):
    """
    Main function
    """
    #making sure the main process catches SIGINT
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = multiprocessing.Pool(processes=num_proc)
    signal.signal(signal.SIGINT, orig_sigint)

    aln_by_ctg = defaultdict(list)
    for aln in parse_alignment(alignment_path):
        aln_by_ctg[aln.trg_id].append(aln)

    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    arguments = []
    for ctg_id in contigs_info:
        arguments.append((ctg_id, contigs_info,
                          aln_by_ctg[ctg_id], results_queue))

    try:
        res = pool.map_async(_thread_worker, arguments)
        res.get(99999999999)
    except KeyboardInterrupt:
        pool.terminate()
    else:
        pool.close()
    pool.join()

    out_fasta = {}
    while not results_queue.empty():
        ctg_id, ctg_seq = results_queue.get()
        out_fasta[ctg_id] = ctg_seq

    return out_fasta


def _contig_profile(alignment, genome_len):
    """
    Computes alignment profile
    """
    MIN_ALIGNMENT = config.vals["min_alignment_length"]

    profile = [Profile() for _ in xrange(genome_len)]
    for aln in alignment:
        if aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT:
            continue

        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in izip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            prof_elem = profile[trg_pos]
            if trg_nuc == "-":
                prof_elem.insertions[aln.qry_id] += qry_nuc
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    return profile


def _flattern_profile(profile):
    growing_seq = []
    ins_group = defaultdict(int)

    for elem in profile:
        ins_group.clear()
        for ins_str in elem.insertions.values():
            ins_group[ins_str] += 1

        coverage = sum(elem.matches.values())

        max_match = elem.nucl
        if elem.matches:
            max_match = max(elem.matches, key=elem.matches.get)
        max_insert = None
        if ins_group:
            max_insert = max(ins_group, key=ins_group.get)

        if max_match != "-":
            growing_seq += max_match
        if max_insert and ins_group[max_insert] > coverage / 2:
            growing_seq += max_insert

    return "".join(growing_seq)
