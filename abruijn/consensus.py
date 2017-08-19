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

from abruijn.alignment import shift_gaps, SynchronizedReader
import abruijn.config as config

logger = logging.getLogger()

class Profile:
    __slots__ = ("insertions", "matches", "nucl")

    def __init__(self):
        self.insertions = defaultdict(str)
        self.matches = defaultdict(int)
        self.nucl = "-"
        #self.length = length
        #self.insertions = [defaultdict(str) for _ in xrange(length)]
        #self.matches = [defaultdict(int) for _ in xrange(length)]
        #self.nucl = ["-" for _ in xrange(length)]


def _thread_worker(blasr_reader, contigs_info, results_queue,
                   error_queue):
    try:
        blasr_reader.init_reading()

        while not blasr_reader.is_eof():
            ctg_id, ctg_aln = blasr_reader.get_chunk()
            if ctg_id is None:
                break

            profile = _contig_profile(ctg_aln, contigs_info[ctg_id].length)
            sequence = _flattern_profile(profile)
            results_queue.put((ctg_id, sequence))

    except Exception as e:
        error_queue.put(e)


def get_consensus(alignment_path, contigs_info, num_proc):
    """
    Main function
    """
    blasr_reader = SynchronizedReader(alignment_path)
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    error_queue = manager.Queue()

    #making sure the main process catches SIGINT
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    threads = []
    for _ in xrange(num_proc):
        threads.append(multiprocessing.Process(target=_thread_worker,
                                               args=(blasr_reader, contigs_info,
                                                     results_queue, error_queue)))
    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()

    if not error_queue.empty():
        raise error_queue.get()

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

    #profile = Profile(genome_len)
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
                #profile.insertions[trg_pos][aln.qry_id] += qry_nuc
                prof_elem.insertions[aln.qry_id] += qry_nuc
            else:
                #profile.nucl[trg_pos] = trg_nuc
                prof_elem.nucl = trg_nuc
                #profile.matches[trg_pos][qry_nuc] += 1
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    return profile


def _flattern_profile(profile):
    growing_seq = []
    ins_group = defaultdict(int)

    #for i in xrange(profile.length):
    for elem in profile:
        pos_matches = elem.matches
        pos_insertions = elem.insertions
        pos_nucl = elem.nucl

        ins_group.clear()
        for ins_str in pos_insertions.values():
            ins_group[ins_str] += 1

        coverage = sum(pos_matches.values())

        max_match = pos_nucl
        if len(pos_matches):
            max_match = max(pos_matches, key=pos_matches.get)
        max_insert = None
        if ins_group:
            max_insert = max(ins_group, key=ins_group.get)

        if max_match != "-":
            growing_seq.append(max_match)
        if max_insert and ins_group[max_insert] > coverage / 2:
            growing_seq.append(max_insert)

    return "".join(growing_seq)
