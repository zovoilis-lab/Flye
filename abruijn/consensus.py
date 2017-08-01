#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Quick and dirty alignment consensus
"""

import logging
from collections import defaultdict
from itertools import izip

from abruijn.alignment import shift_gaps
import abruijn.config as config

logger = logging.getLogger()

class Profile:
    def __init__(self):
        self.insertions = defaultdict(str)
        #self.deletions = 0
        self.matches = {"A" : 0, "C": 0, "G" : 0, "T" : 0, "N" : 0, "-" : 0}
        self.nucl = None


def get_consensus(alignment, contigs_info):
    out_fasta = {}
    aln_by_ctg = defaultdict(list)
    for aln in alignment:
        aln_by_ctg[aln.trg_id].append(aln)

    for ctg_id, ctg_aln in aln_by_ctg.iteritems():
        logger.debug("Processing {0}".format(ctg_id))
        profile = _contig_profile(ctg_aln, contigs_info[ctg_id].length)
        sequence = _flattern_profile(profile)
        out_fasta[ctg_id] = sequence

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
