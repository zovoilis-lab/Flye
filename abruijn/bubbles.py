#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Separates alignment into small bubbles for further correction
"""

import logging
from collections import defaultdict, namedtuple
from bisect import bisect
from itertools import izip
import math

import abruijn.fasta_parser as fp
import abruijn.config as config
from abruijn.alignment import shift_gaps


logger = logging.getLogger()


class ProfileInfo:
    def __init__(self):
        self.nucl = ""
        self.num_inserts = 0
        self.num_deletions = 0
        self.num_missmatch = 0
        self.coverage = 0


class Bubble:
    def __init__(self, contig_id, position):
        self.contig_id = contig_id
        self.position = position
        self.branches = []
        self.consensus = ""


def get_bubbles(alignment, contigs_info, err_mode):
    """
    The main function: takes an alignment and returns bubbles
    """
    aln_by_ctg = defaultdict(list)

    for aln in alignment:
        aln_by_ctg[aln.trg_id].append(aln)

    bubbles = []
    for ctg_id, ctg_aln in aln_by_ctg.iteritems():
        logger.debug("Processing {0}".format(ctg_id))
        profile = _compute_profile(ctg_aln, contigs_info[ctg_id].length)
        partition = _get_partition(profile, err_mode)
        bubbles.extend(_get_bubble_seqs(ctg_aln, profile, partition,
                                        contigs_info[ctg_id]))

    bubbles = _filter_outliers(bubbles)
    return bubbles


def output_bubbles(bubbles, out_file):
    """
    Outputs list of bubbles into file
    """
    with open(out_file, "w") as f:
        for bubble in bubbles:
            f.write(">{0} {1} {2}\n".format(bubble.contig_id,
                                            bubble.position,
                                            len(bubble.branches)))
            f.write(bubble.consensus + "\n")
            for branch_id, branch in enumerate(bubble.branches):
                f.write(">{0}\n".format(branch_id))
                f.write(branch + "\n")


def _filter_outliers(bubbles):
    new_bubbles = []
    for bubble in bubbles:
        if len(bubble.branches) == 0:
            logger.debug("Empty bubble {0}".format(bubble.position))
            continue

        new_branches = []
        median_branch = (sorted(bubble.branches, key=len)
                                [len(bubble.branches) / 2])
        if len(median_branch) == 0:
            continue

        for branch in bubble.branches:
            incons_rate = float(abs(len(branch) -
                                len(median_branch))) / len(median_branch)
            if incons_rate < 0.5:
                if len(branch) == 0:
                    branch = "A"
                    logger.debug("Zero branch")
                new_branches.append(branch)

        new_bubbles.append(Bubble(bubble.contig_id, bubble.position))
        new_bubbles[-1].consensus = bubble.consensus
        new_bubbles[-1].branches = new_branches

    return new_bubbles


def _is_solid_kmer(profile, position, err_mode):
    """
    Checks if the kmer at given position is solid
    """
    MISSMATCH_RATE = config.vals["err_modes"][err_mode]["solid_missmatch"]
    INS_RATE = config.vals["err_modes"][err_mode]["solid_indel"]
    SOLID_LEN = config.vals["solid_kmer_length"]

    for i in xrange(position, position + SOLID_LEN):
        if profile[i].coverage == 0:
            return False
        local_missmatch = float(profile[i].num_missmatch +
                                profile[i].num_deletions) / profile[i].coverage
        local_ins = float(profile[i].num_inserts) / profile[i].coverage
        if local_missmatch > MISSMATCH_RATE or local_ins > INS_RATE:
            return False
    return True


def _is_simple_kmer(profile, position):
    """
    Checks if the kmer with center at the given position is simple
    """
    SIMPLE_LEN = config.vals["simple_kmer_length"]

    extended_len = SIMPLE_LEN * 2
    nucl_str = map(lambda p: p.nucl, profile[position - extended_len / 2 :
                                             position + extended_len / 2])

    #single nucleotide homopolymers
    for i in xrange(extended_len / 2 - SIMPLE_LEN / 2,
                    extended_len / 2 + SIMPLE_LEN / 2 - 1):
        if nucl_str[i] == nucl_str[i + 1]:
            return False

    #dinucleotide homopolymers
    for shift in [0, 1]:
        for i in xrange(SIMPLE_LEN - shift - 1):
            pos = extended_len / 2 - SIMPLE_LEN + shift + i * 2
            if (nucl_str[pos : pos + 2] == nucl_str[pos + 2 : pos + 4]):
                return False

    """
    #trinucleotide homopolymers
    for shift in [0, 1, 2]:
        for i in xrange(SIMPLE_LEN - shift - 1):
            pos = shift + i * 3
            if (nucl_str[pos : pos + 3] == nucl_str[pos + 3 : pos + 6]):
                #logger.debug("tri" + "".join(nucl_str))
                return False
    """

    return True


def _compute_profile(alignment, genome_len):
    """
    Computes alignment profile
    """
    MIN_ALIGNMENT = config.vals["min_alignment_length"]
    logger.debug("Computing profile")

    profile = [ProfileInfo() for _ in xrange(genome_len)]
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
                prof_elem.num_inserts += 1
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.coverage += 1

                if qry_nuc == "-":
                    prof_elem.num_deletions += 1
                elif trg_nuc != qry_nuc:
                    prof_elem.num_missmatch += 1

            trg_pos += 1

    return profile


def _get_partition(profile, err_mode):
    """
    Partitions genome into sub-alignments at solid regions / simple kmers
    """
    logger.debug("Partitioning genome")
    SOLID_LEN = config.vals["solid_kmer_length"]
    SIMPLE_LEN = config.vals["simple_kmer_length"]
    MAX_BUBBLE = config.vals["max_bubble_length"]

    solid_flags = [False for _ in xrange(len(profile))]
    prof_pos = 0
    while prof_pos < len(profile) - SOLID_LEN:
        if _is_solid_kmer(profile, prof_pos, err_mode):
            for i in xrange(prof_pos, prof_pos + SOLID_LEN):
                solid_flags[i] = True
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    partition = []
    prev_partition = SOLID_LEN

    long_bubbles = 0
    prof_pos = SOLID_LEN
    while prof_pos < len(profile) - SOLID_LEN:
        cur_partition = prof_pos + SIMPLE_LEN / 2
        landmark = (all(solid_flags[prof_pos : prof_pos + SIMPLE_LEN]) and
                    _is_simple_kmer(profile, cur_partition))

        if prof_pos - prev_partition > MAX_BUBBLE:
            long_bubbles += 1

        if landmark or prof_pos - prev_partition > MAX_BUBBLE:
            partition.append(cur_partition)
            prev_partition = cur_partition
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    logger.debug("Partitioned into {0} segments".format(len(partition) + 1))
    logger.debug("Long bubbles: {0}".format(long_bubbles))

    return partition


def _get_bubble_seqs(alignment, profile, partition, contig_info):
    """
    Given genome landmarks, forms bubble sequences
    """
    logger.debug("Forming bubble sequences")
    MIN_ALIGNMENT = config.vals["min_alignment_length"]

    bubbles = []
    ext_partition = [0] + partition + [contig_info.length]
    for p_left, p_right in zip(ext_partition[:-1], ext_partition[1:]):
        bubbles.append(Bubble(contig_info.id, p_left))
        consensus = map(lambda p: p.nucl, profile[p_left : p_right])
        bubbles[-1].consensus = "".join(consensus)

    for aln in alignment:
        if aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT:
            continue

        bubble_id = bisect(partition, aln.trg_start % contig_info.length)
        next_bubble_start = ext_partition[bubble_id + 1]
        chromosome_start = (bubble_id == 0 and
                            not contig_info.type == "circular")
        chromosome_end = (aln.trg_end > partition[-1] and not
                          contig_info.type == "circular")

        branch_start = None
        first_segment = True
        trg_pos = aln.trg_start
        for i, trg_nuc in enumerate(aln.trg_seq):
            if trg_nuc == "-":
                continue
            if trg_pos >= contig_info.length:
                trg_pos -= contig_info.length

            if trg_pos >= next_bubble_start or trg_pos == 0:
                if not first_segment or chromosome_start:
                    branch_seq = aln.qry_seq[branch_start : i].replace("-", "")
                    bubbles[bubble_id].branches.append(branch_seq)

                first_segment = False
                bubble_id = bisect(partition, trg_pos)
                next_bubble_start = ext_partition[bubble_id + 1]
                branch_start = i

            trg_pos += 1

        if chromosome_end:
            branch_seq = aln.qry_seq[branch_start:].replace("-", "")
            bubbles[-1].branches.append(branch_seq)

    return bubbles
