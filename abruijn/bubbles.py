#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Separates alignment into small bubbles for further correction
"""

import logging
from collections import defaultdict, namedtuple
from bisect import bisect
from copy import deepcopy
from itertools import izip
import math

import abruijn.fasta_parser as fp
import abruijn.config as config

ContigInfo = namedtuple("ContigInfo", ["id", "length", "type"])

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


def get_bubbles(alignment):
    """
    The main function: takes an alignment and returns bubbles
    """
    logger.info("Separating draft genome into bubbles")
    aln_by_ctg = defaultdict(list)
    contigs_info = {}
    circular_window = config.vals["circular_window"]

    for aln in alignment:
        aln_by_ctg[aln.trg_id].append(aln)
        if aln.trg_id not in contigs_info:
            contig_type = aln.trg_id.split("_")[0]
            contig_len = aln.trg_len
            if contig_type == "circular" and contig_len > circular_window:
                contig_len -= circular_window
            contigs_info[aln.trg_id] = ContigInfo(aln.trg_id, contig_len,
                                                  contig_type)

    bubbles = []
    for ctg_id, ctg_aln in aln_by_ctg.iteritems():
        logger.debug("Processing {0}".format(ctg_id))
        profile = _compute_profile(ctg_aln, contigs_info[ctg_id].length)
        partition = _get_partition(profile)
        bubbles.extend(_get_bubble_seqs(ctg_aln, profile, partition,
                                        contigs_info[ctg_id].length,
                                        ctg_id))

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
            #else:
            #    logger.warning("Branch inconsistency with rate {0}, id {1}"
            #                    .format(incons_rate, bubble.position))

        new_bubbles.append(deepcopy(bubble))
        new_bubbles[-1].branches = new_branches

    return new_bubbles


def _is_solid_kmer(profile, position):
    """
    Checks if the kmer at given position is solid
    """
    MISSMATCH_RATE = config.vals["solid_missmatch"]
    INS_RATE = config.vals["solid_indel"]
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
            #logger.debug("single" + str(nucl_str))
            return False

    #dinucleotide homopolymers
    for shift in [0, 1]:
        for i in xrange(SIMPLE_LEN - shift - 1):
            pos = extended_len / 2 - SIMPLE_LEN + shift + i * 2
            if (nucl_str[pos : pos + 2] == nucl_str[pos + 2 : pos + 4]):
                #logger.debug("di" + "".join(nucl_str))
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



def _shift_gaps(seq_trg, seq_qry):
    """
    Shifts all ambigious query gaps to the right
    """
    lst_trg, lst_qry = list("$" + seq_trg + "$"), list("$" + seq_qry + "$")
    is_gap = False
    gap_start = 0
    for i in xrange(len(lst_trg)):
        if is_gap and lst_qry[i] != "-":
            is_gap = False
            swap_left = gap_start - 1
            swap_right = i - 1

            while (swap_left > 0 and swap_right >= gap_start and
                   lst_qry[swap_left] == lst_trg[swap_right]):
                lst_qry[swap_left], lst_qry[swap_right] = \
                            lst_qry[swap_right], lst_qry[swap_left]
                swap_left -= 1
                swap_right -= 1

        if not is_gap and lst_qry[i] == "-":
            is_gap = True
            gap_start = i

    return "".join(lst_qry[1 : -1])


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

        if aln.trg_sign == "+":
            trg_seq, qry_seq = aln.trg_seq, aln.qry_seq
        else:
            trg_seq = fp.reverse_complement(aln.trg_seq)
            qry_seq = fp.reverse_complement(aln.qry_seq)
        qry_seq = _shift_gaps(trg_seq, qry_seq)

        trg_offset = 0
        for i, (trg_nuc, qry_nuc) in enumerate(izip(trg_seq, qry_seq)):
            if trg_nuc == "-":
                trg_offset -= 1
            trg_pos = (aln.trg_start + i + trg_offset) % genome_len
            prof_elem = profile[trg_pos]

            if trg_nuc == "-":
                prof_elem.num_inserts += 1
            else:
                prof_elem.nucl = trg_nuc
                if prof_elem.nucl == "N" and qry_nuc != "-":
                    prof_elem.nucl = qry_nuc

                prof_elem.coverage += 1

                if qry_nuc == "-":
                    prof_elem.num_deletions += 1
                elif trg_nuc != qry_nuc:
                    prof_elem.num_missmatch += 1

    return profile


def _get_partition(profile):
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
        if _is_solid_kmer(profile, prof_pos):
            for i in xrange(prof_pos, prof_pos + SOLID_LEN):
                solid_flags[i] = True
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    partition = []
    prev_partition = SOLID_LEN

    prof_pos = SOLID_LEN
    while prof_pos < len(profile) - SOLID_LEN:
        cur_partition = prof_pos + SIMPLE_LEN / 2
        landmark = (all(solid_flags[prof_pos : prof_pos + SIMPLE_LEN]) and
                    _is_simple_kmer(profile, cur_partition))

        #if landmark or prof_pos - prev_part > MAX_BUBBLE:
        if landmark:
            partition.append(cur_partition)
            prev_partition = cur_partition
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    logger.debug("Partitioned into {0} segments".format(len(partition) + 1))

    return partition


def _get_bubble_seqs(alignment, profile, partition, genome_len, ctg_id):
    """
    Given genome landmarks, forms bubble sequences
    """
    logger.debug("Forming bubble sequences")
    MIN_ALIGNMENT = config.vals["min_alignment_length"]

    bubbles = []
    ext_partition = [0] + partition + [genome_len]
    for p_left, p_right in zip(ext_partition[:-1], ext_partition[1:]):
        bubbles.append(Bubble(ctg_id, p_left))
        consensus = map(lambda p: p.nucl, profile[p_left : p_right])
        bubbles[-1].consensus = "".join(consensus)

    for aln in alignment:
        if aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT:
            continue

        if aln.trg_sign == "+":
            trg_seq, qry_seq = aln.trg_seq, aln.qry_seq
        else:
            trg_seq = fp.reverse_complement(aln.trg_seq)
            qry_seq = fp.reverse_complement(aln.qry_seq)

        trg_offset = 0
        prev_bubble_id = bisect(partition, aln.trg_start % genome_len)
        first_segment = True
        branch_start = None
        for i, trg_nuc in enumerate(trg_seq):
            if trg_nuc == "-":
                trg_offset -= 1
                continue
            trg_pos = (aln.trg_start + i + trg_offset) % genome_len

            bubble_id = bisect(partition, trg_pos)
            if bubble_id != prev_bubble_id:
                if not first_segment:
                    branch_seq = qry_seq[branch_start:i].replace("-", "")
                    bubbles[prev_bubble_id].branches.append(branch_seq)

                first_segment = False
                prev_bubble_id = bubble_id
                branch_start = i

    return bubbles
