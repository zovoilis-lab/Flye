#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Separates alignment into small bubbles for further correction
"""

import bisect
import logging
from collections import defaultdict, namedtuple

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
    def __init__(self, contig_id, bubble_id):
        self.contig_id = contig_id
        self.bubble_id = bubble_id
        self.branches = []
        self.consensus = ""


def get_bubbles(alignment, extended_simple):
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
        partition = _get_partition(profile, extended_simple)
        bubbles.extend(_get_bubble_seqs(ctg_aln, partition,
                                        contigs_info[ctg_id].length,
                                        ctg_id))
    bubbles = _add_consensus(bubbles)

    return bubbles



def output_bubbles(bubbles, out_file):
    """
    Outputs list of bubbles into file
    """
    with open(out_file, "w") as f:
        for bubble in bubbles:
            f.write(">{0} {1} {2}\n".format(bubble.contig_id,
                                            bubble.bubble_id,
                                            len(bubble.branches)))
            f.write(bubble.consensus + "\n")
            for branch_id, branch in enumerate(bubble.branches):
                f.write(">{0}\n".format(branch_id))
                f.write(branch + "\n")


def _add_consensus(bubbles):
    """
    Adds consensus sequences and filters outliers
    """
    new_bubbles = []
    for bubble in bubbles:
        if len(bubble.branches) == 0:
            logger.debug("Empty bubble {0}".format(bubble.bubble_id))
            continue

        consensus = sorted(bubble.branches,
                           key=len)[len(bubble.branches) / 2]
        bubble.consensus = consensus

        new_branches = []
        for branch in bubble.branches:
            incons_rate = float(abs(len(branch) -
                                    len(consensus))) / len(consensus)
            if incons_rate < 0.5:
                new_branches.append(branch)
            #else:
            #    logger.warning("Branch inconsistency with rate {0}, id {1}"
            #                    .format(incons_rate, bubble.bubble_id))

        bubble.branches = new_branches
        new_bubbles.append(bubble)

    return new_bubbles


def _is_solid_kmer(profile, position):
    """
    Checks if the kmer at given position is solid
    """
    SOLID_LEN = config.vals["solid_kmer_length"]
    MISSMATCH_RATE = config.vals["solid_missmatch"]
    INS_RATE = config.vals["solid_indel"]
    for i in xrange(position, position + SOLID_LEN):
        if profile[i].coverage == 0:
            return False
        local_missmatch = float(profile[i].num_missmatch +
                                profile[i].num_deletions) / profile[i].coverage
        local_ins = float(profile[i].num_inserts) / profile[i].coverage
        if local_missmatch > MISSMATCH_RATE or local_ins > INS_RATE:
            return False
    return True


def _is_simple_kmer(profile, position, extended_simple):
    """
    Checks if the kmer with center at the given position is simple
    """
    SIMPLE_LEN = config.vals["simple_kmer_length"]

    ext_simple = SIMPLE_LEN * 3
    nucl_str = map(lambda p: p.nucl, profile[position - ext_simple / 2 :
                                             position + ext_simple / 2])

    #single nucleotide homopolymers
    for i in xrange(ext_simple / 2 - SIMPLE_LEN / 2,
                    ext_simple / 2 + SIMPLE_LEN / 2 - 1):
        if nucl_str[i] == nucl_str[i + 1]:
            #logger.debug("single" + str(nucl_str))
            return False

    if extended_simple:
        #dinucleotide homopolymers
        for shift in [0, 1]:
            for i in xrange(SIMPLE_LEN - shift - 1):
                pos = ext_simple / 2 - SIMPLE_LEN + shift + i * 2
                if (nucl_str[pos : pos + 2] == nucl_str[pos + 2 : pos + 4]):
                    #logger.debug("di" + "".join(nucl_str))
                    return False

        #trinucleotide homopolymers
        for shift in [0, 1, 2]:
            for i in xrange(SIMPLE_LEN - shift - 1):
                pos = shift + i * 3
                if (nucl_str[pos : pos + 3] == nucl_str[pos + 3 : pos + 6]):
                    #logger.debug("tri" + "".join(nucl_str))
                    return False

    #logger.debug("".join(nucl_str))
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
        for i in xrange(len(trg_seq)):
            if trg_seq[i] == "-":
                trg_offset -= 1
            trg_pos = (aln.trg_start + i + trg_offset) % genome_len

            if trg_seq[i] == "-":
                profile[trg_pos].num_inserts += 1
            else:
                profile[trg_pos].nucl = trg_seq[i]
                if profile[trg_pos].nucl == "N" and qry_seq[i] != "-":
                    profile[trg_pos].nucl = qry_seq[i]

                profile[trg_pos].coverage += 1

                if qry_seq[i] == "-":
                    profile[trg_pos].num_deletions += 1
                elif trg_seq[i] != qry_seq[i]:
                    profile[trg_pos].num_missmatch += 1

    return profile


def _get_partition(profile, extended_simple):
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
    prof_pos = SOLID_LEN
    prev_part = SOLID_LEN
    while prof_pos < len(profile) - SOLID_LEN:
        landmark = (all(solid_flags[prof_pos : prof_pos + SIMPLE_LEN]) and
                    _is_simple_kmer(profile, prof_pos + SIMPLE_LEN / 2,
                                    extended_simple))

        if landmark or prof_pos - prev_part > MAX_BUBBLE:
            partition.append(prof_pos + SIMPLE_LEN / 2)
            prev_part = partition[-1]
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    logger.debug("Partitioned into {0} segments".format(len(partition) + 1))

    return partition


def _get_bubble_seqs(alignment, partition, genome_len, ctg_id):
    """
    Given genome landmarks, forms bubble sequences
    """
    logger.debug("Forming bubble sequences")
    MIN_ALIGNMENT = config.vals["min_alignment_length"]

    bubbles = [Bubble(ctg_id, x) for x in xrange(len(partition) + 1)]
    for aln in alignment:
        if aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT:
            continue

        if aln.trg_sign == "+":
            trg_seq, qry_seq = aln.trg_seq, aln.qry_seq
        else:
            trg_seq = fp.reverse_complement(aln.trg_seq)
            qry_seq = fp.reverse_complement(aln.qry_seq)

        trg_offset = 0
        prev_bubble_id = bisect.bisect(partition, aln.trg_start % genome_len)
        first_segment = True
        branch_start = None
        for i in xrange(len(trg_seq)):
            if trg_seq[i] == "-":
                trg_offset -= 1
                continue
            trg_pos = (aln.trg_start + i + trg_offset) % genome_len

            bubble_id = bisect.bisect(partition, trg_pos)
            if bubble_id != prev_bubble_id:
                if not first_segment:
                    branch_seq = qry_seq[branch_start:i].replace("-", "")
                    if len(branch_seq):
                        bubbles[prev_bubble_id].branches.append(branch_seq)

                first_segment = False
                prev_bubble_id = bubble_id
                branch_start = i

    return bubbles
