#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Separates alignment into small bubbles for further correction
"""

import bisect
import logging

import abruijn.fasta_parser as fp

logger = logging.getLogger()

class ProfileInfo:
    def __init__(self):
        self.nucl = ""
        self.num_inserts = 0
        self.num_deletions = 0
        self.num_missmatch = 0
        self.coverage = 0


class Bubble:
    def __init__(self, bubble_id):
        self.bubble_id = bubble_id
        self.branches = []


def get_bubbles(alignment, genome_len):
    """
    The main function: takes an alignment and returns bubbles
    """
    logger.info("Separating draft genome into bubbles")
    profile = _compute_profile(alignment, genome_len)
    partition = _get_partition(profile)
    return _get_bubble_seqs(alignment, partition, genome_len)


def output_bubbles(bubbles, out_file):
    """
    Outputs list of bubbles into file
    """
    with open(out_file, "w") as f:
        for bubble in bubbles:
            if len(bubble.branches) == 0:
                logger.warning("Empty bubble {0}".format(bubble.bubble_id))
                continue

            consensus = sorted(bubble.branches,
                               key=len)[len(bubble.branches) / 2]
            f.write(">bubble {0} {1}\n".format(bubble.bubble_id,
                                               len(bubble.branches)))
            f.write(consensus + "\n")
            for branch_id, branch in enumerate(bubble.branches):
                #rate = float(abs(len(branch) - len(consensus))) / len(consensus)
                #if rate > 0.5:
                #    logger.warning("Branch inconsistency with rate {0}, id {1}"
                #                    .format(rate, bubble.bubble_id))
                f.write(">{0}\n".format(branch_id))
                f.write(branch + "\n")


def _is_solid_kmer(profile, position, kmer_length):
    """
    Checks if the kmer at given position is solid
    """
    MISSMATCH_RATE = 0.2
    INS_RATE = 0.2
    for i in xrange(position, position + kmer_length):
        local_missmatch = float(profile[i].num_missmatch +
                                profile[i].num_deletions) / profile[i].coverage
        local_ins = float(profile[i].num_inserts) / profile[i].coverage
        if local_missmatch > MISSMATCH_RATE or local_ins > INS_RATE:
            return False
    return True


def _is_simple_kmer(profile, position, kmer_length):
    """
    Checks if the kmer at given position is simple
    """
    for i in xrange(position, position + kmer_length - 1):
        if profile[i].nucl == profile[i + 1].nucl:
            return False
    return True


def _compute_profile(alignment, genome_len):
    """
    Computes alignment profile
    """
    logger.debug("Computing profile")
    MIN_ALIGNMENT = 5000
    profile = [ProfileInfo() for _ in xrange(genome_len)]
    for aln in alignment:
        if aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT:
            continue

        if aln.trg_sign == "+":
            trg_seq, qry_seq = aln.trg_seq, aln.qry_seq
        else:
            trg_seq = fp.reverse_complement(aln.trg_seq)
            qry_seq = fp.reverse_complement(aln.qry_seq)

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


def _get_partition(profile):
    """
    Partitions genome into sub-alignments at solid regions / simple kmers
    """
    logger.debug("Partitioning genome")
    SOLID_LEN = 10
    GOLD_LEN = 4
    solid_flags = [False for _ in xrange(len(profile))]
    prof_pos = 0
    num_solid = 0
    while prof_pos < len(profile) - SOLID_LEN:
        if _is_solid_kmer(profile, prof_pos, SOLID_LEN):
            for i in xrange(prof_pos, prof_pos + SOLID_LEN):
                solid_flags[i] = True
            prof_pos += SOLID_LEN
            num_solid += 1
        else:
            prof_pos += 1

    partition = []
    prev_part = 0
    for prof_pos in xrange(0, len(profile) - GOLD_LEN):
        if solid_flags[prof_pos]:
            if (_is_simple_kmer(profile, prof_pos, GOLD_LEN) and
                prof_pos + GOLD_LEN / 2 - prev_part > SOLID_LEN):
                prev_part = prof_pos + GOLD_LEN / 2
                partition.append(prof_pos + GOLD_LEN / 2)
    logger.debug("Partitioned into {0} segments".format(len(partition) + 1))

    return partition


def _get_bubble_seqs(alignment, partition, genome_len):
    """
    Given genome landmarks, forms bubble sequences
    """
    logger.debug("Forming bubble sequences")
    MIN_ALIGNMENT = 5000
    bubbles = [Bubble(x) for x in xrange(len(partition) + 1)]
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
                    #if prev_bubble_id == 0:
                    #    logger.debug("id {0}, b {1}, e {2}, seq_len {3}, rb {4}"
                    #                    .format(prev_bubble_id, branch_start,
                    #                            i, len(branch_seq),
                    #                            aln.trg_start + i + trg_offset))
                    if len(branch_seq):
                        bubbles[prev_bubble_id].branches.append(branch_seq)

                first_segment = False
                prev_bubble_id = bubble_id
                branch_start = i

    return bubbles
