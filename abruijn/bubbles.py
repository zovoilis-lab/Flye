#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Separates alignment into small bubbles for further correction
"""

import bisect
import logging
from collections import defaultdict, namedtuple
from copy import deepcopy
from itertools import izip
import math

import abruijn.fasta_parser as fp
import abruijn.config as config

ContigInfo = namedtuple("ContigInfo", ["id", "length", "type"])
DiscordandRead = namedtuple("DiscordandRead", ["alignment", "loop"])
AlignmentCluster = namedtuple("AlignmentCluster", ["reads", "start", "end"])

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

    #bubbles = _filter_outliers(bubbles)
    return bubbles


def patch_genome(alignment, reference_file, out_patched):
    aln_by_ctg = defaultdict(list)
    for aln in alignment:
        aln_by_ctg[aln.trg_id].append(aln)

    ref_fasta = fp.read_fasta_dict(reference_file)
    fixed_fasta = {}

    for ctg_id, ctg_aln in aln_by_ctg.iteritems():
        patches = _get_patching_alignmemnts(ctg_aln)
        fixed_sequence = _apply_patches(patches, ref_fasta[ctg_id])
        fixed_fasta[ctg_id] = fixed_sequence

    fp.write_fasta_dict(fixed_fasta, out_patched)


def _get_patching_alignmemnts(alignment):
    MIN_ALIGNMENT = config.vals["min_alignment_length"]

    discordand_reads = []
    for aln in alignment:
        if (aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT or
            aln.qry_start > 500 or aln.qry_len - aln.qry_end > 500):
            continue

        trg_gaps = 0
        qry_gaps = 0
        trg_strip = 0
        qry_strip = 0

        for trg_nuc, qry_nuc in izip(aln.trg_seq, aln.qry_seq):
            if trg_nuc == "-":
                trg_strip += 1
            else:
                if trg_strip > 100:
                    trg_gaps += trg_strip
                trg_strip = 0

            if qry_nuc == "-":
                qry_strip += 1
            else:
                if qry_strip > 100:
                    qry_gaps += qry_strip
                qry_strip = 0

        if abs(trg_gaps - qry_gaps) > 500:
            discordand_reads.append(DiscordandRead(aln, trg_gaps - qry_gaps))


    #create clusters of overlapping discordand reads
    overlap_clusters = []
    discordand_reads.sort(key = lambda t: t.alignment.trg_start)
    right_end = 0
    for read in discordand_reads:
        if not overlap_clusters or right_end < read.alignment.trg_start:
            overlap_clusters.append([])
        right_end = max(read.alignment.trg_end, right_end)
        overlap_clusters[-1].append(read)
    ##

    #filter outliers
    overlap_clusters = filter(lambda c: len(c) > 10, overlap_clusters)

    ##split clusters based on loop rate and compute common overlap
    nonzero_clusters = []
    split_clusters = sum([_split_reads(cl) for cl in overlap_clusters], [])
    for cl in split_clusters:
        left_ends = []
        right_ends = []
        for read in cl:
            logger.debug("\t{0}".format(read.loop))
            left_ends.append(read.alignment.trg_start)
            right_ends.append(read.alignment.trg_end)

        common_left = sorted(left_ends)[3 * len(left_ends) / 4]
        common_right = sorted(right_ends)[len(right_ends) / 4]
        logger.debug("-----{0} {1}".format(common_left, common_right))
        if common_right > common_left:
            nonzero_clusters.append(AlignmentCluster(cl, common_left,
                                                     common_right))
    ###

    ##TODO: filter conflicting clusters
    non_conflicting_clusters = []
    for cl in nonzero_clusters:
        if (not non_conflicting_clusters or
            cl.start > non_conflicting_clusters[-1].end):

            non_conflicting_clusters.append(cl)
        elif len(cl.reads) > len(non_conflicting_clusters[-1].reads):
            non_conflicting_clusters[-1] = cl

    patches = []
    for cluster in non_conflicting_clusters:
        logger.debug("Cluster {0} {1}".format(cluster.start, cluster.end))
        for aln, loop in cluster.reads:
            logger.debug("\t{0}\t{1}\t{2}"
                         .format(aln.trg_start, aln.trg_end, loop))

        #num_positive = sum(1 if d > 0 else 0
        #                   for (a, d) in cluster.reads)
        #cluster_positive = num_positive > len(cluster.reads) / 2
        #filtered_alignments = [r for r in cluster.reads
        #                       if (r.loop > 0) == cluster_positive]
        #chosen_patch = max(filtered_alignments, key=lambda r: abs(r.loop))

        cluster.reads.sort(key=lambda (a, l): l)
        chosen_patch = cluster.reads[len(cluster.reads) / 2]

        patches.append(chosen_patch.alignment)
        logger.debug("Chosen: {0}".format(chosen_patch.loop))

    return patches


def mean(lst):
    return (float(sum(lst)) / len(lst))


def stddev(lst):
    if len(lst) < 2:
        return 0.0

    l_sum = 0
    l_mean = mean(lst)
    for val in lst:
        l_sum += float(val - l_mean) * float(val - l_mean)
    return math.sqrt(l_sum / (len(lst) - 1))


def _split_reads(reads):
    new_clusters = [reads]
    while True:
        any_split = False
        for cl in new_clusters:
            res = _split_two(cl)

            if len(res) > 1:
                new_clusters.remove(cl)
                new_clusters.extend(res)
                any_split = True
                break

        if not any_split:
            break

    return new_clusters


def _split_two(cluster_reads):
    MIN_CLUSTER = 5
    if len(cluster_reads) < 2 * MIN_CLUSTER:
        return [cluster_reads]

    best_diff = 0
    best_point = None
    for split_point in xrange(MIN_CLUSTER, len(cluster_reads) - MIN_CLUSTER):
        left_cluster = cluster_reads[:split_point]
        right_cluster = cluster_reads[split_point:]
        left_std = stddev(map(lambda r: r.loop, left_cluster))
        right_std = stddev(map(lambda r: r.loop, right_cluster))

        if best_point is None or left_std + right_std < best_diff:
            best_diff = left_std + right_std
            best_point = split_point

    left_cluster = cluster_reads[:best_point]
    right_cluster = cluster_reads[best_point:]
    left_std = stddev(map(lambda r: r.loop, left_cluster))
    right_std = stddev(map(lambda r: r.loop, right_cluster))
    left_mean = mean(map(lambda r: r.loop, left_cluster))
    right_mean = mean(map(lambda r: r.loop, right_cluster))

    #TODO: replace with t-test
    need_split = abs(left_mean - right_mean) > 2 * left_std

    if not need_split:
        return [cluster_reads]

    #check if SV coordinates are overlapping
    return [left_cluster, right_cluster]


def _apply_patches(patches, sequence):
    prev_cut = 0
    out_sequence = ""
    for patch in patches:
        out_sequence += sequence[prev_cut : patch.trg_start]
        patched_sequence = patch.qry_seq
        if patch.trg_sign == "-":
            patched_sequence = fp.reverse_complement(patched_sequence)

        out_sequence += patched_sequence.replace("-", "")
        prev_cut = patch.trg_end

    out_sequence += sequence[prev_cut:]
    return out_sequence


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
        for branch in bubble.branches:
            incons_rate = float(abs(len(branch) -
                                len(bubble.consensus))) / len(bubble.consensus)
            if len(branch) == 0:
                branch = "A"
            if incons_rate < 0.5:
                new_branches.append(branch)
            #else:
            #    logger.warning("Branch inconsistency with rate {0}, id {1}"
            #                    .format(incons_rate, bubble.position))

        new_bubbles.append(deepcopy(bubble))
        new_bubbles[-1].branches = new_branches

    return new_bubbles


def _is_solid_kmer(profile, position, kmer_length):
    """
    Checks if the kmer at given position is solid
    """
    MISSMATCH_RATE = config.vals["solid_missmatch"]
    INS_RATE = config.vals["solid_indel"]
    for i in xrange(position, position + kmer_length):
        if profile[i].coverage == 0:
            return False
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
    if (profile[position].nucl == profile[position + 2].nucl and
        profile[position + 1].nucl == profile[position +3].nucl):
        return False
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
        if _is_solid_kmer(profile, prof_pos, SOLID_LEN):
            for i in xrange(prof_pos, prof_pos + SOLID_LEN):
                solid_flags[i] = True
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    partition = []
    prev_part = 0
    for prof_pos in xrange(0, len(profile) - SIMPLE_LEN):
        if all(solid_flags[prof_pos : prof_pos + SIMPLE_LEN]):
            if (_is_simple_kmer(profile, prof_pos, SIMPLE_LEN) and
                    prof_pos + SIMPLE_LEN / 2 - prev_part > SOLID_LEN):

                if prof_pos + SIMPLE_LEN / 2 - prev_part > 1000:
                    logger.info("Long bubble {0}, at {1}"
                            .format(prof_pos + SIMPLE_LEN / 2 - prev_part,
                                    prof_pos + SIMPLE_LEN / 2))

                prev_part = prof_pos + SIMPLE_LEN / 2
                partition.append(prof_pos + SIMPLE_LEN / 2)

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
                    #if len(branch_seq):
                    bubbles[prev_bubble_id].branches.append(branch_seq)

                first_segment = False
                prev_bubble_id = bubble_id
                branch_start = i

    return bubbles
