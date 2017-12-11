#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Patches large structural variations in assembly by replacing
them with reads
"""

from collections import defaultdict, namedtuple
from itertools import izip
import math
import logging

import flye.fasta_parser as fp
import flye.config as config


DiscordandRead = namedtuple("DiscordandRead", ["alignment", "loop"])
AlignmentCluster = namedtuple("AlignmentCluster", ["reads", "start", "end"])
ReadPatch = namedtuple("ReadPatch", ["alignment", "start", "end"])

logger = logging.getLogger()


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
    MIN_LOOP = 500
    MIN_GAP = 100
    MAX_OVERHANG = 500
    MIN_CLUSTER = 10

    discordand_reads = []
    for aln in alignment:
        if (aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT or
            aln.qry_start > MAX_OVERHANG or
            aln.qry_len - aln.qry_end > MAX_OVERHANG):
            continue

        trg_gaps = 0
        qry_gaps = 0
        trg_strip = 0
        qry_strip = 0

        for trg_nuc, qry_nuc in izip(aln.trg_seq, aln.qry_seq):
            if trg_nuc == "-":
                trg_strip += 1
            else:
                if trg_strip > MIN_GAP:
                    trg_gaps += trg_strip
                trg_strip = 0

            if qry_nuc == "-":
                qry_strip += 1
            else:
                if qry_strip > MIN_GAP:
                    qry_gaps += qry_strip
                qry_strip = 0

        #minus = insertion, plus = deletion
        if abs(trg_gaps - qry_gaps) > MIN_LOOP:
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
    overlap_clusters = filter(lambda c: len(c) > MIN_CLUSTER, overlap_clusters)

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
            logger.debug("\t{0}\t{1}\t{2}\t{3}"
                         .format(aln.trg_start, aln.trg_end, aln.trg_sign, loop))

        cluster.reads.sort(key=lambda (a, l): l)
        chosen_patch = cluster.reads[len(cluster.reads) / 2]

        patches.append(ReadPatch(chosen_patch.alignment, cluster.start,
                                 cluster.end))
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
        out_sequence += sequence[prev_cut : patch.start]

        trg_seq = patch.alignment.trg_seq
        qry_seq = patch.alignment.qry_seq
        if patch.alignment.trg_sign == "-":
            trg_seq = fp.reverse_complement(trg_seq)
            qry_seq = fp.reverse_complement(qry_seq)

        seq_start = None
        seq_end = None
        nucl_count = 0

        for pos, nucl in enumerate(trg_seq):
            if nucl != "-":
                nucl_count += 1
            if nucl_count == patch.start - patch.alignment.trg_start:
                seq_start = pos
            if nucl_count == patch.end - patch.alignment.trg_start:
                seq_end = pos
                break

        out_sequence += qry_seq[seq_start : seq_end].replace("-", "")
        prev_cut = patch.end

    out_sequence += sequence[prev_cut:]
    return out_sequence
