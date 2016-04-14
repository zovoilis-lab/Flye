from __future__ import print_function
import abruijn.fasta_parser as fp
import bisect

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
    profile = _compute_profile(alignment, genome_len)
    partition = _get_partition(profile)
    return _get_bubble_seqs(alignment, partition, genome_len)


def output_bubbles(bubbles, out_file):
    with open(out_file, "w") as f:
        for bubble in bubbles:
            if len(bubble.branches) == 0:
                print("Warrning: empty bubble {0}".format(bubble.bubble_id))
                continue

            consensus = sorted(bubble.branches,
                               key=len)[len(bubble.branches) / 2]
            f.write(">bubble {0} {1}\n".format(bubble.bubble_id,
                                               len(bubble.branches)))
            f.write(consensus + "\n")
            for branch_id, branch in enumerate(bubble.branches):
                f.write(">{0}\n".format(branch_id))
                f.write(branch + "\n")


def _is_solid_kmer(profile, position, kmer_length):
    MISSMATCH_RATE = 0.2
    INS_RATE = 0.2
    for i in xrange(position, position + kmer_length):
        local_missmatch = float(profile[i].num_missmatch +
                                profile[i].num_deletions) / profile[i].coverage
        local_ins = float(profile[i].num_inserts) / profile[i].coverage
        if local_missmatch > MISSMATCH_RATE or local_ins > INS_RATE:
            return False
    return True


def _is_gold_kmer(profile, position, kmer_length):
    for i in xrange(position, position + kmer_length - 1):
        if profile[i].nucl == profile[i + 1].nucl:
            return False
    return True


def _compute_profile(alignment, genome_len):
    print("Computing profile")
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
    print("Partitioning genome")
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
    #divided = False
    prev_part = -SOLID_LEN
    for prof_pos in xrange(0, len(profile) - GOLD_LEN):
        if solid_flags[prof_pos]:
            #if is_gold_kmer(profile, prof_pos, GOLD_LEN) and not divided:
            if (_is_gold_kmer(profile, prof_pos, GOLD_LEN) and
                prof_pos + GOLD_LEN / 2 - prev_part > SOLID_LEN):
                #divided = True
                prev_part = prof_pos + GOLD_LEN / 2
                partition.append(prof_pos + GOLD_LEN / 2)
        #else:
        #    divided = False

    return partition


def _get_bubble_seqs(alignment, partition, genome_len):
    print("Forming bubble sequences")
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
