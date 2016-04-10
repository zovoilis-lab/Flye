#!/usr/bin/env python2.7

from __future__ import print_function
import sys
import os
import subprocess
from collections import namedtuple

import wrapper.fasta_parser as fp


BLASR_BIN = "/home/fenderglass/Bioinf/tools/blasr/blasr.sh"
CIRCULAR_WINDOW = 50000
NUM_PROC = 8


Alignment = namedtuple("Alignment", ["qry_start", "qry_end", "qry_sign", "qry_len",
                                     "trg_start", "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate"])
class ProfileInfo:
    def __init__(self):
        self.num_inserts = 0
        self.num_deletions = 0
        self.num_missmatch = 0
        self.coverage = 0


def parse_blasr(filename):
    print("Parsing blasr")
    alignments = []
    errors = []
    with open(filename, "r") as f:
        for line in f:
            tokens = line.strip().split()
            err_rate = 1 - float(tokens[17].count("|")) / len(tokens[17])
            alignments.append(Alignment(int(tokens[2]), int(tokens[3]),
                                        tokens[4], int(tokens[1]),
                                        int(tokens[7]), int(tokens[8]),
                                        tokens[9], int(tokens[6]),
                                        tokens[16], tokens[18], err_rate))
            errors.append(err_rate)

    mean_err = float(sum(errors)) / len(errors)
    print("Read error rate: {0:5.2f}".format(mean_err))
    return alignments


def compose_raw_genome(contig_parts, out_file):
    """
    Concatenates contig parts and appends suffix from the beginning
    """
    fragment_index = {}
    genome_framents = fp.read_fasta_dict(contig_parts)
    for h, s in genome_framents.iteritems():
        idx = int(h.split("_")[1])
        fragment_index[idx] = s
    genome_seq = "".join(map(fragment_index.get, sorted(fragment_index)))
    genome_len = len(genome_seq)
    assert len(genome_seq) > CIRCULAR_WINDOW
    #genome_seq += genome_seq[:CIRCULAR_WINDOW]
    fp.write_fasta_dict({"contig_1" : genome_seq}, out_file)
    return genome_len


def run_blasr(reference_file, reads_file, num_proc, out_file):
    subprocess.check_call([BLASR_BIN, reads_file, reference_file, "-bestn", "1",
                           "-minMatch", "15", "-maxMatch", "25", "-m", "5",
                           "-nproc", str(num_proc), "-out", out_file])


def compute_profile(alignment, genome_len):
    print("Computing profile")
    MIN_ALIGNMENT = 5000
    profile = [ProfileInfo() for _ in xrange(genome_len)]
    num_accepted = 0
    for aln in alignment:
        #print(aln.err_rate, aln.trg_end - aln.trg_start)
        if aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT:
            continue

        num_accepted += 1

        if aln.trg_sign == "+":
            trg_seq, qry_seq = aln.trg_seq, aln.qry_seq
        else:
            trg_seq = fp.reverse_complement(aln.trg_seq)
            qry_seq = fp.reverse_complement(aln.qry_seq)
        trg_start = aln.trg_start % genome_len

        trg_offset = 0
        for i in xrange(len(trg_seq)):
            if trg_seq[i] == "-":
                trg_offset -= 1
            trg_pos = trg_start + i + trg_offset

            if trg_seq[i] == "-":
                profile[trg_pos].num_inserts += 1
            else:
                profile[trg_pos].coverage += 1
                if qry_seq[i] == "-":
                    profile[trg_pos].num_deletions += 1
                elif trg_seq[i] != qry_seq[i]:
                    profile[trg_pos].num_missmatch += 1

    print(num_accepted)
    return profile


def get_solid_regions(profile):
    for i in xrange(len(profile)):
        print(i, profile[i].coverage, profile[i].num_inserts,
              profile[i].num_deletions, profile[i].num_missmatch)


def main():
    if len(sys.argv) != 3:
        print("Usage: abruijn.py pre_assembly reads_file")
        return 1

    RAW_GENOME = "raw_genome.fasta"
    BLASR_AILNMENT = "blasr.m5"
    genome_len = compose_raw_genome(sys.argv[1], RAW_GENOME)
    #run_blasr(RAW_GENOME, sys.argv[2], NUM_PROC, BLASR_AILNMENT)
    alignment = parse_blasr(BLASR_AILNMENT)
    profile = compute_profile(alignment, genome_len)
    get_solid_regions(profile)

    return 0


if __name__ == "__main__":
    sys.exit(main())
