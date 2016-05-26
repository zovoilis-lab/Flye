#!/usr/bin/env python

from __future__ import print_function
import sys
from collections import namedtuple

from Bio import SeqIO

Alignment = namedtuple("Alignment", ["qry_id", "trg_id", "qry_start", "qry_end",
                                     "qry_sign", "qry_len", "trg_start",
                                     "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate"])

def _parse_blasr(filename):
    """
    Parse Blasr output
    """
    alignments = []
    errors = []
    with open(filename, "r") as f:
        for line in f:
            tokens = line.strip().split()
            read_name = tokens[0].rsplit("/", 1)[0]    #Blasr, whyyy??
            err_rate = 1 - float(tokens[17].count("|")) / len(tokens[17])
            alignments.append(Alignment(read_name, tokens[5], int(tokens[2]),
                                        int(tokens[3]), tokens[4],
                                        int(tokens[1]), int(tokens[7]),
                                        int(tokens[8]), tokens[9],
                                        int(tokens[6]), tokens[16],
                                        tokens[18], err_rate))
            errors.append(err_rate)

    mean_err = float(sum(errors)) / len(errors)
    return alignments, mean_err


def _compute_profile(alignment):
    """
    Computes alignment profile
    """
    MIN_ALIGNMENT = 1000

    genome_len = alignment[0].trg_len
    coverage = [0 for _ in xrange(genome_len)]
    for aln in alignment:
        if aln.err_rate > 0.5 or aln.trg_end - aln.trg_start < MIN_ALIGNMENT:
            continue

        trg_seq, qry_seq = aln.trg_seq, aln.qry_seq
        trg_offset = 0
        for i in xrange(len(trg_seq)):
            if trg_seq[i] == "-":
                trg_offset -= 1
            trg_pos = (aln.trg_start + i + trg_offset) % genome_len

            if trg_seq[i] == "-":
                pass
            else:

                coverage[trg_pos] += 1

    mean_coverage = sum(coverage) / len(coverage)
    return coverage, mean_coverage


def label_reads(alignment, coverage):
    renaming = {}
    alignment.sort(key=lambda a: a.trg_start)
    for read_id, aln in enumerate(alignment):
        chimeric = abs(float(aln.qry_end - aln.qry_start -
                             aln.qry_len)) / aln.qry_len > 0.2
        read_cov = sum(coverage[aln.trg_start:aln.trg_end]) / (aln.trg_end -
                                                               aln.trg_start)
        new_name = "id:{0}_len:{1}_chim:{2}_cov:{3}_s:{4}_e:{5}" \
                        .format(read_id, aln.qry_len, int(chimeric),
                                read_cov, aln.trg_start, aln.trg_end)
        renaming[aln.qry_id] = new_name

    return renaming


def rename(reads_in, rename_dict, out_file):
    unaln_id = 0
    with open(out_file, "w") as f:
        for seq in SeqIO.parse(reads_in, "fasta"):
            if seq.id in rename_dict:
                seq.id = rename_dict[seq.id]
            else:
                seq.id = "unaln_{0}".format(unaln_id)
                unaln_id += 1
            seq.description = ""
            SeqIO.write(seq, f, "fasta")


def main():
    if len(sys.argv) != 4:
        print("Usage: label-reads.py reads alignment fasta_out")
        return 1

    alignment, mean_err = _parse_blasr(sys.argv[2])
    cov, mean_cov = _compute_profile(alignment)

    print("Mean alignment error", mean_err)
    print("Mean coverage", mean_cov)

    rename_dict = label_reads(alignment, cov)
    rename(sys.argv[1], rename_dict, sys.argv[3])


if __name__ == "__main__":
    main()
