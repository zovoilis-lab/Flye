#!/usr/bin/env python

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
            err_rate = 1 - float(tokens[17].count("|")) / len(tokens[17])
            alignments.append(Alignment(tokens[0], tokens[5], int(tokens[2]),
                                        int(tokens[3]), tokens[4],
                                        int(tokens[1]), int(tokens[7]),
                                        int(tokens[8]), tokens[9],
                                        int(tokens[6]), tokens[16],
                                        tokens[18], err_rate))
            errors.append(err_rate)

    mean_err = float(sum(errors)) / len(errors)
    return alignments, mean_err


def label_reads(alignment):
    renaming = {}
    alignment.sort(key=lambda a: a.trg_start)
    for read_id, aln in enumerate(alignment):
        chimeric = abs(float(aln.qry_end - aln.qry_start -
                             aln.qry_len)) / aln.qry_len > 0.2
        new_name = "id:{0}_len:{1}_chim:{2}_s:{3}_e:{4}" \
                        .format(read_id, aln.qry_len, int(chimeric),
                                aln.trg_start, aln.trg_end)
        renaming[aln.qry_id] = new_name

    return renaming


def rename(reads_in, rename_dict):
    unaln_id = 0
    for seq in SeqIO.parse(reads_in, "fasta"):
        print seq.id, type(seq.id)
        if seq.id in rename_dict:
            seq.id = rename_dict[seq.id]
        else:
            seq.id = "unaln_{0}".format(unaln_id)
            unaln_id += 1
        seq.description = ""
        SeqIO.write(seq, sys.stdout, "fasta")


def main():
    if len(sys.argv) != 3:
        print("Usage: label-reads.py read alignment")
        return 1

    alignment, err = _parse_blasr(sys.argv[2])
    rename_dict = label_reads(alignment)
    rename(sys.argv[1], rename_dict)


if __name__ == "__main__":
    main()
