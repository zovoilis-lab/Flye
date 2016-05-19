#!/usr/bin/env python

import sys
from Bio import SeqIO

for e_num, seq in enumerate(SeqIO.parse(sys.argv[1], "fasta")):
    seq.id = "seq_{0}".format(e_num)
    seq.description = ""
    SeqIO.write(seq, sys.stdout, "fasta")
