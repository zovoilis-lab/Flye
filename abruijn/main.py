from __future__ import print_function
import sys
import os

import abruijn.alignment as aln
import abruijn.bubbles as bbl
import abruijn.polish as pol
import abruijn.fasta_parser as fp

NUM_PROC = 8

def main():
    if len(sys.argv) != 3:
        print("Usage: abruijn.py pre_assembly reads_file")
        return 1

    workdir = "workdir"
    alignment, genome_len = aln.get_alignment(sys.argv[1], sys.argv[2],
                                              NUM_PROC, workdir)
    bubbles = bbl.get_bubbles(alignment, genome_len)
    polished_seq = pol.polish(bubbles, NUM_PROC, workdir)
    fp.write_fasta_dict({"genome": polished_seq}, "genome.fasta")
    return 0


if __name__ == "__main__":
    sys.exit(main())
