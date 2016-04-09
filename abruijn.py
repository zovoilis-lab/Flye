#!/usr/bin/env python2.7

from __future__ import print_function
import sys
import os
import subprocess

import wrapper.fasta_parser as fp


BLASR_BIN = "/home/fenderglass/Bioinf/tools/blasr/blasr.sh"
CIRCULAR_WINDOW = 50000
NUM_PROC = 8


def compose_raw_genome(contig_parts, out_file):
    """
    Concatenates contig parts and appends suffix from the beginning
    """
    genome_framents = fp.read_fasta_dict(contig_parts)
    genome_seq = "".join(genome_framents.values())
    assert len(genome_seq) > CIRCULAR_WINDOW
    genome_seq += genome_seq[:CIRCULAR_WINDOW]
    fp.write_fasta_dict({"contig_1" : genome_seq}, out_file)


def run_blasr(reference_file, reads_file, num_proc, out_file):
    subprocess.check_call([BLASR_BIN, reads_file, reference_file, "-bestn", "1",
                           "-minMatch", "15", "-maxMatch", "25", "-m", "5",
                           "-nproc", str(num_proc), "-out", out_file])


def main():
    if len(sys.argv) != 3:
        print("Usage: abruijn.py pre_assembly reads_file")
        return 1

    RAW_GENOME = "raw_genome.fasta"
    BLASR_AILNMENT = "blasr.m5"
    compose_raw_genome(sys.argv[1], RAW_GENOME)
    run_blasr(RAW_GENOME, sys.argv[2], NUM_PROC, BLASR_AILNMENT)

    return 0


if __name__ == "__main__":
    sys.exit(main())
