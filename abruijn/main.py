#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Main logic of the package
"""

from __future__ import print_function
import sys
import os
import logging

import abruijn.alignment as aln
import abruijn.bubbles as bbl
import abruijn.polish as pol
import abruijn.fasta_parser as fp
import abruijn.assemble as asm


logger = logging.getLogger()


def run(reads_file, work_dir, num_proc, debug):
    log_file = os.path.join(work_dir, "abruijn.log")
    enable_logging(log_file, debug)

    logger.info("Running ABruijn")
    aln.check_binaries()
    pol.check_binaries()
    asm.check_binaries()

    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    preassembly = os.path.join(work_dir, "read_edges.fasta")
    logger.info("Assembling reads")
    asm.assemble(reads_file, preassembly)
    alignment, genome_len = aln.get_alignment(preassembly, reads_file,
                                              num_proc, work_dir)
    bubbles = bbl.get_bubbles(alignment, genome_len)
    logger.info("Polishing draft assembly")
    polished_seq = pol.polish(bubbles, num_proc, work_dir)
    out_genome = os.path.join(work_dir, "contigs.fasta")
    fp.write_fasta_dict({"contig_1": polished_seq}, out_genome)
    logger.info("Done! Your assembly is in file: " + out_genome)


def enable_logging(log_file, debug):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def main():
    DEBUG = True
    NUM_PROC = 8
    if len(sys.argv) != 3:
        print("Usage: abruijn.py reads_file out_dir")
        return 1

    try:
        run(sys.argv[1], sys.argv[2], NUM_PROC, DEBUG)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException) as e:
        print("Error: ", e, file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
