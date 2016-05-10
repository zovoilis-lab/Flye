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
import argparse

import abruijn.alignment as aln
import abruijn.bubbles as bbl
import abruijn.polish as pol
import abruijn.fasta_parser as fp
import abruijn.assemble as asm
from abruijn.__version__ import __version__


logger = logging.getLogger()


def run(args):
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    work_dir = os.path.abspath(args.out_dir)

    log_file = os.path.join(work_dir, "abruijn.log")
    enable_logging(log_file, args.debug)

    logger.info("Running ABruijn")
    aln.check_binaries()
    pol.check_binaries()
    asm.check_binaries()

    preassembly = os.path.join(work_dir, "read_edges.fasta")
    asm.assemble(args.reads, preassembly, args.kmer_size, args.min_cov,
                 args.max_cov, args.coverage)
    alignment, contigs_info, profile = \
            aln.get_alignment(preassembly, args.reads,
                              args.threads, work_dir)
    bubbles = bbl.get_bubbles(alignment, contigs_info)
    polished_seqs = pol.polish(bubbles, args.threads, work_dir, profile)
    out_genome = os.path.join(work_dir, "contigs.fasta")
    fp.write_fasta_dict(polished_seqs, out_genome)
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
    parser = argparse.ArgumentParser(description="ABruijn: assembly of long and"
                                     "error-prone reads")

    parser.add_argument("reads", metavar="reads",
                        help="path to a file with reads in FASTA format")
    parser.add_argument("out_dir", metavar="out_dir",
                        help="output directory")
    parser.add_argument("coverage", metavar="coverage", type=int,
                        help="estimated assembly coverage")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=1, help="number of parallel threads "
                        "(default: 1)")
    parser.add_argument("-k", "--kmer-size", dest="kmer_size", type=int,
                        default=15, help="kmer size (default: 15)")
    parser.add_argument("-m", "--min-cov", dest="min_cov", type=int,
                        default=None, help="minimum kmer coverage "
                        "(default: auto)")
    parser.add_argument("-x", "--max-cov", dest="max_cov", type=int,
                        default=None, help="maximum kmer coverage "
                        "(default: auto)")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    try:
        run(args)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException) as e:
        logger.error("Error: {0}".format(e))
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
