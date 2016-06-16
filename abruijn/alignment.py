#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs Blasr aligner and parses its output
"""

import os
import sys
from collections import namedtuple, defaultdict
import subprocess
import logging

import abruijn.fasta_parser as fp
from abruijn.utils import which
import abruijn.config as config


logger = logging.getLogger()
BLASR_BIN = "blasr"

Alignment = namedtuple("Alignment", ["qry_id", "trg_id", "qry_start", "qry_end",
                                     "qry_sign", "qry_len", "trg_start",
                                     "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate"])
ContigInfo = namedtuple("ContigInfo", ["id", "length", "type"])

class AlignmentException(Exception):
    pass


def check_binaries():
    if not which(BLASR_BIN):
        raise AlignmentException("Blasr is not installed")


def compose_raw_genome(contig_parts, out_file):
    """
    Concatenates contig parts and appends suffix from the beginning
    """
    CIRCULAR_WINDOW = 50000

    genome_framents = fp.read_fasta_dict(contig_parts)
    contig_types = {}
    by_contig = defaultdict(list)
    for h, s in genome_framents.iteritems():
        tokens = h.split("_")
        cont_type = tokens[0]
        contig_id = tokens[0] + "_" + tokens[1]
        part_id = int(tokens[3])
        contig_types[contig_id] = cont_type
        by_contig[contig_id].append((part_id, s))

    contigs_info = {}
    contigs_fasta = {}
    for contig_id, contig_seqs in by_contig.iteritems():
        seqs_sorted = sorted(contig_seqs, key=lambda p: p[0])
        contig_concat = "".join(map(lambda p: p[1], seqs_sorted))
        contig_len = len(contig_concat)
        if contig_types[contig_id] == "circular":
            if len(contig_concat) > CIRCULAR_WINDOW:
                contig_concat += contig_concat[:CIRCULAR_WINDOW]
        contigs_fasta[contig_id] = contig_concat
        contigs_info[contig_id] = ContigInfo(contigs_info, contig_len,
                                             contig_types[contig_id])

    fp.write_fasta_dict(contigs_fasta, out_file)
    return contigs_info


def get_alignment(reference_file, reads_file, num_proc,
                  work_dir, iter_id):
    """
    Runs blasr and return parsed alignment
    """
    logger.info("Running BLASR")
    blasr_file = os.path.join(work_dir, "alignment_" + iter_id + ".m5")
    _run_blasr(reference_file, reads_file, num_proc, blasr_file)
    alignment, mean_err = _parse_blasr(blasr_file)
    profile = _choose_profile(mean_err)
    return alignment, profile


def _choose_profile(err_rate):
    if err_rate < config.vals["err_rate_threshold"]:
        profile = "pacbio"
    else:
        profile = "nano"

    logger.debug("Alignment error rate: {0}".format(err_rate))
    logger.info("Chosen '{0}' polishing profile".format(profile))
    return profile


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


def _run_blasr(reference_file, reads_file, num_proc, out_file):
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([BLASR_BIN, reads_file, reference_file,
                               "-bestn", "1", "-minMatch", "15",
                               "-maxMatch", "25", "-m", "5",
                               "-nproc", str(num_proc), "-out", out_file],
                               stderr=devnull)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error("While running blasr: " + str(e))
        raise AlignmentException("Error in alignment module, exiting")
