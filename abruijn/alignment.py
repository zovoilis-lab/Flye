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
        raise AlignmentException("BLASR is not installed")


def concatenate_contigs(contigs_file):
    """
    Concatenates contig parts output by assembly module
    """
    genome_framents = fp.read_fasta_dict(contigs_file)
    contig_types = {}
    by_contig = defaultdict(list)
    for h, s in genome_framents.iteritems():
        tokens = h.split("_")
        cont_type = tokens[0]
        contig_id = tokens[0] + "_" + tokens[1]
        part_id = int(tokens[3])
        contig_types[contig_id] = cont_type
        by_contig[contig_id].append((part_id, s))

    contigs_fasta = {}
    for contig_id, contig_seqs in by_contig.iteritems():
        seqs_sorted = sorted(contig_seqs, key=lambda p: p[0])
        contig_concat = "".join(map(lambda p: p[1], seqs_sorted))
        contig_len = len(contig_concat)
        contigs_fasta[contig_id] = contig_concat

    return contigs_fasta


def make_blasr_reference(contigs_fasta, out_file):
    """
    Outputs 'reference' for BLASR run, appends a suffix to circular contigs
    """
    circular_window = config.vals["circular_window"]
    for contig_id in contigs_fasta:
        contig_type = contig_id.split("_")[0]
        if (contig_type == "circular" and
            len(contigs_fasta[contig_id]) > circular_window):
            contigs_fasta[contig_id] += \
                    contigs_fasta[contig_id][:circular_window]

    fp.write_fasta_dict(contigs_fasta, out_file)


def make_alignment(reference_file, reads_file, num_proc,
                   out_alignment):
    """
    Runs BLASR
    """
    logger.info("Running BLASR")
    _run_blasr(reference_file, reads_file, num_proc, out_alignment)


def parse_alignment(alignment_file):
    """
    Parses BLASR alignment and choses error profile base on error rate
    """
    circular_window = config.vals["circular_window"]
    alignment, mean_error = _parse_blasr(alignment_file, change_strand=True)
    contigs_info = {}

    for aln in alignment:
        if aln.trg_id not in contigs_info:
            contig_type = aln.trg_id.split("_")[0]
            contig_len = aln.trg_len
            if contig_type == "circular" and contig_len > circular_window:
                contig_len -= circular_window
            contigs_info[aln.trg_id] = ContigInfo(aln.trg_id, contig_len,
                                                  contig_type)

    return alignment, contigs_info, mean_error


def choose_error_mode(err_rate):
    """
    Choses error mode base on error rate
    """
    if err_rate < config.vals["err_rate_threshold"]:
        profile = "pacbio"
    else:
        profile = "nano"

    logger.info("Chosen '{0}' error mode".format(profile))
    return profile


def shift_gaps(seq_trg, seq_qry):
    """
    Shifts all ambigious query gaps to the right
    """
    lst_trg, lst_qry = list("$" + seq_trg + "$"), list("$" + seq_qry + "$")
    is_gap = False
    gap_start = 0
    for i in xrange(len(lst_trg)):
        if is_gap and lst_qry[i] != "-":
            is_gap = False
            swap_left = gap_start - 1
            swap_right = i - 1

            while (swap_left > 0 and swap_right >= gap_start and
                   lst_qry[swap_left] == lst_trg[swap_right]):
                lst_qry[swap_left], lst_qry[swap_right] = \
                            lst_qry[swap_right], lst_qry[swap_left]
                swap_left -= 1
                swap_right -= 1

        if not is_gap and lst_qry[i] == "-":
            is_gap = True
            gap_start = i

    return "".join(lst_qry[1 : -1])


def _parse_blasr(filename, change_strand):
    """
    Parse Blasr output
    """
    alignments = []
    errors = []
    with open(filename, "r") as f:
        for line in f:
            tokens = line.strip().split()
            err_rate = 1 - float(tokens[17].count("|")) / len(tokens[17])
            if tokens[9] == "+" and change_strand:
                trg_seq, qry_seq = tokens[16], tokens[18]
            else:
                trg_seq = fp.reverse_complement(tokens[16])
                qry_seq = fp.reverse_complement(tokens[18])

            alignments.append(Alignment(tokens[0], tokens[5], int(tokens[2]),
                                        int(tokens[3]), tokens[4],
                                        int(tokens[1]), int(tokens[7]),
                                        int(tokens[8]), tokens[9],
                                        int(tokens[6]), trg_seq,
                                        qry_seq, err_rate))
            errors.append(err_rate)

    mean_err = float(sum(errors)) / len(errors)
    logger.debug("Alignment error rate: {0}".format(mean_err))
    return alignments, mean_err


def _guess_blasr_version():
    """
    Tries to guess whether we need one or two dashed in command line
    """
    try:
        devnull = open(os.devnull, "w")
        stdout = subprocess.check_output([BLASR_BIN, "-version"])
    except subprocess.CalledProcessError as e:
        return "pacbio_new"

    version = stdout.splitlines()[0].split()[-1]
    tokens = version.split(".", 1)
    if len(tokens) == 1:    #probably original one
        return "original"

    major, minor = tokens[0], tokens[1].split(".", 1)[0]
    if int(major) < 2:
        return "original"
    if int(major) < 5 or int(minor) < 1:
        return "pacbio_old"

    return "pacbio_new"


def _run_blasr(reference_file, reads_file, num_proc, out_file):
    cmdline = [BLASR_BIN, reads_file, reference_file,
               "-bestn", "1", "-minMatch", "15",
               "-maxMatch", "20", "-m", "5",
               "-nproc", str(num_proc), "-out", out_file,
               "-advanceHalf", "-advanceExactMatches", "10"]

    blasr_version = _guess_blasr_version()
    if blasr_version in ["pacbio_old", "pacbio_new"]:
        cmdline.extend(["-fastSDP", "-aggressiveIntervalCut"])
    if blasr_version == "pacbio_new":
        cmdline = map(lambda cmd: cmd.replace("-", "--")
                      if cmd.startswith("-") and len(cmd) > 2 else cmd,
                      cmdline)

    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call(cmdline, stderr=devnull)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error("While running blasr: " + str(e))
        raise AlignmentException("Error in alignment module, exiting")
