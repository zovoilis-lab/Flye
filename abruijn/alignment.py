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
import multiprocessing
import ctypes

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


class SynchronizedReader(object):
    def __init__(self, blasr_alignment):
        #will not be changed during exceution
        self.aln_path = blasr_alignment
        self.change_strand = True

        #will be shared between processes
        self.lock = multiprocessing.Lock()
        self.eof = multiprocessing.Value(ctypes.c_bool, False)
        self.position = multiprocessing.Value(ctypes.c_longlong, 0)

    def init_reading(self):
        """
        Call from the reading process, initializing local variables
        """
        if not os.path.exists(self.aln_path):
            raise AlignmentException("Can't open {0}".format(self.aln_path))
        self.aln_file = open(self.aln_path, "r")
        self.processed_contigs = set()

    def is_eof(self):
        return self.eof.value

    def get_chunk(self):
        """
        Alignment file is expected to be sorted!
        """
        #print os.getpid(), "Waiting for lock"
        with self.lock:
            #print os.getpid(), "Reading from ", self.position.value
            self.aln_file.seek(self.position.value)
            if self.eof.value:
                return None, []

            current_contig = None
            buffer = []
            while True:
                self.position.value = self.aln_file.tell()
                line = self.aln_file.readline()
                if not line:
                    break

                tokens = line.strip().split()
                if len(tokens) < 18:
                    raise AlignmentException("Error reading BLASR file")

                read_contig = tokens[5]
                if read_contig in self.processed_contigs:
                    raise AlignmentException("Alignment file is not sorted")

                err_rate = 1 - float(tokens[17].count("|")) / len(tokens[17])
                #self.errors.append(err_rate)

                if tokens[9] == "+" and self.change_strand:
                    trg_seq, qry_seq = tokens[16], tokens[18]
                else:
                    trg_seq = fp.reverse_complement(tokens[16])
                    qry_seq = fp.reverse_complement(tokens[18])
                aln = Alignment(tokens[0], tokens[5], int(tokens[2]),
                                int(tokens[3]), tokens[4],
                                int(tokens[1]), int(tokens[7]),
                                int(tokens[8]), tokens[9],
                                int(tokens[6]), trg_seq,
                                qry_seq, err_rate)

                if read_contig != current_contig:
                    prev_contig = current_contig
                    current_contig = read_contig

                    if prev_contig is not None:
                        self.processed_contigs.add(prev_contig)
                        #print os.getpid(), "Read", prev_contig, len(buffer)
                        return prev_contig, buffer
                    else:
                        buffer = [aln]
                else:
                    buffer.append(aln)

            #mean_err = float(sum(self.errors)) / len(self.errors)
            #logger.debug("Alignment error rate: {0}".format(mean_err))
            self.eof.value = True
            return current_contig, buffer


def check_binaries():
    if not which(BLASR_BIN):
        raise AlignmentException("BLASR is not installed")
    if not which("sort"):
        raise AlignmentException("UNIX sort utility is not available")


"""
def concatenate_contigs(contigs_file):
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
"""


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
    Runs BLASR and sort its output
    """
    _run_blasr(reference_file, reads_file, num_proc, out_alignment)
    logger.debug("Sorting alignment file")
    temp_file = out_alignment + "_sorted"
    subprocess.check_call(["sort", "-k", "6", out_alignment],
                          stdout=open(temp_file, "w"))
    os.remove(out_alignment)
    os.rename(temp_file, out_alignment)


def get_contigs_info(contigs_file):
    contigs_info = {}
    contigs_fasta = fp.read_fasta_dict(contigs_file)
    for ctg_id, ctg_seq in contigs_fasta.iteritems():
        contig_type = ctg_id.split("_")[0]
        contigs_info[ctg_id] = ContigInfo(ctg_id, len(ctg_seq),
                                          contig_type)

    return contigs_info


"""
def parse_alignment(alignment_file, ctg_id=None):
    circular_window = config.vals["circular_window"]
    alignment, _mean_error = _parse_blasr(alignment_file, change_strand=True,
                                          ctg_id=ctg_id)

    return alignment
"""


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


"""
def _parse_blasr(filename, change_strand, ctg_id):
    alignments = []
    errors = []
    with open(filename, "r") as f:
        for line in f:
            tokens = line.strip().split()

            if ctg_id is not None and tokens[5] != ctg_id:
                continue

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
"""


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
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AlignmentException(str(e))
