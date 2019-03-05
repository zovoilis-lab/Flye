#(c) 2013-2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some basic FASTA I/O
"""

import logging
import gzip
import io
import os

from string import maketrans

logger = logging.getLogger()


class FastaError(Exception):
    pass


def read_sequence_dict(filename):
    """
    Reads Fasta/q file (could be gzip'ed) into a dictionary
    """
    try:
        seq_dict = {}
        gzipped, fastq = _is_fastq(filename)

        if not gzipped:
            handle = open(filename, "r")
        else:
            handle = os.popen("zcat {0}".format(filename))
            #gz = gzip.open(filename, "rb")
            #handle = io.BufferedReader(gz)

        if fastq:
            for hdr, seq, qual in _read_fastq(handle):
                if not _validate_seq(seq):
                    raise FastaError("Invalid char while reading {0}"
                                     .format(filename))
                seq_dict[hdr] = to_acgt(seq)
        else:
            for hdr, seq in _read_fasta(handle):
                if not _validate_seq(seq):
                    raise FastaError("Invalid char while reading {0}"
                                     .format(filename))
                seq_dict[hdr] = to_acgt(seq)

        return seq_dict

    except IOError as e:
        raise FastaError(e)


def read_sequence_lengths(filename):
    """
    Reads length of Fasta/q sequences (could be gzip'ed) into a dictionary
    """
    try:
        length_dict = {}
        gzipped, fastq = _is_fastq(filename)

        if not gzipped:
            handle = open(filename, "r")
        else:
            handle = os.popen("zcat {0}".format(filename))
            #gz = gzip.open(filename, "rb")
            #handle = io.BufferedReader(gz)

        if fastq:
            for hdr, seq, qual in _read_fastq(handle):
                if not _validate_seq(seq):
                    raise FastaError("Invalid char while reading {0}"
                                     .format(filename))
                length_dict[hdr] = len(seq)
        else:
            for hdr, seq in _read_fasta(handle):
                if not _validate_seq(seq):
                    raise FastaError("Invalid char while reading {0}"
                                     .format(filename))
                length_dict[hdr] = len(seq)

        return length_dict

    except IOError as e:
        raise FastaError(e)


def _is_fastq(filename):
    suffix = filename.rsplit(".")[-1]
    without_gz = filename
    gzipped = False

    if suffix == "gz":
        gzipped = True
        without_gz = filename.rstrip(".gz")

    suffix = without_gz.rsplit(".")[-1]
    if suffix in ["fasta", "fa"]:
        return gzipped, False

    if suffix in ["fastq", "fq"]:
        return gzipped, True

    raise FastaError("Unkown file extension: " + filename)


def _read_fasta(file_handle):
    header = None
    seq = []

    for line in file_handle:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if header:
                yield header, "".join(seq)
                seq = []
            header = line[1:].split(" ")[0]
        else:
            seq.append(line)

    if header and len(seq):
        yield header, "".join(seq)


def _read_fastq(file_handle):
    seq = None
    qual = None
    header = None
    state_counter = 0

    for no, line in enumerate(file_handle):
        line = line.strip()
        if not line:
            continue

        if state_counter == 0:
            if line[0] != "@":
                raise FastaError("Fastq format error: {0} at line {1}"
                                    .format(file_handle.name, no))
            header = line[1:].split(" ")[0]

        if state_counter == 1:
            seq = line

        if state_counter == 2:
            if line[0] != "+":
                raise FastaError("Fastq format error: {0} at line {1}"
                                    .format(file_handle.name, no))

        if state_counter == 3:
            qual = line
            yield header, seq, qual

        state_counter = (state_counter + 1) % 4


def write_fasta_dict(fasta_dict, filename):
    """
    Writes dictionary with fasta to file
    """
    with open(filename, "w") as f:
        for header in sorted(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(fasta_dict[header]), 60):
                f.write(fasta_dict[header][i:i + 60] + "\n")


COMPL = maketrans("ATGCURYKMSWBVDHNXatgcurykmswbvdhnx",
                  "TACGAYRMKSWVBHDNXtacgayrmkswvbhdnx")
def reverse_complement(string):
    return string[::-1].translate(COMPL)


VALID_CHARS = "ACGTURYKMSWBDHVNXatgcurykmswbvdhnx"
def _validate_seq(sequence):
    if len(sequence.translate(None, VALID_CHARS)):
        return False
    return True


ACGT_CHARS = "ACGTacgt"
TO_ACGT = maketrans("URYKMSWBVDHNXurykmswbvdhnx",
                    "ACGTACGTACGTAacgtacgtacgta")
def to_acgt(dna_str):
    """
    assumes tha all characters are valid
    """
    if len(dna_str.translate(None, ACGT_CHARS)) == 0:
        return dna_str
    else:
        if not to_acgt.ACGT_WARN:
            to_acgt.ACGT_WARN = True
            logger.warning("Input contain non-ACGT characters - "
                           "they will be converted to arbitrary ACGTs")
        return dna_str.translate(TO_ACGT)
to_acgt.ACGT_WARN = False
