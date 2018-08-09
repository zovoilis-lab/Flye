#(c) 2013-2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some basic FASTA I/O
"""

import logging

from string import maketrans

logger = logging.getLogger()

class FastaError(Exception):
    pass

def read_fasta_dict(filename):
    """
    Reads fasta file into dictionary. Also performs some validation
    """
    #logger.debug("Reading contigs file")

    header = None
    seq = []
    fasta_dict = {}

    try:
        with open(filename, "r") as f:
            for lineno, line in enumerate(f):
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        fasta_dict[header] = "".join(seq)
                        seq = []
                    header = line[1:].split(" ")[0]
                else:
                    if not _validate_seq(line):
                        raise FastaError("Invalid char in \"{0}\" at line {1}"
                                         .format(filename, lineno))
                    seq.append(line)

            if header and len(seq):
                fasta_dict[header] = "".join(seq)

    except IOError as e:
        raise FastaError(e)

    return fasta_dict


def write_fasta_dict(fasta_dict, filename):
    """
    Writes dictionary with fasta to file
    """
    with open(filename, "w") as f:
        for header in sorted(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(fasta_dict[header]), 60):
                f.write(fasta_dict[header][i:i + 60] + "\n")

def read_fast_file(filename):
    fasta_bool = False
    fastq_bool = False
    try:
        with open(filename, "r") as f:
            line = f.readline()
            if line.startswith(">"):
                fasta_bool = True
            elif line.startswith("@"):
                fastq_bool = True
            else:
                raise FastaError("{0} not FASTA or FASTQ file".format(filename))
    except IOError as e:
        raise FastaError(e)
    
    if fasta_bool:
        return read_fasta_dict(filename)
    elif fastq_bool:
        fastq_dict = read_fastq_dict(filename)
        return {header:fastq_dict[header][0] for header in fastq_dict}

def read_fastq_dict(filename):
    """
    Reads fastq file into dictionary. Also performs some validation
    """
    #logger.debug("Reading contigs file")

    header = None
    seq = []
    qual = []
    fastq_dict = {}
    seq_bool = False
    qual_bool = False

    try:
        with open(filename, "r") as f:
            for lineno, line in enumerate(f):
                line = line.strip()
                if line.startswith("@"):
                    if header:
                        fastq_dict[header] = ("".join(seq), "".join(qual))
                        seq = []
                        qual = []
                    header = line[1:].split(" ")[0]
                    seq_bool = True
                    qual_bool = False
                elif seq_bool and _validate_seq(line):
                    seq.append(line)
                elif line.startswith("+"):
                    seq_bool = False
                    qual_bool = True
                elif qual_bool and _validate_qual(line):
                    qual.append(line)
                else:
                    raise FastaError("Invalid char in \"{0}\" at line {1}"
                                         .format(filename, lineno))

            if header and len(seq) and len(qual):
                fastq_dict[header] = ("".join(seq), "".join(qual))

    except IOError as e:
        raise FastaError(e)

    return fastq_dict

COMPL = maketrans("ATGCURYKMSWBVDHNXatgcurykmswbvdhnx",
                  "TACGAYRMKSWVBHDNXtacgayrmkswvbhdnx")
def reverse_complement(string):
    return string[::-1].translate(COMPL)


def _validate_seq(sequence):
    VALID_CHARS = "ACGTURYKMSWBDHVNXatgcurykmswbvdhnx"
    if len(sequence.translate(None, VALID_CHARS)):
        return False
    return True

def _validate_qual(qualities):
    VALID_CHARS = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    if len(qualities.translate(None, VALID_CHARS)):
        return False
    return True