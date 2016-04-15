from __future__ import print_function
import os
import sys
from collections import namedtuple
import subprocess

import abruijn.fasta_parser as fp
from abruijn.utils import which


BLASR_BIN = "blasr"
CIRCULAR_WINDOW = 50000

Alignment = namedtuple("Alignment", ["qry_start", "qry_end", "qry_sign", "qry_len",
                                     "trg_start", "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate"])

class AlignmentException(Exception):
    pass


def check_binaries():
    if not which(BLASR_BIN):
        raise AlignmentException("Blasr is not installed")


def get_alignment(draft_file, reads_file, num_proc, work_dir):
    print("Aligning reads and formin mini-alignments", file=sys.stderr)
    raw_genome = os.path.join(work_dir, "draft_assembly.fasta")
    blasr_aln = os.path.join(work_dir, "alignment.m5")
    genome_len = _compose_raw_genome(draft_file, raw_genome)
    if not _run_blasr(raw_genome, reads_file, num_proc, blasr_aln):
        raise AlignmentException("Error in alignment module, exiting")
    return _parse_blasr(blasr_aln), genome_len


def _parse_blasr(filename):
    alignments = []
    errors = []
    with open(filename, "r") as f:
        for line in f:
            tokens = line.strip().split()
            err_rate = 1 - float(tokens[17].count("|")) / len(tokens[17])
            alignments.append(Alignment(int(tokens[2]), int(tokens[3]),
                                        tokens[4], int(tokens[1]),
                                        int(tokens[7]), int(tokens[8]),
                                        tokens[9], int(tokens[6]),
                                        tokens[16], tokens[18], err_rate))
            errors.append(err_rate)

    mean_err = float(sum(errors)) / len(errors)
    #print("Read error rate: {0:5.2f}".format(mean_err))
    return alignments


def _compose_raw_genome(contig_parts, out_file):
    """
    Concatenates contig parts and appends suffix from the beginning
    """
    fragment_index = {}
    genome_framents = fp.read_fasta_dict(contig_parts)
    for h, s in genome_framents.iteritems():
        idx = int(h.split("_")[1])
        fragment_index[idx] = s
    genome_seq = "".join(map(fragment_index.get, sorted(fragment_index)))
    genome_len = len(genome_seq)
    assert len(genome_seq) > CIRCULAR_WINDOW
    genome_seq += genome_seq[:CIRCULAR_WINDOW]
    fp.write_fasta_dict({"contig_1" : genome_seq}, out_file)
    return genome_len


def _run_blasr(reference_file, reads_file, num_proc, out_file):
    try:
        subprocess.check_call([BLASR_BIN, reads_file, reference_file,
                               "-bestn", "1", "-minMatch", "15",
                               "-maxMatch", "25", "-m", "5",
                               "-nproc", str(num_proc), "-out", out_file])
    except (subprocess.CalledProcessError, OSError) as e:
        print("Error while running blasr: " + str(e))
        return False
    return True
