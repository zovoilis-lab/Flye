#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs Minimap2 and parses its output
"""

import os
import re
import sys
from collections import namedtuple, defaultdict
import subprocess
import logging
import multiprocessing
import ctypes

import flye.utils.fasta_parser as fp
from flye.utils.utils import which
import flye.config.py_cfg as cfg


logger = logging.getLogger()
MINIMAP_BIN = "flye-minimap2"

Alignment = namedtuple("Alignment", ["qry_id", "trg_id", "qry_start", "qry_end",
                                     "qry_sign", "qry_len", "trg_start",
                                     "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate",
                                     "is_secondary"])

ContigInfo = namedtuple("ContigInfo", ["id", "length", "type"])


class AlignmentException(Exception):
    pass


class PafHit:
    """
    Stores paf alignment
    """
    __slots__ = ("query", "query_length", "query_start", "query_end",
                 "target", "target_length", "target_start", "target_end")
    def __init__(self, raw_hit):
        hit = raw_hit.split()

        self.query = hit[0]
        self.query_length = int(hit[1])
        self.query_start = int(hit[2])
        self.query_end = int(hit[3])

        self.target = hit[5]
        self.target_length = int(hit[6])
        self.target_start = int(hit[7])
        self.target_end = int(hit[8])

    def query_mapping_length(self):
        return self.query_end - self.query_start + 1

    def target_mapping_length(self):
        return self.target_end - self.target_start + 1

    def query_left_overhang(self):
        return self.query_start

    def query_right_overhang(self):
        return self.query_length - self.query_end + 1

    def target_left_overhang(self):
        return self.target_start

    def target_right_overhang(self):
        return self.target_length - self.target_end + 1


def read_paf(filename):
    """
    Streams out paf alignments
    """
    hits = []
    with open(filename) as f:
        for raw_hit in f:
            yield PafHit(raw_hit)


def read_paf_grouped(filename):
    """
    Outputs chunks of alignments for each (query, target)pair.
    Assumes that PAF alignment is already sorted by query.
    """
    prev_hit = None
    target_hits = defaultdict(list)
    for hit in read_paf(filename):
        if prev_hit is not None and hit.query != prev_hit.query:
            for trg in sorted(target_hits):
                yield target_hits[trg]
            target_hits = defaultdict(list)

        target_hits[hit.target].append(hit)
        prev_hit = hit

    if len(target_hits):
        for trg in sorted(target_hits):
            yield target_hits[trg]


class SynchronizedSamReader(object):
    """
    Parses SAM file in multiple threads.
    """
    def __init__(self, sam_alignment, reference_fasta,
                 max_coverage=None, use_secondary=False):
        #will not be changed during exceution, each process has its own copy
        self.aln_path = sam_alignment
        self.aln_file = None
        self.ref_fasta = reference_fasta
        self.change_strand = True
        self.max_coverage = max_coverage
        self.seq_lengths = {}
        self.use_secondary = use_secondary

        #reading SAM header
        if not os.path.exists(self.aln_path):
            raise AlignmentException("Can't open {0}".format(self.aln_path))

        with open(self.aln_path, "r") as f:
            for line in f:
                if not line or not _is_sam_header(line):
                    break
                if line.startswith("@SQ"):
                    seq_name = None
                    seq_len = None
                    for tag in line.split():
                        if tag.startswith("SN"):
                            seq_name = tag[3:]
                        if tag.startswith("LN"):
                            seq_len = int(tag[3:])
                    if seq_name and seq_len:
                        self.seq_lengths[seq_name] = seq_len

        #will be shared between processes
        self.lock = multiprocessing.Lock()
        self.eof = multiprocessing.Value(ctypes.c_bool, False)
        self.position = multiprocessing.Value(ctypes.c_longlong, 0)

    def init_reading(self):
        """
        Call from the reading process, initializing local variables
        """
        self.aln_file = open(self.aln_path, "r")
        self.processed_contigs = set()
        self.cigar_parser = re.compile("[0-9]+[MIDNSHP=X]")

    def stop_reading(self):
        """
        Call when the reading is done
        """
        self.aln_file.close()

    def is_eof(self):
        return self.eof.value

    def parse_cigar(self, cigar_str, read_str, ctg_name, ctg_pos):
        ctg_str = self.ref_fasta[ctg_name]
        trg_seq = []
        qry_seq = []
        trg_start = ctg_pos - 1
        trg_pos = ctg_pos - 1
        qry_start = 0
        qry_pos = 0

        left_hard = True
        left_soft = True
        hard_clipped_left = 0
        hard_clipped_right = 0
        soft_clipped_left = 0
        soft_clipped_right = 0
        for token in self.cigar_parser.findall(cigar_str):
            size, op = int(token[:-1]), token[-1]
            if op == "H":
                if left_hard:
                    qry_start += size
                    hard_clipped_left += size
                else:
                    hard_clipped_right += size
            elif op == "S":
                qry_pos += size
                if left_soft:
                    soft_clipped_left += size
                else:
                    soft_clipped_right += size
            elif op == "M":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                qry_pos += size
                trg_pos += size
            elif op == "I":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append("-" * size)
                qry_pos += size
            elif op == "D":
                qry_seq.append("-" * size)
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                trg_pos += size
            else:
                raise AlignmentException("Unsupported CIGAR operation: " + op)
            left_hard = False
            if op != "H":
                left_soft = False

        trg_seq = "".join(trg_seq)
        qry_seq = "".join(qry_seq)
        matches = 0
        for i in xrange(len(trg_seq)):
            if trg_seq[i] == qry_seq[i]:
                matches += 1
        err_rate = 1 - float(matches) / len(trg_seq)

        trg_end = trg_pos
        qry_end = qry_pos + hard_clipped_left
        qry_len = qry_end + hard_clipped_right
        qry_start += soft_clipped_left
        qry_end -= soft_clipped_right

        return (trg_start, trg_end, len(ctg_str), trg_seq,
                qry_start, qry_end, qry_len, qry_seq, err_rate)

    def get_chunk(self):
        """
        Alignment file is expected to be sorted!
        """

        buffer = []
        parsed_contig = None

        with self.lock:
            self.aln_file.seek(self.position.value)
            if self.eof.value:
                return None, []

            current_contig = None
            while True:
                self.position.value = self.aln_file.tell()
                line = self.aln_file.readline()
                if not line: break
                if _is_sam_header(line): continue

                tokens = line.strip().split()
                if len(tokens) < 11:
                    continue
                    #raise AlignmentException("Error reading SAM file")

                read_contig = tokens[2]
                flags = int(tokens[1])
                is_unmapped = flags & 0x4
                is_secondary = flags & 0x100
                is_supplementary = flags & 0x800    #allow supplementary

                #if is_unmapped or is_secondary: continue
                if is_unmapped: continue
                if is_secondary and not self.use_secondary: continue
                if read_contig in self.processed_contigs:
                    raise AlignmentException("Alignment file is not sorted")

                if read_contig != current_contig:
                    prev_contig = current_contig
                    current_contig = read_contig

                    if prev_contig is not None:
                        self.processed_contigs.add(prev_contig)
                        parsed_contig = prev_contig
                        break
                    else:
                        buffer = [tokens]
                else:
                    buffer.append(tokens)

            if not parsed_contig:
                self.eof.value = True
                parsed_contig = current_contig
        #end with

        sequence_length = 0
        alignments = []
        for tokens in buffer:
            read_id = tokens[0]
            read_contig = tokens[2]
            cigar_str = tokens[5]
            read_str = tokens[9]
            ctg_pos = int(tokens[3])
            flags = int(tokens[1])
            is_reversed = flags & 0x16
            is_secondary = flags & 0x100

            if read_str == "*":
                raise Exception("Error parsing SAM: record without read sequence")

            (trg_start, trg_end, trg_len, trg_seq,
            qry_start, qry_end, qry_len, qry_seq, err_rate) = \
                    self.parse_cigar(cigar_str, read_str, read_contig, ctg_pos)

            #OVERHANG = cfg.vals["read_aln_overhang"]
            #if (float(qry_end - qry_start) / qry_len > self.min_aln_rate or
            #        trg_start < OVERHANG or trg_len - trg_end < OVERHANG):
            aln = Alignment(read_id, read_contig, qry_start,
                            qry_end, "-" if is_reversed else "+",
                            qry_len, trg_start, trg_end, "+", trg_len,
                            qry_seq, trg_seq, err_rate, is_secondary)
            alignments.append(aln)

            sequence_length += qry_end - qry_start
            #In rare cases minimap2 does not output SQ tag, so need to check
            if parsed_contig in self.seq_lengths:
                contig_length = self.seq_lengths[parsed_contig]
                if sequence_length / contig_length > self.max_coverage:
                    break

        return parsed_contig, alignments


def check_binaries():
    if not which(MINIMAP_BIN):
        raise AlignmentException("Minimap2 is not installed")
    if not which("sort"):
        raise AlignmentException("UNIX sort utility is not available")


def make_alignment(reference_file, reads_file, num_proc,
                   work_dir, platform, out_alignment, reference_mode,
                   sam_output):
    """
    Runs minimap2 and sorts its output
    """
    minimap_ref_mode = {False: "ava", True: "map"}
    minimap_reads_mode = {"nano": "ont", "pacbio": "pb"}
    mode = minimap_ref_mode[reference_mode] + "-" + minimap_reads_mode[platform]

    _run_minimap(reference_file, reads_file, num_proc, mode,
                 out_alignment, sam_output)

    if sam_output:
        _preprocess_sam(out_alignment, work_dir)


def _preprocess_sam(sam_file, work_dir):
    """
    Proprocesses minimap2 output by adding SEQ
    to secondary alignments, removing
    unaligned reads and then sorting
    file by reference sequence id
    """
    expanded_sam = sam_file + "_expanded"
    merged_file = sam_file + "_merged"
    sorted_file = sam_file + "_sorted"

    #puting SAM headers to the final postprocessed file first
    with open(sam_file, "r") as hdr_in, open(merged_file, "w") as fout:
        for line in hdr_in:
            if not _is_sam_header(line):
                break
            fout.write(line)

    #adding SEQ fields to secondary alignments
    with open(sam_file, "r") as fin, open(expanded_sam, "w") as fout:
        prev_id = None
        prev_seq = None
        primary_reversed = None
        for line in fin:
            if _is_sam_header(line):
                continue

            tokens = line.strip().split()
            flags = int(tokens[1])
            is_unmapped = flags & 0x4
            is_secondary = flags & 0x100
            is_supplementary = flags & 0x800
            is_reversed = flags & 0x16

            if is_unmapped:
                continue

            read_id, cigar_str, read_seq = tokens[0], tokens[5], tokens[9]
            has_hard_clipped = "H" in cigar_str

            #Checking format assumptions
            if has_hard_clipped:
                if is_secondary:
                    raise Exception("Secondary alignment with hard-clipped bases")
                if not is_supplementary:
                    raise Exception("Primary alignment with hard-clipped bases")
            if not is_secondary and read_seq == "*":
                raise Exception("Missing SEQ for non-secondary alignment")

            if read_seq == "*":
                if read_id != prev_id:
                    raise Exception("SAM file is not sorted by read names")
                if is_reversed == primary_reversed:
                    tokens[9] = prev_seq
                else:
                    tokens[9] = fp.reverse_complement(prev_seq)

            #Assuming that the first read alignmnent in SAM is primary
            elif prev_id != read_id:
                if has_hard_clipped:
                    raise Exception("Hard clipped bases in the primamry read")
                prev_id = read_id
                prev_seq = read_seq
                primary_reversed = is_reversed

            fout.write("\t".join(tokens) + "\n")

    #don't need the original SAM anymore, cleaning up space
    os.remove(sam_file)

    #logger.debug("Sorting alignment file")
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    subprocess.check_call(["sort", "-k", "3,3", "-T", work_dir, expanded_sam],
                          stdout=open(sorted_file, "w"), env=env)

    #don't need the expanded file anymore
    os.remove(expanded_sam)

    #appending to the final file, that already contains headers
    with open(sorted_file, "r") as sort_in, open(merged_file, "a") as fout:
        for line in sort_in:
            if not _is_sam_header(line):
                fout.write(line)

    os.remove(sorted_file)
    os.rename(merged_file, sam_file)


def get_contigs_info(contigs_file):
    contigs_info = {}
    contigs_fasta = fp.read_sequence_dict(contigs_file)
    for ctg_id, ctg_seq in contigs_fasta.iteritems():
        contig_type = ctg_id.split("_")[0]
        contigs_info[ctg_id] = ContigInfo(ctg_id, len(ctg_seq),
                                          contig_type)

    return contigs_info


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


def get_uniform_alignments(alignments, seq_len):
    """
    Leaves top alignments for each position within contig
    assuming uniform coverage distribution
    """
    def _get_median(lst):
        if not lst:
            raise ValueError("_get_median() arg is an empty sequence")
        sorted_list = sorted(lst)
        if len(lst) % 2 == 1:
            return sorted_list[len(lst)/2]
        else:
            mid1 = sorted_list[(len(lst)/2) - 1]
            mid2 = sorted_list[(len(lst)/2)]
            return float(mid1 + mid2) / 2

    WINDOW = 100
    MIN_COV = 10
    COV_RATE = 1.25

    #split contig into windows, get median read coverage over all windows and
    #determine the quality threshold cutoffs for each window
    wnd_primary_cov = [0 for _ in xrange(seq_len / WINDOW + 1)]
    wnd_aln_quality = [[] for _ in xrange(seq_len / WINDOW + 1)]
    wnd_qual_thresholds = [1.0 for _ in xrange(seq_len / WINDOW + 1)]
    for aln in alignments:
        for i in xrange(aln.trg_start / WINDOW, aln.trg_end / WINDOW):
            if not aln.is_secondary:
                wnd_primary_cov[i] += 1
            wnd_aln_quality[i].append(aln.err_rate)

    #for each window, select top X alignmetns, where X is the median read coverage
    cov_threshold = max(int(COV_RATE * _get_median(wnd_primary_cov)), MIN_COV)
    for i in xrange(len(wnd_aln_quality)):
        if len(wnd_aln_quality[i]) > cov_threshold:
            wnd_qual_thresholds[i] = sorted(wnd_aln_quality[i])[cov_threshold]

    #for each alignment, count in how many windows it passes the threshold
    filtered_alignments = []
    total_sequence = 0
    filtered_sequence = 0
    for aln in alignments:
        good_windows = 0
        total_windows = aln.trg_end / WINDOW - aln.trg_start / WINDOW
        total_sequence += aln.trg_end - aln.trg_start
        for i in xrange(aln.trg_start / WINDOW, aln.trg_end / WINDOW):
            if aln.err_rate <= wnd_qual_thresholds[i]:
                good_windows += 1

        if good_windows > total_windows / 2:
            filtered_alignments.append(aln)
            filtered_sequence += aln.trg_end - aln.trg_start

    filtered_reads_rate = 1 - float(len(filtered_alignments)) / len(alignments)
    filtered_seq_rate = 1 - float(filtered_sequence) / total_sequence
    #logger.debug("Filtered {0:7.2f}% reads, {1:7.2f}% sequence"
    #                .format(filtered_reads_rate * 100, filtered_seq_rate * 100))

    return filtered_alignments


def split_into_chunks(fasta_in, chunk_size):
    out_dict = {}
    for header, seq in fasta_in.iteritems():
        #print len(seq)
        for i in xrange(0, max(len(seq) / chunk_size, 1)):
            chunk_hdr = "{0}$chunk_{1}".format(header, i)
            start = i * chunk_size
            end = (i + 1) * chunk_size
            if len(seq) - end < chunk_size:
                end = len(seq)

            #print(start, end)
            out_dict[chunk_hdr] = seq[start : end]

    return out_dict


def merge_chunks(fasta_in, fold_function=lambda l: "".join(l)):
    """
    Merges sequence chunks. Chunk names are in format `orig_name$chunk_id`.
    Each chunk is as dictionary entry. Value type is arbitrary and
    one can supply a custom fold function
    """
    def name_split(h):
        orig_hdr, chunk_id = h.rsplit("$", 1)
        return orig_hdr, int(chunk_id.rsplit("_", 1)[1])

    out_dict = {}
    cur_seq = []
    cur_contig = None
    for hdr in sorted(fasta_in, key=name_split):
        orig_name, chunk_id = name_split(hdr)
        if orig_name != cur_contig:
            if cur_contig != None:
                out_dict[cur_contig] = fold_function(cur_seq)
            cur_seq = []
            cur_contig = orig_name
        cur_seq.append(fasta_in[hdr])

    if cur_seq:
        out_dict[cur_contig] = fold_function(cur_seq)

    return out_dict


def _run_minimap(reference_file, reads_files, num_proc, mode, out_file,
                 sam_output):
    cmdline = [MINIMAP_BIN, reference_file]
    cmdline.extend(reads_files)
    cmdline.extend(["-x", mode, "-t", str(num_proc)])
    if sam_output:
        #a = SAM output, p = min primary-to-seconday score
        #N = max secondary alignments
        cmdline.extend(["-a", "-p", "0.7", "-N", "10"])

    try:
        devnull = open(os.devnull, "w")
        #logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline, stderr=devnull,
                              stdout=open(out_file, "w"))
    except (subprocess.CalledProcessError, OSError) as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AlignmentException(str(e))


def _is_sam_header(line):
    return line[:3] in ["@PG", "@HD", "@SQ", "@RG", "@CO"]
