import os
import logging
import subprocess

import flye.fasta_parser as fp


logger = logging.getLogger()
MINIMAP_BIN = "flye-minimap2"


class ShortPlasmidsAssemblyException(Exception):
    pass


class Hit:
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


class Segment:
    def __init__(self, begin, end):
        self.begin = begin
        self.end = end


def unite_segments(segments):
    segments.sort(key=lambda segment: segment.begin)
    united_segments = [segments[0]]

    for i in xrange(1, len(segments)):
        if segments[i].begin <= united_segments[-1].end:
            if segments[i].end > united_segments[-1].end:
                united_segments[-1].end = segments[i].end
        else:
            united_segments.append(segments[i])

    return united_segments


def calc_alignment_rate(hit, mapping_segments):
    query_coverage = 0
    united_segments = unite_segments(mapping_segments)

    for segment in united_segments:
        query_coverage += segment.end - segment.begin

    return round(float(query_coverage) / hit.query_length, 5)


def calc_alignment_rates(paf_alignment):
    hits = []

    with open(paf_alignment) as f:
        for raw_hit in f:
            hits.append(Hit(raw_hit))

    hits.sort(key=lambda hit: (hit.query, hit.target))

    alignment_rates = dict()
    current_hit = None
    mapping_segments = []

    for hit in hits:
        if current_hit is None or current_hit.query != hit.query \
                               or current_hit.target != hit.target:
            if current_hit is not None:
                aln_rate = calc_alignment_rate(current_hit, mapping_segments)

                if current_hit.query not in alignment_rates:
                    alignment_rates[current_hit.query] = dict()

                alignment_rates[current_hit.query][current_hit.target] = aln_rate

            current_hit = hit
            mapping_segments = []

        mapping_segments.append(Segment(hit.query_start, hit.query_end))

    return alignment_rates


def find_unmapped_reads(alignment_rates, reads_files, aln_rate_threshold):
    unmapped_reads = dict()

    for file in reads_files:
        fasta_dict = fp.read_fasta_dict(file)

        for read, sequence in fasta_dict.items():
            contigs = alignment_rates.get(read)

            if contigs is None:
                unmapped_reads[read] = sequence
            else:
                is_unmapped = True
                for contig, aln_rate in contigs.items():
                    if aln_rate >= aln_rate_threshold:
                        is_unmapped = False
                        
                if is_unmapped:
                    unmapped_reads[read] = sequence


    return unmapped_reads


def run_minimap(preset, contigs_file, reads_files, num_proc, out_file):
    cmdline = [MINIMAP_BIN]
    cmdline.extend(['-x', preset])
    cmdline.extend([contigs_file])
    cmdline.extend(reads_files)
    cmdline.extend(['-t', str(num_proc)])

    try:
        devnull = open(os.devnull, "w")
        logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline, stderr=devnull,
                              stdout=open(out_file, "w"))
    except (subprocess.CalledProcessError, OSError) as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise ShortPlasmidsAssemblyException(str(e))

