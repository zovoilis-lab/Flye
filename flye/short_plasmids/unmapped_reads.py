#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)


import flye.utils.fasta_parser as fp
from flye.polishing.alignment import read_paf_grouped
import logging
from collections import defaultdict

logger = logging.getLogger()


class MappingSegment:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def length(self):
        return self.end - self.start + 1


def unite_mapping_segments(segments):
    segments.sort(key=lambda segment: segment.start)
    united_segments = [segments[0]]

    for i in xrange(1, len(segments)):
        if segments[i].start <= united_segments[-1].end:
            if segments[i].end > united_segments[-1].end:
                united_segments[-1].end = segments[i].end
        else:
            united_segments.append(segments[i])

    return united_segments


def calc_mapping_rate(read_length, mapping_segments):
    read_coverage = 0
    united_segments = unite_mapping_segments(mapping_segments)

    for mapping_segment in united_segments:
        read_coverage += mapping_segment.length()

    return round(float(read_coverage) / read_length, 3)


def calc_mapping_rates(reads2contigs_mapping):
    mapping_rates = defaultdict(dict)
    target_hits = defaultdict(list)
    prev_hit = None

    #read_paf_grouped assumes that hits are sorted in query order.
    #It returns chunks of alignments for all (query, target) combinations
    for hit_group in read_paf_grouped(reads2contigs_mapping):
        map_segs = map(lambda h: MappingSegment(h.query_start, h.query_end),
                       hit_group)
        mapping_rate = calc_mapping_rate(hit_group[0].query_length, map_segs)
        mapping_rates[hit_group[0].query][hit_group[0].target] = mapping_rate

    return mapping_rates


def extract_unmapped_reads(args, reads2contigs_mapping, unmapped_reads_path,
                           mapping_rate_threshold):
    mapping_rates = calc_mapping_rates(reads2contigs_mapping)
    total_bases = 0
    unmapped_bases = 0

    with open(unmapped_reads_path, "w") as fout:
        for file in args.reads:
            for hdr, sequence in fp.stream_sequence(file):
                total_bases += len(sequence)

                is_unmapped = True
                contigs = mapping_rates.get(hdr)
                if contigs is not None:
                    is_unmapped = True
                    for contig, mapping_rate in contigs.iteritems():
                        if mapping_rate >= mapping_rate_threshold:
                            is_unmapped = False

                if is_unmapped:
                    unmapped_bases += len(sequence)
                    fout.write(">{0}\n{1}\n".format(hdr, sequence))

    logger.debug("Unmapped sequence: {0} / {1} ({2})"
                 .format(unmapped_bases, total_bases,
                         float(unmapped_bases) / total_bases))
