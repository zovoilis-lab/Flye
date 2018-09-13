#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)


import flye.utils.fasta_parser as fp
from flye.polishing.alignment import read_paf


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
    hits = read_paf(reads2contigs_mapping)
    hits.sort(key=lambda hit: (hit.query, hit.target))

    mapping_rates = dict()
    current_hit = None
    mapping_segments = []

    for hit in hits:
        if current_hit is None or hit.query != current_hit.query or \
                                  hit.target != current_hit.target:
            if current_hit is not None:
                mapping_rate = calc_mapping_rate(current_hit.query_length,
                                                 mapping_segments)
                if current_hit.query not in mapping_rates:
                    mapping_rates[current_hit.query] = dict()

                mapping_rates[current_hit.query][current_hit.target] = mapping_rate

            current_hit = hit
            mapping_segments = []

        mapping_segments.append(MappingSegment(hit.query_start, hit.query_end))

    return mapping_rates


def extract_unmapped_reads(args, reads2contigs_mapping, mapping_rate_threshold):
    mapping_rates = calc_mapping_rates(reads2contigs_mapping)
    unmapped_reads = dict()
    n_processed_reads = 0

    for file in args.reads:
        fasta_dict = fp.read_sequence_dict(file)
        for read, sequence in fasta_dict.items():
            contigs = mapping_rates.get(read)
            if contigs is None:
                unmapped_reads[read] = sequence
            else:
                is_unmapped = True
                for contig, mapping_rate in contigs.items():
                    if mapping_rate >= mapping_rate_threshold:
                        is_unmapped = False

                if is_unmapped:
                    unmapped_reads[read] = sequence

        n_processed_reads += len(fasta_dict)

    return unmapped_reads, n_processed_reads
