import os
import logging
import subprocess
import flye.utils.fasta_parser as fp

from flye.short_plasmids.utils import find_connected_components


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


class MappingSegment:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def length(self):
        return self.end - self.start + 1


def _unite_mapping_segments(segments):
    segments.sort(key=lambda segment: segment.start)
    united_segments = [segments[0]]

    for i in xrange(1, len(segments)):
        if segments[i].start <= united_segments[-1].end:
            if segments[i].end > united_segments[-1].end:
                united_segments[-1].end = segments[i].end
        else:
            united_segments.append(segments[i])

    return united_segments


def _calc_mapping_rate(read_length, mapping_segments):
    read_coverage = 0
    united_segments = _unite_mapping_segments(mapping_segments)

    for mapping_segment in united_segments:
        read_coverage += mapping_segment.length()

    return round(float(read_coverage) / read_length, 3)


def _calc_mapping_rates(reads2contigs_mapping):
    hits = []

    with open(reads2contigs_mapping) as f:
        for raw_hit in f:
            hits.append(Hit(raw_hit))

    hits.sort(key=lambda hit: (hit.query, hit.target))

    mapping_rates = dict()
    current_hit = None
    mapping_segments = []

    for hit in hits:
        if current_hit is None or hit.query != current_hit.query or \
                                  hit.target != current_hit.target:
            if current_hit is not None:
                mapping_rate = _calc_mapping_rate(current_hit.query_length,
                                                  mapping_segments)
                if current_hit.query not in mapping_rates:
                    mapping_rates[current_hit.query] = dict()

                mapping_rates[current_hit.query][current_hit.target] = mapping_rate

            current_hit = hit
            mapping_segments = []

        mapping_segments.append(MappingSegment(hit.query_start, hit.query_end))

    return mapping_rates


def _extract_unmapped_reads(mapping_rates, reads_paths, mapping_rate_threshold):
    unmapped_reads = dict()
    n_processed_reads = 0

    for file in reads_paths:
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


def _is_circular_read(hit, max_overhang=150):
    if hit.query != hit.target:
        return False

    if not hit.query_start < hit.query_end < hit.target_start < hit.target_end:
        return False

    if not hit.query_left_overhang() < max_overhang:
        return False

    if not hit.target_right_overhang() < max_overhang:
        return False

    return True


def _extract_circular_reads(unmapped_reads_mapping, max_overhang=150):
    circular_reads = dict()

    with open(unmapped_reads_mapping) as f:
        for raw_hit in f:
            current_hit = Hit(raw_hit)
            if _is_circular_read(current_hit, max_overhang):
                hit = circular_reads.get(current_hit.query)
                if hit is None or current_hit.query_mapping_length() > \
                   hit.query_mapping_length():
                    circular_reads[current_hit.query] = current_hit

    return circular_reads


def _trim_circular_reads(circular_reads, unmapped_reads):
    trimmed_circular_reads = dict()

    i = 0
    for read, hit in circular_reads.items():
        sequence = unmapped_reads[read][:hit.target_start].upper()
        trimmed_circular_reads["circular_read" + str(i)] = sequence
        i += 1

    return trimmed_circular_reads


def _mapping_segments_without_intersection(circular_pair):
    if not circular_pair[1].query_start < circular_pair[1].query_end < \
           circular_pair[0].query_start < circular_pair[0].query_end:
        return False

    if not circular_pair[0].target_start < circular_pair[0].target_end < \
           circular_pair[1].target_start < circular_pair[1].target_end:
        return False

    return True


def _extract_circular_pairs(unmapped_reads_mapping, max_overhang=300):
    hits = []

    with open(unmapped_reads_mapping) as f:
        for raw_hit in f:
            hits.append(Hit(raw_hit))

    hits.sort(key=lambda hit: (hit.query, hit.target))

    circular_pairs = []
    circular_pair = [None, None]
    previous_hit = None
    has_overlap = False
    is_circular = False

    for hit in hits:
        if hit.query == hit.target:
            continue

        if previous_hit is None or \
           hit.query != previous_hit.query or \
           hit.target != previous_hit.target:
            if previous_hit is not None and has_overlap and is_circular:
                if _mapping_segments_without_intersection(circular_pair):
                    circular_pairs.append(circular_pair)

            circular_pair = [None, None]
            has_overlap = False
            is_circular = False
            previous_hit = hit

        if not has_overlap:
            if hit.query_right_overhang() < max_overhang and \
               hit.target_left_overhang() < max_overhang:
                has_overlap = True
                circular_pair[0] = hit
                continue

        if not is_circular:
            if hit.query_left_overhang() < max_overhang and \
               hit.target_right_overhang() < max_overhang:
                is_circular = True
                circular_pair[1] = hit

    return circular_pairs


def _trim_circular_pairs(circular_pairs, unmapped_reads):
    trimmed_circular_pairs = dict()

    for i, pair in enumerate(circular_pairs):
        lhs = unmapped_reads[pair[0].query]
        rhs = unmapped_reads[pair[0].target]
        trimmed_seq = lhs[pair[1].query_end:pair[0].query_end]
        trimmed_seq += rhs[pair[0].target_end:]
        trimmed_circular_pairs["circular_pair" + str(i)] = trimmed_seq.upper()

    return trimmed_circular_pairs


def _extract_unique_plasmids(trimmed_reads_mappings, trimmed_reads_path,
                             mapping_rate_threashold=0.8,
                             max_length_difference=500,
                             min_sequence_length=1000):
    trimmed_reads = set()
    hits = []

    with open(trimmed_reads_mappings) as f:
        for raw_hit in f:
            hit = Hit(raw_hit)
            hits.append(hit)
            trimmed_reads.add(hit.query)
            trimmed_reads.add(hit.target)

    trimmed_reads = list(trimmed_reads)
    n_trimmed_reads = len(trimmed_reads)
    read2int = dict()
    int2read = dict()

    for i in range(n_trimmed_reads):
        read2int[trimmed_reads[i]] = i
        int2read[i] = trimmed_reads[i]

    similarity_graph = [[] for _ in range(n_trimmed_reads)]
    hits.sort(key=lambda hit: (hit.query, hit.target))

    current_hit = None
    query_mapping_segments = []
    target_mapping_segments = []

    for hit in hits:
        if hit.query == hit.target:
            continue

        if current_hit is None or \
           hit.query != current_hit.query or \
           hit.target != current_hit.target:
            if current_hit is not None:
                query_length = current_hit.query_length
                target_length = current_hit.target_length
                query_mapping_rate = \
                    _calc_mapping_rate(query_length, query_mapping_segments)
                target_mapping_rate = \
                    _calc_mapping_rate(target_length, target_mapping_segments)

                if query_mapping_rate > mapping_rate_threashold and \
                   target_mapping_rate > mapping_rate_threashold and \
                   abs(query_length - target_length) < max_length_difference:
                    vertex1 = read2int[current_hit.query]
                    vertex2 = read2int[current_hit.target]
                    similarity_graph[vertex1].append(vertex2)
                    similarity_graph[vertex2].append(vertex1)

            query_mapping_segments = []
            target_mapping_segments = []
            current_hit = hit

        query_mapping_segments.append(
            MappingSegment(hit.query_start, hit.query_end))
        target_mapping_segments.append(
            MappingSegment(hit.target_start, hit.target_end))

    connected_components, n_components = \
        find_connected_components(similarity_graph)

    groups = [[] for _ in range(n_components)]
    for i in range(len(connected_components)):
        groups[connected_components[i]].append(int2read[i])

    groups = [group for group in groups if len(group) > 1]
    trimmed_reads_dict = fp.read_sequence_dict(trimmed_reads_path)
    unique_plasmids = dict()

    for group in groups:
        sequence = trimmed_reads_dict[group[0]]
        if len(sequence) >= min_sequence_length:
            unique_plasmids[group[0]] = sequence

    return unique_plasmids


def _run_minimap(preset, contigs_path, reads_paths, num_proc, out_file):
    cmdline = [MINIMAP_BIN]
    cmdline.extend(["-x", preset])
    cmdline.extend([contigs_path])
    cmdline.extend(reads_paths)
    cmdline.extend(["-t", str(num_proc)])

    try:
        devnull = open(os.devnull, "w")
        logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline, stderr=devnull,
                              stdout=open(out_file, "w"))
    except (subprocess.CalledProcessError, OSError) as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise ShortPlasmidsAssemblyException(str(e))


def assemble_short_plasmids(args, work_dir, contigs_path):
    logger.debug("Assembling short plasmids")

    reads2contigs_mapping = os.path.join(work_dir, "reads2contigs.paf")
    preset = ["map-pb", "map-ont"][args.platform == "nano"]

    if not os.path.isfile(reads2contigs_mapping):
        logger.debug("Mapping reads to contigs")
        _run_minimap(preset, contigs_path, args.reads,
                     args.threads, reads2contigs_mapping)

    logger.debug("Calculating mapping rates")
    mapping_rates = _calc_mapping_rates(reads2contigs_mapping)

    logger.debug("Extracting unmapped reads")
    unmapped_reads, n_processed_reads = \
        _extract_unmapped_reads(mapping_rates, args.reads,
                                mapping_rate_threshold=0.5)

    n_unmapped_reads = len(unmapped_reads)
    unmapped_reads_ratio = 100 * float(len(unmapped_reads)) / n_processed_reads
    unmapped_reads_ratio = round(unmapped_reads_ratio, 1)
    logger.debug("Extracted {} unmapped reads ({} %)".format(
        n_unmapped_reads, unmapped_reads_ratio))

    unmapped_reads_path = os.path.join(work_dir, "unmapped_reads.fasta")
    fp.write_fasta_dict(unmapped_reads, unmapped_reads_path)

    unmapped_reads_mapping = os.path.join(work_dir, "unmapped_ava.paf")

    if not os.path.isfile(unmapped_reads_mapping):
        logger.debug("Finding self-mappings for unmapped reads")
        preset = ["ava-pb", "ava-ont"][args.platform == "nano"]
        _run_minimap(preset, unmapped_reads_path, [unmapped_reads_path],
                     args.threads, unmapped_reads_mapping)

    logger.debug("Extracting circular reads")
    circular_reads = _extract_circular_reads(unmapped_reads_mapping)
    logger.debug("Extracted {} circular reads".format(len(circular_reads)))

    logger.debug("Extracing circular pairs")
    circular_pairs = _extract_circular_pairs(unmapped_reads_mapping)
    logger.debug("Extracted {} circular pairs".format(len(circular_pairs)))

    logger.debug("Extracting unique plasmids from circular sequences")
    trimmed_circular_reads = _trim_circular_reads(circular_reads, unmapped_reads)
    trimmed_circular_pairs = _trim_circular_pairs(circular_pairs, unmapped_reads)
    trimmed_sequences_path = os.path.join(work_dir, "trimmed_sequences.fasta")

    fp.write_fasta_dict(trimmed_circular_reads, trimmed_sequences_path)
    fp.write_fasta_dict(trimmed_circular_pairs, trimmed_sequences_path, "a")

    trimmed_sequences_mappings = os.path.join(work_dir, "trimmed.paf")

    if not os.path.isfile(trimmed_sequences_mappings):
        preset = ["ava-pb", "ava-ont"][args.platform == "nano"]
        _run_minimap(preset, trimmed_sequences_path, [trimmed_sequences_path],
                     args.threads, trimmed_sequences_mappings)

    plasmids = _extract_unique_plasmids(trimmed_sequences_mappings,
                                        trimmed_sequences_path)
    logger.debug("Extracted {} unique plasmids".format(len(plasmids)))
    return plasmids
