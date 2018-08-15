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

    def query_hit_length(self):
        return self.query_end - self.query_start

    def target_hit_length(self):
        return self.target_end - self.target_start

    def query_left_overhang(self):
        return self.query_start

    def query_right_overhang(self):
        return self.query_length - self.query_end

    def target_left_overhang(self):
        return self.target_start

    def target_right_overhang(self):
        return self.target_length - self.target_end


class Segment:
    def __init__(self, begin, end):
        self.begin = begin
        self.end = end

    def length(self):
        return self.end - self.begin


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


def calc_alignment_rate(sequence_length, mapping_segments):
    sequence_coverage = 0
    united_segments = unite_segments(mapping_segments)

    for segment in united_segments:
        sequence_coverage += segment.length()

    return round(float(sequence_coverage) / sequence_length, 5)


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
                aln_rate = calc_alignment_rate(current_hit.query_length,
                                               mapping_segments)
                if current_hit.query not in alignment_rates:
                    alignment_rates[current_hit.query] = dict()

                alignment_rates[current_hit.query][current_hit.target] = aln_rate

            current_hit = hit
            mapping_segments = []

        mapping_segments.append(Segment(hit.query_start, hit.query_end))

    return alignment_rates


def represents_circular_read(hit):
    if hit.query != hit.target:
        return False

    if not hit.query_start < hit.query_end < hit.target_start < hit.target_end:
        return False

    max_overhang = 150
    if not hit.query_left_overhang() < max_overhang:
        return False

    if not hit.target_right_overhang() < max_overhang:
        return False

    return True


def find_unmapped_reads(alignment_rates, reads_files, aln_rate_threshold):
    unmapped_reads = dict()
    n_reads = 0

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

        n_reads += len(fasta_dict)

    return unmapped_reads, n_reads


def find_circular_reads(paf_unmapped_reads):
    circular_reads = dict()

    with open(paf_unmapped_reads) as f:
        for raw_hit in f:
            current_hit = Hit(raw_hit)
            if represents_circular_read(current_hit):
                hit = circular_reads.get(current_hit.query)
                if hit is None or current_hit.query_hit_length() > \
                        hit.query_hit_length():
                    circular_reads[current_hit.query] = current_hit

    return circular_reads


def trim_circular_reads(circular_reads, unmapped_reads):
    trimmed_reads = dict()

    ctr = 0
    for read, hit in circular_reads.items():
        sequence = unmapped_reads[read]
        trimmed_reads['circular_seq' + str(ctr)] = sequence[:hit.target_start]
        ctr += 1

    return trimmed_reads


def find_circular_pairs(paf_unmapped_reads):
    hits = []

    with open(paf_unmapped_reads) as f:
        for raw_hit in f:
            hits.append(Hit(raw_hit))

    hits.sort(key=lambda hit: (hit.query, hit.target))

    circular_pairs = []
    circular_pair = [None, None]
    previous_hit = None
    has_overlap = False
    is_circular = False
    max_overhang = 300

    for hit in hits:
        if hit.query == hit.target:
            continue

        if previous_hit is None or previous_hit.query != hit.query or \
                previous_hit.target != hit.target:
            if previous_hit is not None and has_overlap and is_circular:
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


def trim_circular_pairs(circular_pairs, unmapped_reads):
    trimmed_reads = dict()

    for i, pair in enumerate(circular_pairs):
        left_sequence = unmapped_reads[pair[0].query]
        right_sequence = unmapped_reads[pair[0].target]
        trimmed_sequence = left_sequence[pair[1].query_end:pair[0].query_end]
        trimmed_sequence += right_sequence[pair[0].target_end:]
        trimmed_reads['circular_pair' + str(i)] = trimmed_sequence

    return trimmed_reads


def find_connected_components(graph):
    def dfs(vertex, connected_components_counter):
        connected_components[vertex] = connected_components_counter
        used[vertex] = True
        for neighbour in graph[vertex]:
            if not used[neighbour]:
                dfs(neighbour, connected_components_counter)

    n_vertices = len(graph)
    connected_components = [0 for _ in range(n_vertices)]
    connected_components_counter = 0
    used = [0 for _ in range(n_vertices)]

    for i in range(n_vertices):
        if not used[i]:
            dfs(i, connected_components_counter)
            connected_components_counter += 1

    return connected_components, connected_components_counter


def extract_unique_plasmids(paf_trimmed_reads, trimmed_reads_fasta):
    trimmed_reads = set()
    hits = []

    with open(paf_trimmed_reads) as f:
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

        if current_hit is None or hit.query != current_hit.query or \
                hit.target != current_hit.target:
            if current_hit is not None:
                query_length = current_hit.query_length
                target_length = current_hit.target_length
                query_aln_rate = calc_alignment_rate(query_length,
                                                     query_mapping_segments)
                target_aln_rate = calc_alignment_rate(target_length,
                                                      target_mapping_segments)

                if query_aln_rate >= 0.8 and target_aln_rate >= 0.8 and \
                        abs(query_length - target_length) <= 500:
                    vertex1 = read2int[current_hit.query]
                    vertex2 = read2int[current_hit.target]
                    similarity_graph[vertex1].append(vertex2)
                    similarity_graph[vertex2].append(vertex1)

            query_mapping_segments = []
            target_mapping_segments = []
            current_hit = hit

        query_mapping_segments.append(Segment(hit.query_start, hit.query_end))
        target_mapping_segments.append(Segment(hit.target_start, hit.target_end))

    connected_components, n_components = find_connected_components(similarity_graph)

    groups = [[] for _ in range(n_components)]
    for i in range(len(connected_components)):
        groups[connected_components[i]].append(int2read[i])

    groups = [group for group in groups if len(group) > 1]
    trimmed_reads_dict = fp.read_fasta_dict(trimmed_reads_fasta)
    unique_plasmids = dict()

    for group in groups:
        unique_plasmids[group[0]] = trimmed_reads_dict[group[0]]

    return unique_plasmids


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
