import os
import logging
import subprocess


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

    return mapping_rates
