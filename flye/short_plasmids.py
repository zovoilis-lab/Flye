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


def calc_alignment_rates(paf_alignment, paf_alignment_sorted):
    hits = []

    with open(paf_alignment) as f:
        for raw_hit in f:
            hits.append(Hit(raw_hit))

    hits.sort(key=lambda hit: (hit.query, hit.target))

    with open(paf_alignment_sorted, 'w') as f:
        for hit in hits:
            f.write(hit.query + ' ' + hit.target + '\n')



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

