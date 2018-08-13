import os
import logging
import subprocess


logger = logging.getLogger()
MINIMAP_BIN = "flye-minimap2"


class ShortPlasmidsAssemblyException(Exception):
    pass


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