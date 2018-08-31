import os
import logging
import subprocess

import flye.utils.fasta_parser as fp
import flye.short_plasmids.unmapped_reads as unmapped
import flye.short_plasmids.circular_sequences as circular


logger = logging.getLogger()
MINIMAP_BIN = "flye-minimap2"


class ShortPlasmidsAssemblyException(Exception):
    pass


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
    _run_minimap(preset, contigs_path, args.reads,
                 args.threads, reads2contigs_mapping)

    logger.debug("Extracting unmapped reads")
    unmapped_reads, n_processed_reads = \
        unmapped.extract_unmapped_reads(args, reads2contigs_mapping,
                                        mapping_rate_threshold=0.5)

    n_unmapped_reads = len(unmapped_reads)
    unmapped_reads_ratio = 100 * float(len(unmapped_reads)) / n_processed_reads
    unmapped_reads_ratio = round(unmapped_reads_ratio, 1)
    logger.debug("Extracted {} unmapped reads ({} %)".format(
        n_unmapped_reads, unmapped_reads_ratio))

    unmapped_reads_path = os.path.join(work_dir, "unmapped_reads.fasta")
    fp.write_fasta_dict(unmapped_reads, unmapped_reads_path)

    unmapped_reads_mapping = os.path.join(work_dir, "unmapped_ava.paf")

    logger.debug("Finding self-mappings for unmapped reads")
    preset = ["ava-pb", "ava-ont"][args.platform == "nano"]
    _run_minimap(preset, unmapped_reads_path, [unmapped_reads_path],
                 args.threads, unmapped_reads_mapping)

    logger.debug("Extracting circular reads")
    circular_reads = circular.extract_circular_reads(unmapped_reads_mapping)
    logger.debug("Extracted {} circular reads".format(len(circular_reads)))

    logger.debug("Extracing circular pairs")
    circular_pairs = circular.extract_circular_pairs(unmapped_reads_mapping)
    logger.debug("Extracted {} circular pairs".format(len(circular_pairs)))

    logger.debug("Extracting unique plasmids from circular sequences")
    trimmed_circular_reads = \
        circular.trim_circular_reads(circular_reads, unmapped_reads)
    trimmed_circular_pairs = \
        circular.trim_circular_pairs(circular_pairs, unmapped_reads)
    trimmed_sequences_path = os.path.join(work_dir, "trimmed_sequences.fasta")

    fp.write_fasta_dict(trimmed_circular_reads, trimmed_sequences_path)
    fp.write_fasta_dict(trimmed_circular_pairs, trimmed_sequences_path, "a")

    trimmed_sequences_mapping = os.path.join(work_dir, "trimmed.paf")

    preset = ["ava-pb", "ava-ont"][args.platform == "nano"]
    _run_minimap(preset, trimmed_sequences_path, [trimmed_sequences_path],
                 args.threads, trimmed_sequences_mapping)

    plasmids = \
        circular.extract_unique_plasmids(trimmed_sequences_mapping,
                                         trimmed_sequences_path)
    logger.debug("Extracted {} unique plasmids".format(len(plasmids)))

    # remove all unnecesarry files
    os.remove(reads2contigs_mapping)
    os.remove(unmapped_reads_path)
    os.remove(unmapped_reads_mapping)
    os.remove(trimmed_sequences_path)
    os.remove(trimmed_sequences_mapping)

    return plasmids
