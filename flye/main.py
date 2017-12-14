#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Main logic of the package
"""

from __future__ import print_function
import sys
import os
import logging
import argparse
import json
import shutil
import subprocess

import flye.alignment as aln
import flye.bubbles as bbl
import flye.polish as pol
import flye.fasta_parser as fp
import flye.assemble as asm
import flye.repeat_graph as repeat
import flye.consensus as cons
import flye.scaffolder as scf
from flye.__version__ import __version__
import flye.config as config
from flye.bytes2human import human2bytes

logger = logging.getLogger()

class ResumeException(Exception):
    pass

class Job(object):
    """
    Describes an abstract list of jobs with persistent
    status that can be resumed
    """
    run_description = {"stage_name" : ""}

    def __init__(self):
        self.name = None
        self.args = None
        self.work_dir = None
        self.out_files = {}
        self.log_file = None

    def run(self):
        pass

    def save(self, save_file):
        Job.run_description["stage_name"] = self.name

        with open(save_file, "w") as fp:
            json.dump(Job.run_description, fp)

    def load(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)
            Job.run_description = data

    def completed(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)

            for file in self.out_files.values():
                if not os.path.exists(file):
                    return False

            return True


class JobAssembly(Job):
    def __init__(self, args, work_dir, log_file):
        super(JobAssembly, self).__init__()
        #self.out_assembly = out_assembly
        self.args = args
        self.work_dir = work_dir
        self.log_file = log_file

        self.name = "assembly"
        self.assembly_dir = os.path.join(self.work_dir, "0-assembly")
        self.assembly_filename = os.path.join(self.assembly_dir,
                                              "draft_assembly.fasta")
        self.out_files["assembly"] = self.assembly_filename

    def run(self):
        if not os.path.isdir(self.assembly_dir):
            os.mkdir(self.assembly_dir)
        asm.assemble(self.args, self.assembly_filename, self.log_file,
                     self.args.asm_config)
        if os.path.getsize(self.assembly_filename) == 0:
            raise asm.AssembleException("No contigs were assembled - "
                                        "are you using corrected input instead of raw?")


class JobRepeat(Job):
    def __init__(self, args, work_dir, log_file, in_assembly):
        super(JobRepeat, self).__init__()

        self.args = args
        self.in_assembly = in_assembly
        self.log_file = log_file
        self.name = "repeat"

        self.repeat_dir = os.path.join(work_dir, "2-repeat")
        contig_sequences = os.path.join(self.repeat_dir, "graph_paths.fasta")
        assembly_graph = os.path.join(self.repeat_dir, "graph_final.dot")
        contigs_stats = os.path.join(self.repeat_dir, "contigs_stats.txt")
        self.out_files["contigs"] = contig_sequences
        self.out_files["scaffold_links"] = os.path.join(self.repeat_dir,
                                                        "scaffolds_links.txt")
        self.out_files["assembly_graph"] = assembly_graph
        self.out_files["stats"] = contigs_stats

    def run(self):
        if not os.path.isdir(self.repeat_dir):
            os.mkdir(self.repeat_dir)
        logger.info("Performing repeat analysis")
        repeat.analyse_repeats(self.args, self.in_assembly, self.repeat_dir,
                               self.log_file, self.args.asm_config)


class JobFinalize(Job):
    def __init__(self, args, work_dir, log_file,
                 contigs_file, graph_file, repeat_stats,
                 polished_stats, scaffold_links):
        super(JobFinalize, self).__init__()

        self.args = args
        self.log_file = log_file
        self.name = "finalize"
        self.contigs_file = contigs_file
        self.graph_file = graph_file
        self.repeat_stats = repeat_stats
        self.polished_stats = polished_stats
        self.scaffold_links = scaffold_links

        self.out_files["contigs"] = os.path.join(work_dir, "contigs.fasta")
        self.out_files["scaffolds"] = os.path.join(work_dir, "scaffolds.fasta")
        self.out_files["stats"] = os.path.join(work_dir, "contigs_info.txt")
        self.out_files["graph"] = os.path.join(work_dir, "assembly_graph.dot")

    def run(self):
        shutil.copy2(self.contigs_file, self.out_files["contigs"])
        shutil.copy2(self.graph_file, self.out_files["graph"])

        scf.generate_scaffolds(self.contigs_file, self.scaffold_links,
                               self.out_files["scaffolds"])

        contigs_length = {}
        contigs_coverage = {}
        repeat_stats = {}
        table_lines = open(self.repeat_stats, "r").readlines()
        header_line = table_lines[0]
        for line in table_lines[1:]:
            tokens = line.strip().split("\t")
            repeat_stats[tokens[0]] = tokens[1:]
            if self.polished_stats is None:
                contigs_length[tokens[0]] = int(tokens[1])
                contigs_coverage[tokens[0]] = int(tokens[2])

        if self.polished_stats is not None:
            for line in open(self.polished_stats, "r").readlines()[1:]:
                tokens = line.strip().split("\t")
                contigs_length[tokens[0]] = int(tokens[1])
                contigs_coverage[tokens[0]] = int(tokens[2])

        all_contigs = sorted(contigs_length.keys(),
                             key=lambda c: int(c.split("_")[1]))

        with open(self.out_files["stats"], "w") as f:
            f.write(header_line)
            for ctg_id in all_contigs:
                fields = repeat_stats[ctg_id]
                fields[0] = str(contigs_length[ctg_id])
                fields[1] = str(contigs_coverage[ctg_id])
                f.write("{0}\t{1}\n".format(ctg_id, "\t".join(fields)))

        total_length = sum(contigs_length.values())
        if total_length > 0:
            largest_len = contigs_length[max(contigs_length,
                                             key=contigs_length.get)]
            ctg_n50 = _calc_n50(contigs_length.values(), total_length)
            mean_read_cov = 0
            for ctg_id in contigs_coverage:
                mean_read_cov += contigs_length[ctg_id] * contigs_coverage[ctg_id]
            mean_read_cov /= total_length

            logger.info("Assembly statistics:\n\n"
                        "\tContigs:\t{0}\n"
                        "\tTotal length:\t{1}\n"
                        "\tContigs N50:\t{2}\n"
                        "\tLargest contig:\t{3}\n"
                        "\tMean coverage:\t{4}\n"
                        .format(len(contigs_length), total_length,
                                ctg_n50, largest_len, mean_read_cov))

        logger.info("Done! your assembly is in {0} file"
                    .format(self.out_files["contigs"]))


class JobConsensus(Job):
    def __init__(self, args, work_dir, in_contigs):
        super(JobConsensus, self).__init__()

        self.args = args
        self.in_contigs = in_contigs
        self.consensus_dir = os.path.join(work_dir, "1-consensus")
        self.out_consensus = os.path.join(self.consensus_dir, "consensus.fasta")
        self.name = "consensus"
        self.out_files["consensus"] = self.out_consensus

    def run(self):
        if not os.path.isdir(self.consensus_dir):
            os.mkdir(self.consensus_dir)

        logger.info("Running Minimap2")
        out_alignment = os.path.join(self.consensus_dir, "minimap.sam")
        aln.make_alignment(self.in_contigs, self.args.reads, self.args.threads,
                           self.consensus_dir, self.args.platform, out_alignment)

        contigs_info = aln.get_contigs_info(self.in_contigs)
        logger.info("Computing consensus")
        consensus_fasta = cons.get_consensus(out_alignment, self.in_contigs,
                                             contigs_info, self.args.threads,
                                             self.args.platform,
                                             self.args.min_overlap)
        fp.write_fasta_dict(consensus_fasta, self.out_consensus)


class JobPolishing(Job):
    def __init__(self, args, work_dir, log_file, in_contigs):
        super(JobPolishing, self).__init__()

        self.args = args
        self.log_file = log_file
        self.in_contigs = in_contigs
        self.polishing_dir = os.path.join(work_dir, "3-polishing")

        self.name = "polishing"
        final_conitgs = os.path.join(self.polishing_dir,
                                     "polished_{0}.fasta".format(args.num_iters))
        self.out_files["contigs"] = final_conitgs
        self.out_files["stats"] = os.path.join(self.polishing_dir,
                                               "contigs_stats.txt")

    def run(self):
        if not os.path.isdir(self.polishing_dir):
            os.mkdir(self.polishing_dir)

        prev_assembly = self.in_contigs
        contig_lengths = None
        for i in xrange(self.args.num_iters):
            logger.info("Polishing genome ({0}/{1})".format(i + 1,
                                                self.args.num_iters))

            alignment_file = os.path.join(self.polishing_dir,
                                          "minimap_{0}.sam".format(i + 1))
            logger.info("Running Minimap2")
            aln.make_alignment(prev_assembly, self.args.reads, self.args.threads,
                               self.polishing_dir, self.args.platform,
                               alignment_file)

            logger.info("Separating alignment into bubbles")
            contigs_info = aln.get_contigs_info(prev_assembly)
            bubbles_file = os.path.join(self.polishing_dir,
                                        "bubbles_{0}.fasta".format(i + 1))
            coverage_stats = \
                bbl.make_bubbles(alignment_file, contigs_info, prev_assembly,
                                 self.args.platform, self.args.threads,
                                 self.args.min_overlap, bubbles_file)

            logger.info("Correcting bubbles")
            polished_file = os.path.join(self.polishing_dir,
                                         "polished_{0}.fasta".format(i + 1))
            contig_lengths = pol.polish(bubbles_file, self.args.threads,
                                        self.args.platform, self.polishing_dir,
                                        i + 1, polished_file)
            prev_assembly = polished_file

        with open(self.out_files["stats"], "w") as f:
            f.write("contig_id\tlength\tcoverage\n")
            for ctg_id in contig_lengths:
                f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                        contig_lengths[ctg_id], coverage_stats[ctg_id]))


def _create_job_list(args, work_dir, log_file):
    """
    Build pipeline as a list of consecutive jobs
    """
    jobs = []

    #Assembly job
    jobs.append(JobAssembly(args, work_dir, log_file))
    draft_assembly = jobs[-1].out_files["assembly"]

    #Consensus
    if args.read_type == "raw":
        jobs.append(JobConsensus(args, work_dir, draft_assembly))
        draft_assembly = jobs[-1].out_files["consensus"]

    #Repeat analysis
    jobs.append(JobRepeat(args, work_dir, log_file, draft_assembly))
    raw_contigs = jobs[-1].out_files["contigs"]
    scaffold_links = jobs[-1].out_files["scaffold_links"]
    graph_file = jobs[-1].out_files["assembly_graph"]
    repeat_stats = jobs[-1].out_files["stats"]

    #Polishing
    contigs_file = raw_contigs
    polished_stats = None
    if args.num_iters > 0:
        jobs.append(JobPolishing(args, work_dir, log_file, raw_contigs))
        contigs_file = jobs[-1].out_files["contigs"]
        polished_stats = jobs[-1].out_files["stats"]

    #Report results
    jobs.append(JobFinalize(args, work_dir, log_file, contigs_file,
                            graph_file, repeat_stats, polished_stats,
                            scaffold_links))

    return jobs


def _set_kmer_size(args):
    """
    Select k-mer size based on the target genome size
    """
    if args.genome_size.isdigit():
        args.genome_size = int(args.genome_size)
    else:
        args.genome_size = human2bytes(args.genome_size.upper())

    """
    multiplier = 1
    suffix = args.reads.rsplit(".")[-1]
    if suffix == "gz":
        suffix = args.reads.rsplit(".")[-2]
        multiplier = 2

    if suffix in ["fasta", "fa"]:
        reads_size = os.path.getsize(args.reads) * multiplier
    elif suffix in ["fastq", "fq"]:
        reads_size = os.path.getsize(args.reads) * multiplier / 2
    else:
        raise ResumeException("Uknown input reads format: " + suffix)
    """

    #genome_coverage = reads_size / args.genome_size
    logger.debug("Genome size: {0}".format(args.genome_size))
    #logger.debug("Estimated genome coverage: {0}".format(genome_coverage))
    args.kmer_size = config.vals["small_kmer"]
    if args.genome_size > config.vals["big_genome"]:
        args.kmer_size = config.vals["big_kmer"]
    logger.debug("Chosen k-mer size: {0}".format(args.kmer_size))


def _set_read_attributes(args):
    root = os.path.dirname(__file__)
    if args.read_type == "raw":
        args.asm_config = os.path.join(root, "resource", config.vals["raw_cfg"])
    elif args.read_type in ["corrected", "subassemblies"]:
        args.asm_config = os.path.join(root, "resource",
                                       config.vals["corrected_cfg"])


def _run(args):
    """
    Runs the pipeline
    """
    logger.info("Running Flye " + _version())
    logger.debug("Cmd: {0}".format(" ".join(sys.argv)))

    for read_file in args.reads:
        if not os.path.exists(read_file):
            raise ResumeException("Can't open " + read_file)

    save_file = os.path.join(args.out_dir, "flye.save")
    jobs = _create_job_list(args, args.out_dir, args.log_file)

    current_job = 0
    if args.resume or args.resume_from:
        if not os.path.exists(save_file):
            raise ResumeException("Can't find save file")

        logger.info("Resuming previous run")
        if args.resume_from:
            job_to_resume = args.resume_from
        else:
            job_to_resume = json.load(open(save_file, "r"))["stage_name"]

        can_resume = False
        for i in xrange(len(jobs)):
            if jobs[i].name == job_to_resume:
                jobs[i].load(save_file)
                current_job = i
                if not jobs[i - 1].completed(save_file):
                    raise ResumeException("Can't resume: stage {0} incomplete"
                                          .format(jobs[i].name))
                can_resume = True
                break

        if not can_resume:
            raise ResumeException("Can't resume: stage {0} does not exist"
                                  .format(job_to_resume))

    for i in xrange(current_job, len(jobs)):
        jobs[i].save(save_file)
        jobs[i].run()


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len / 2:
            n50 = l
            break
    return n50


def _usage():
    return ("flye (--pacbio-raw | --pacbio-corr | --nano-raw |\n"
            "\t     --nano-corr | --subassemblies) file1 [file_2 ...]\n"
            "\t     --genome-size size --out-dir dir_path [--threads int]\n"
            "\t     [--iterations int] [--min-overlap int] [--resume]\n"
            "\t     [--debug] [--version] [--help]")


def _epilog():
    return ("Input reads could be in FASTA or FASTQ format, uncompressed\n"
            "or compressed with gz. Currenlty, raw and corrected reads\n"
            "from PacBio and ONT are supported. Additionally, --subassemblies\n"
            "option does a consensus assembly of high-quality input contigs.\n"
            "You may specify multiple fles with reads (separated by spaces).\n"
            "Mixing different read types is not yet supported.\n\n"
            "You must provide an estimate of the genome size as input,\n"
            "which is used for solid k-mers selection. The estimate could\n"
            "be rough (e.g. withing 0.5x-2x range) and does not affect\n"
            "the other assembly stages. Standard size modificators are\n"
            "supported (e.g. 5m or 2.6g)")


def _version():
    repo_root = os.path.dirname((os.path.dirname(__file__)))
    try:
        git_label = subprocess.check_output(["git", "-C", repo_root,
                                            "describe", "--tags"],
                                            stderr=open(os.devnull, "w"))
        commit_id = git_label.strip("\n").split("-", 1)[-1]
        return __version__ + "-" + commit_id
    except (subprocess.CalledProcessError, OSError):
        pass
    return __version__ + "-release"


def main():
    def check_int_range(value, min_val, max_val, require_odd=False):
        ival = int(value)
        if ival < min_val or ival > max_val:
             raise argparse.ArgumentTypeError("value should be in "
                            "range [{0}, {1}]".format(min_val, max_val))
        if require_odd and ival % 2 == 0:
            raise argparse.ArgumentTypeError("should be an odd number")
        return ival

    parser = argparse.ArgumentParser \
        (description="Assembly of long and error-prone reads",
         formatter_class=argparse.RawDescriptionHelpFormatter,
         usage=_usage(), epilog=_epilog())

    read_group = parser.add_mutually_exclusive_group(required=True)
    read_group.add_argument("--pacbio-raw", dest="pacbio_raw",
                        default=None, metavar="path", nargs="+",
                        help="PacBio raw reads")
    read_group.add_argument("--pacbio-corr", dest="pacbio_corrected",
                        default=None, metavar="path", nargs="+",
                        help="PacBio corrected reads")
    read_group.add_argument("--nano-raw", dest="nano_raw", nargs="+",
                        default=None, metavar="path",
                        help="ONT raw reads")
    read_group.add_argument("--nano-corr", dest="nano_corrected", nargs="+",
                        default=None, metavar="path",
                        help="ONT corrected reads")
    read_group.add_argument("--subassemblies", dest="subassemblies", nargs="+",
                        default=None, metavar="path",
                        help="high-quality contig-like input")

    parser.add_argument("-g", "--genome-size", dest="genome_size",
                        metavar="size", required=True,
                        help="estimated genome size (for example, 5m or 2.6g)")
    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")

    parser.add_argument("-t", "--threads", dest="threads",
                        type=lambda v: check_int_range(v, 1, 128),
                        default=1, metavar="int", help="number of parallel threads "
                        "(default: 1)")
    parser.add_argument("-i", "--iterations", dest="num_iters",
                        type=lambda v: check_int_range(v, 0, 10),
                        default=1, help="number of polishing iterations "
                        "(default: 1)", metavar="int")
    parser.add_argument("-m", "--min-overlap", dest="min_overlap", metavar="int",
                        type=lambda v: check_int_range(v, 1000, 10000),
                        default=5000, help="minimum overlap between reads "
                        "(default: 5000)")

    parser.add_argument("--resume", action="store_true",
                        dest="resume", default=False,
                        help="resume from the last completed stage")
    parser.add_argument("--resume-from", dest="resume_from", metavar="stage_name",
                        default=None, help="resume from a custom stage")
    #parser.add_argument("--kmer-size", dest="kmer_size",
    #                    type=lambda v: check_int_range(v, 11, 31, require_odd=True),
    #                    default=None, help="kmer size (default: auto)")
    #parser.add_argument("--min-coverage", dest="min_kmer_count",
    #                    type=lambda v: check_int_range(v, 1, 1000),
    #                    default=None, help="minimum kmer coverage "
    #                    "(default: auto)")
    #parser.add_argument("--max-coverage", dest="max_kmer_count",
    #                    type=lambda v: check_int_range(v, 1, 1000),
    #                    default=None, help="maximum kmer coverage "
    #                    "(default: auto)")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("-v", "--version", action="version", version=_version())
    args = parser.parse_args()

    if args.pacbio_raw:
        args.reads = args.pacbio_raw
        args.platform = "pacbio"
        args.read_type = "raw"
    if args.pacbio_corrected:
        args.reads = args.pacbio_corrected
        args.platform = "pacbio"
        args.read_type = "corrected"
    if args.nano_raw:
        args.reads = args.nano_raw
        args.platform = "nano"
        args.read_type = "raw"
    if args.nano_corrected:
        args.reads = args.nano_corrected
        args.platform = "nano"
        args.read_type = "corrected"
    if args.subassemblies:
        args.reads = args.subassemblies
        args.platform = "pacbio"
        args.read_type = "subassemblies"

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    args.out_dir = os.path.abspath(args.out_dir)

    args.log_file = os.path.join(args.out_dir, "flye.log")
    _enable_logging(args.log_file, args.debug,
                    overwrite=not args.resume and not args.resume_from)

    aln.check_binaries()
    pol.check_binaries()
    asm.check_binaries()

    _set_kmer_size(args)
    _set_read_attributes(args)

    try:
        _run(args)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException, repeat.RepeatException,
            ResumeException) as e:
        logger.error(e)
        return 1

    return 0
