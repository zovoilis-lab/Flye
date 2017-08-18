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

import abruijn.alignment as aln
import abruijn.bubbles as bbl
import abruijn.polish as pol
import abruijn.fasta_parser as fp
import abruijn.assemble as asm
import abruijn.repeat_graph as repeat
import abruijn.consensus as cons
from abruijn.__version__ import __version__


logger = logging.getLogger()

class ResumeException(Exception):
    pass

class Job(object):
    """
    Describes an abstract list of jobs with persistent
    status that can be resumed
    """
    run_description = {"stage_name" : "",
                       "stage_id" : 0}

    def __init__(self):
        self.name = None
        self.args = None
        self.stage_id = 0
        self.work_dir = None
        self.out_files = []

    def run(self):
        pass

    def save(self, save_file):
        with open(save_file, "w") as fp:
            json.dump(Job.run_description, fp)

    def load(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)
            Job.run_description = data

    def can_resume(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)

            if (data["stage_name"] != self.name or
                data["stage_id"] != self.stage_id):
                return False

            for file in self.out_files:
                if not os.path.exists(file):
                    return False

            return True


class JobAssembly(Job):
    def __init__(self, out_assembly, log_file):
        super(JobAssembly, self).__init__()
        self.out_assembly = out_assembly
        self.log_file = log_file
        self.name = "assembly"
        self.out_files = [out_assembly]

    def run(self):
        #reads_order = os.path.join(self.work_dir, "reads_order.fasta")
        asm.assemble(self.args, self.out_assembly, self.log_file)
        #contigs_fasta = aln.concatenate_contigs(reads_order)
        #fp.write_fasta_dict(contigs_fasta, self.out_assembly)

        Job.run_description["stage_name"] = self.name


class JobRepeat(Job):
    def __init__(self, in_assembly, out_folder, log_file):
        super(JobRepeat, self).__init__()
        self.in_assembly = in_assembly
        self.log_file = log_file
        self.out_folder = out_folder
        self.name = "repeat"

        edges_sequences = os.path.join(out_folder, "graph_final.fasta")
        repeat_graph = os.path.join(out_folder, "graph_final.gfa")
        self.out_files = [edges_sequences, repeat_graph]

    def run(self):
        logger.info("Performing repeat analysis")
        repeat.analyse_repeats(self.args, self.in_assembly, self.out_folder,
                               self.log_file)
        Job.run_description["stage_name"] = self.name


class JobAlignment(Job):
    def __init__(self, in_reference, out_alignment, stage_id):
        super(JobAlignment, self).__init__()
        self.in_reference = in_reference
        self.out_alignment = out_alignment
        self.name = "alignment"
        self.stage_id = stage_id
        self.out_files = [out_alignment]

    def run(self):
        #logger.info("Polishing genome ({0}/{1})".format(self.stage_id,
        #                                                self.args.num_iters))
        logger.info("Running BLASR")
        contigs_fasta = fp.read_fasta_dict(self.in_reference)
        reference_file = os.path.join(self.work_dir, "blasr_ref_{0}.fasta"
                                                        .format(self.stage_id))
        aln.make_blasr_reference(contigs_fasta, reference_file)
        aln.make_alignment(reference_file, self.args.reads, self.args.threads,
                           self.out_alignment)
        os.remove(reference_file)

        Job.run_description["stage_name"] = self.name
        Job.run_description["stage_id"] = self.stage_id


class JobConsensus(Job):
    def __init__(self, in_contigs, in_alignment, out_consensus):
        super(JobConsensus, self).__init__()

        self.in_contigs = in_contigs
        self.in_alignment = in_alignment
        self.out_consensus = out_consensus
        self.name = "consensus"
        self.out_files = [out_consensus]

    def run(self):
        logger.info("Computing rough consensus")
        contigs_info = aln.get_contigs_info(self.in_contigs)
        consensus_fasta = cons.get_consensus(self.in_alignment, contigs_info,
                                             self.args.threads)
        fp.write_fasta_dict(consensus_fasta, self.out_consensus)

        Job.run_description["stage_name"] = self.name


class JobPolishing(Job):
    def __init__(self, in_contigs, in_alignment, out_consensus,
                 stage_id, seq_platform):
        super(JobPolishing, self).__init__()

        self.in_contigs = in_contigs
        self.in_alignment = in_alignment
        self.out_consensus = out_consensus
        self.name = "polishing"
        self.stage_id = stage_id
        self.seq_platform = seq_platform
        self.out_files = [out_consensus]

    def run(self):
        logger.info("Polishing genome ({0}/{1})".format(self.stage_id,
                                                        self.args.num_iters))
        contigs_info = aln.get_contigs_info(self.in_contigs)

        logger.info("Separating alignment into bubbles")
        bubbles = bbl.get_bubbles(self.in_alignment, contigs_info,
                                  self.seq_platform, self.args.threads)
        logger.info("Correcting bubbles")
        polished_fasta = pol.polish(bubbles, self.args.threads,
                                    self.seq_platform, self.work_dir,
                                    self.stage_id)
        fp.write_fasta_dict(polished_fasta, self.out_consensus)

        Job.run_description["stage_name"] = self.name
        Job.run_description["stage_id"] = self.stage_id


def _create_job_list(args, work_dir, log_file):
    """
    Build pipeline as a list of consecutive jobs
    """
    jobs = []

    #Assembly job
    draft_assembly = os.path.join(work_dir, "draft_assembly.fasta")
    jobs.append(JobAssembly(draft_assembly, log_file))

    #Pre-polishing
    #alignment_file = os.path.join(work_dir, "blasr_0.m5")
    #pre_polished_file = os.path.join(work_dir, "polished_0.fasta")
    #jobs.append(JobAlignment(draft_assembly, alignment_file, 0))
    #jobs.append(JobConsensus(draft_assembly, alignment_file,
    #                         pre_polished_file))

    #Repeat analysis
    edges_sequences = os.path.join(work_dir, "graph_final.fasta")
    jobs.append(JobRepeat(draft_assembly, work_dir, log_file))

    #Full polishing
    #prev_assembly = edges_sequences
    #for i in xrange(args.num_iters):
    #    alignment_file = os.path.join(work_dir, "blasr_{0}.m5".format(i + 1))
    #    polished_file = os.path.join(work_dir,
    #                                 "polished_{0}.fasta".format(i + 1))
    #    jobs.append(JobAlignment(prev_assembly, alignment_file, i + 1))
    #    jobs.append(JobPolishing(prev_assembly, alignment_file, polished_file,
    #                             i + 1, args.sequencing_platform))
    #    prev_assembly = polished_file

    for job in jobs:
        job.args = args
        job.work_dir = work_dir

    return jobs


def _run(args):
    """
    Runs the pipeline
    """
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    work_dir = os.path.abspath(args.out_dir)

    log_file = os.path.join(work_dir, "abruijn.log")
    _enable_logging(log_file, args.debug, not args.resume)

    logger.info("Running ABruijn")
    #aln.check_binaries()
    pol.check_binaries()
    asm.check_binaries()

    save_file = os.path.join(work_dir, "abruijn.save")
    jobs = _create_job_list(args, work_dir, log_file)

    current_job = 0
    can_resume = False
    if args.resume:
        for i in xrange(len(jobs)):
            if jobs[i].can_resume(save_file):
                jobs[i].load(save_file)
                can_resume = True
                current_job = i + 1

    if args.resume:
        if can_resume:
            logger.info("Resuming previous run")
        else:
            raise ResumeException("Can't resume previous run")

    for i in xrange(current_job, len(jobs)):
        jobs[i].run()
        jobs[i].save(save_file)

    logger.info("Done! Your assembly is in file: " + jobs[-1].out_files[0])


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


def main():
    def check_int_range(value, min_val, max_val, require_odd=False):
        ival = int(value)
        if ival < min_val or ival > max_val:
             raise argparse.ArgumentTypeError("value should be in "
                            "range [{0}, {1}]".format(min_val, max_val))
        if require_odd and ival % 2 == 0:
            raise argparse.ArgumentTypeError("should be an odd number")
        return ival


    parser = argparse.ArgumentParser(description="ABruijn: assembly of long and"
                                     " error-prone reads")

    parser.add_argument("reads", metavar="reads",
                        help="path to reads file (FASTA format)")
    parser.add_argument("out_dir", metavar="out_dir",
                        help="output directory")
    parser.add_argument("coverage", metavar="coverage",
                        type=lambda v: check_int_range(v, 1, 1000),
                        help="estimated assembly coverage")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("--resume", action="store_true",
                        dest="resume", default=False,
                        help="try to resume previous assembly")
    parser.add_argument("-t", "--threads", dest="threads",
                        type=lambda v: check_int_range(v, 0, 9999),
                        default=1, help="number of parallel threads "
                        "(default: 1)")
    parser.add_argument("-i", "--iterations", dest="num_iters",
                        type=lambda v: check_int_range(v, 0, 10),
                        default=1, help="number of polishing iterations "
                        "(default: 1)")
    parser.add_argument("-p", "--platform", dest="sequencing_platform",
                        default="pacbio",
                        choices=["pacbio", "nano", "pacbio_hi_err"],
                        help="sequencing platform (default: pacbio)")
    parser.add_argument("-k", "--kmer-size", dest="kmer_size",
                        type=lambda v: check_int_range(v, 11, 31, require_odd=True),
                        default=15, help="kmer size (default: 15)")
    parser.add_argument("-o", "--min-overlap", dest="min_overlap",
                        type=lambda v: check_int_range(v, 2000, 10000),
                        default=5000, help="minimum overlap between reads "
                        "(default: 5000)")
    parser.add_argument("-m", "--min-coverage", dest="min_kmer_count",
                        type=lambda v: check_int_range(v, 1, 1000),
                        default=None, help="minimum kmer coverage "
                        "(default: auto)")
    parser.add_argument("-x", "--max-coverage", dest="max_kmer_count",
                        type=lambda v: check_int_range(v, 1, 1000),
                        default=None, help="maximum kmer coverage "
                        "(default: auto)")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    try:
        _run(args)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException, ResumeException) as e:
        logger.error("Error: {0}".format(e))
        return 1

    return 0
