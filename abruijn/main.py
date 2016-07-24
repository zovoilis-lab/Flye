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
                       "stage_id" : 0,
                       "error_profile" : ""}

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
    def __init__(self, out_file, log_file):
        super(JobAssembly, self).__init__()
        self.out_file = out_file
        self.log_file = log_file
        self.name = "assembly"

    def run(self):
        reads_order = os.path.join(self.work_dir, "reads_order.fasta")
        asm.assemble(self.args.reads, reads_order, self.args.kmer_size,
                     self.args.min_cov, self.args.max_cov, self.args.coverage,
                     self.args.debug, self.log_file, self.args.threads)
        contigs_fasta = aln.concatenate_contigs(reads_order)
        fp.write_fasta_dict(contigs_fasta, self.out_file)

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
        logger.info("Polishing genome ({0}/{1})".format(self.stage_id,
                                                        self.args.num_iters))
        contigs_fasta = fp.read_fasta_dict(self.in_reference)
        reference_file = os.path.join(self.work_dir, "blasr_ref_{0}.fasta"
                                                        .format(self.stage_id))
        aln.make_blasr_reference(contigs_fasta, reference_file)
        aln.make_alignment(reference_file, self.args.reads, self.args.threads,
                           self.out_alignment)
        os.remove(reference_file)

        Job.run_description["stage_name"] = self.name
        Job.run_description["stage_id"] = self.stage_id


class JobPolishing(Job):
    def __init__(self, in_alignment, in_reference, out_consensus, stage_id):
        super(JobPolishing, self).__init__()
        self.in_alignment = in_alignment
        self.in_reference = in_reference
        self.out_consensus = out_consensus
        self.name = "polishing"
        self.stage_id = stage_id
        self.out_files = [out_consensus]

    def run(self):
        alignment, mean_error = aln.parse_alignment(self.in_alignment)
        if Job.run_description["error_profile"] == "":
            Job.run_description["error_profile"] = \
                                    aln.choose_error_profile(mean_error)

        #bubbles = bbl.get_bubbles(alignment)

        out_patched = os.path.join(self.work_dir, "patched.fasta")
        bbl.patch_genome(alignment, self.in_reference, out_patched)

        polished_fasta = pol.polish(bubbles, self.args.threads,
                                    Job.run_description["error_profile"],
                                    self.work_dir, self.stage_id)
        fp.write_fasta_dict(polished_fasta, self.out_consensus)

        Job.run_description["stage_name"] = self.name
        Job.run_description["stage_id"] = self.stage_id


def _create_job_list(args, work_dir, log_file):
    """
    Build pipeline as a list of consecutive jobs
    """
    jobs = []
    draft_assembly = os.path.join(work_dir, "draft_assembly.fasta")
    jobs.append(JobAssembly(draft_assembly, log_file))

    prev_assembly = draft_assembly
    for i in xrange(args.num_iters):
        alignment_file = os.path.join(work_dir, "blasr_{0}.m5".format(i + 1))
        polished_file = os.path.join(work_dir,
                                     "polished_{0}.fasta".format(i + 1))
        jobs.append(JobAlignment(prev_assembly, alignment_file, i + 1))
        jobs.append(JobPolishing(alignment_file, prev_assembly,
                                 polished_file, i + 1))
        prev_assembly = polished_file

    for i, job in enumerate(jobs):
        job.args = args
        job.work_dir = work_dir

    return jobs


def run(args):
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    work_dir = os.path.abspath(args.out_dir)

    log_file = os.path.join(work_dir, "abruijn.log")
    enable_logging(log_file, args.debug, not args.resume)

    logger.info("Running ABruijn")
    aln.check_binaries()
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


def enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
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
    parser = argparse.ArgumentParser(description="ABruijn: assembly of long and"
                                     " error-prone reads")

    parser.add_argument("reads", metavar="reads",
                        help="path to a file with reads in FASTA format")
    parser.add_argument("out_dir", metavar="out_dir",
                        help="output directory")
    parser.add_argument("coverage", metavar="coverage", type=int,
                        help="estimated assembly coverage")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("--resume", action="store_true",
                        dest="resume", default=False,
                        help="try to resume previous assembly")
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=1, help="number of parallel threads "
                        "(default: 1)")
    parser.add_argument("-i", "--iterations", dest="num_iters", type=int,
                        default=2, help="number of polishing iterations "
                        "(default: 2)")
    parser.add_argument("-k", "--kmer-size", dest="kmer_size", type=int,
                        default=15, help="kmer size (default: 15)")
    parser.add_argument("-m", "--min-cov", dest="min_cov", type=int,
                        default=None, help="minimum kmer coverage "
                        "(default: auto)")
    parser.add_argument("-x", "--max-cov", dest="max_cov", type=int,
                        default=None, help="maximum kmer coverage "
                        "(default: auto)")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    try:
        run(args)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException, ResumeException) as e:
        logger.error("Error: {0}".format(e))
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
