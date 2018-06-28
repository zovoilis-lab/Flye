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
from itertools import combinations

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
import flye.divergence as div
import flye.repeat_resolution as rr

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
        self.out_files["stats"] = os.path.join(work_dir, "assembly_info.txt")
        self.out_files["graph"] = os.path.join(work_dir, "assembly_graph.dot")

    def run(self):
        shutil.copy2(self.contigs_file, self.out_files["contigs"])
        shutil.copy2(self.graph_file, self.out_files["graph"])

        scaffolds = scf.generate_scaffolds(self.contigs_file, self.scaffold_links,
                                           self.out_files["scaffolds"])
        scf.generate_stats(self.repeat_stats, self.polished_stats, scaffolds,
                           self.out_files["stats"])

        logger.info("Final assembly: {0}".format(self.out_files["scaffolds"]))


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
                                             config.vals["min_aln_rate"])
        fp.write_fasta_dict(consensus_fasta, self.out_consensus)


class JobPolishing(Job):
    #CHANGED SO THAT READS ARE NOT IN ARGS
    #Also changed the directory to an input
    def __init__(self, args, work_dir, log_file, reads, in_contigs, polish_dir):
        super(JobPolishing, self).__init__()

        self.args = args
        self.log_file = log_file
        self.in_contigs = in_contigs
        self.reads = reads
        self.polishing_dir = os.path.join(work_dir, polish_dir)

        self.name = "polishing"
        final_contigs = os.path.join(self.polishing_dir,
                                     "polished_{0}.fasta".format(args.num_iters))
        self.out_files["contigs"] = final_contigs
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
            aln.make_alignment(prev_assembly, self.reads, self.args.threads,
                               self.polishing_dir, self.args.platform,
                               alignment_file)

            logger.info("Separating alignment into bubbles")
            
            contigs_info = aln.get_contigs_info(prev_assembly)
            bubbles_file = os.path.join(self.polishing_dir,
                                        "bubbles_{0}.fasta".format(i + 1))
            coverage_stats = \
                bbl.make_bubbles(alignment_file, contigs_info, prev_assembly,
                                 self.args.platform, self.args.threads,
                                 config.vals["min_aln_rate"], bubbles_file)

            logger.info("Correcting bubbles")
            polished_file = os.path.join(self.polishing_dir,
                                         "polished_{0}.fasta".format(i + 1))
            contig_lengths = pol.polish(bubbles_file, self.args.threads,
                                        self.args.platform, self.polishing_dir,
                                        i + 1, polished_file)
            prev_assembly = polished_file

        with open(self.out_files["stats"], "w") as f:
            f.write("seq_name\tlength\tcoverage\n")
            for ctg_id in contig_lengths:
                f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                        contig_lengths[ctg_id], coverage_stats[ctg_id]))


class JobTrestle(Job):
    def __init__(self, args, work_dir, log_file):
        super(JobTrestle, self).__init__()

        self.args = args
        self.resolve_dir = os.path.join(work_dir, "4-trestle")
        self.log_file = log_file
        self.max_iter = 10
        self.buffer_count = 0
        self.min_edge_cov = 10
        self.min_bridge_count = 5
        self.min_bridge_diff = 5

        self.name = "trestle"
        self.out_files["reps"] = os.path.join(self.resolve_dir, 
                                              "resolved_repeats.fasta")
        self.out_files["summary"] = os.path.join(self.resolve_dir, 
                                              "resolve_summary.txt")
        
    def run(self):
        if not os.path.isdir(self.resolve_dir):
            os.mkdir(self.resolve_dir)
            
        logger.info("Running Trestle: resolving unbridged repeats")
        
        #1. Process repeats from graph - each repeat should have its own folder
        logger.info("0.Generating repeat directories and files")
        
        template_name = "template.fasta"
        extended_name = "extended_templates.{0}.{1}.fasta"
        repeat_reads_name = "repeat_reads.fasta"
        init_part_name = "partitioning.0.{0}.txt"
        repeat_label = "repeat_{0}"
        orient_labels = ["forward", "reverse"]
        side_labels = ["in", "out"]
        all_resolved_reps_dict = {}
        
        repeat_list, repeat_edges, all_edge_headers = rr.process_repeats(
                                          self.args.reads, 
                                          self.args.repeats_dump, 
                                          self.args.graph_edges, 
                                          self.resolve_dir, 
                                          template_name,
                                          extended_name,
                                          repeat_reads_name, 
                                          init_part_name, 
                                          repeat_label,
                                          orient_labels,
                                          self.args.min_mult, 
                                          self.args.max_mult,
                                          self.args.extend_len)
        
        logger.info("{0} repeats to be resolved".format(len(repeat_list)))
        for rep_id in sorted(repeat_list):
            logger.info("1.Repeat {0}".format(rep_id))
            repeat_dir = os.path.join(self.resolve_dir, 
                                      repeat_label.format(rep_id))
            orient_reps = [rep_id, -rep_id]
            
            for orientation, rep in zip(orient_labels, orient_reps):
                orient_dir = os.path.join(repeat_dir, orientation)
                template = os.path.join(orient_dir, template_name)
                extended = os.path.join(orient_dir, extended_name)
                repeat_reads = os.path.join(orient_dir, repeat_reads_name)
                
                #2. Polish template and extended templates
                logger.info("a.polishing templates")
                pol_temp_dir = "Polishing.Template"
                pol_temp_job = JobPolishing(self.args, orient_dir, 
                                            self.log_file, [repeat_reads], 
                                            template, pol_temp_dir)
                pol_temp_job.run()
                polished_template = pol_temp_job.out_files["contigs"]
                
                pol_ext_dir = "Polishing.Extended.{0}.{1}"
                polished_extended = {}
                for side in side_labels:
                    for edge_id in repeat_edges[rep][side]:
                        pol_ext_job = JobPolishing(self.args, orient_dir, 
                                           self.log_file, [repeat_reads], 
                                           extended.format(side, edge_id), 
                                           pol_ext_dir.format(side, edge_id))
                        pol_ext_job.run()
                        pol_output = pol_ext_job.out_files["contigs"]
                        polished_extended[(side, edge_id)] = pol_output
                
                #3. Find divergent positions
                logger.info("b.estimating divergence")
                frequency_path = os.path.join(orient_dir, "divergence_frequencies.txt")
                position_path = os.path.join(orient_dir, "tentative_positions.txt")
                summary_path = os.path.join(orient_dir, "divergence_summary.txt")
        
                logger.info("running Minimap2")
                alignment_file = os.path.join(orient_dir, "reads.vs.template.minimap.sam")
                aln.make_alignment(polished_template, [repeat_reads], 
                           self.args.threads, orient_dir, 
                           self.args.platform, alignment_file)
                
                logger.info("finding tentative divergence positions")
                template_info = aln.get_contigs_info(polished_template)
                template_len = template_info[str(rep)].length
                div.find_divergence(alignment_file, polished_template, 
                                    template_info, frequency_path, position_path, 
                                    summary_path, config.vals["min_aln_rate"], 
                                    self.args.platform, self.args.threads,
                                    self.args.sub_thresh, self.args.del_thresh,
                                    self.args.ins_thresh) 
                
                #4. Begin iterations
                partitioning = os.path.join(orient_dir, "partitioning.{0}.{1}.txt")
                           
                cons_align = os.path.join(orient_dir, 
                                          "consensus.{0}.{1}.{2}.vs.template.minimap.sam")
                read_align = os.path.join(orient_dir, 
                                          "reads.vs.consensus.{0}.{1}.{2}.minimap.sam")
                confirmed_pos_path = os.path.join(orient_dir, "confirmed_positions.{0}.{1}.txt")
                edge_reads = os.path.join(orient_dir, "edge_reads.{0}.{1}.{2}.fasta")
                polishing_dir = os.path.join(orient_dir, 
                                             "Polishing.Consensus.{0}.{1}.{2}")
                polished_consensus = {}
                cons_vs_cons = os.path.join(orient_dir, 
                                          "consensus.{0}.{1}.{2}.vs.consensus.{3}.{4}.{5}.minimap.sam")
                side_stats = os.path.join(orient_dir, "{0}_stats.txt")                     
                integrated_stats = os.path.join(orient_dir, "integrated_stats.txt")
                int_confirmed_path = os.path.join(orient_dir, "integrated_confirmed_positions.{0}.{1}.txt")
                resolved_rep_name = "".join(["resolved_repeat_{0}".format(rep),".{0}.fasta"])
                resolved_rep_path = os.path.join(orient_dir, resolved_rep_name)
                """"""
                test_pos = os.path.join(orient_dir, "test_pos.{0}.{1}.txt")
                num_test = 10
                zero_it = 0
                side_it = {s:0 for s in side_labels}
                edge_below_cov = {s:False for s in side_labels}
                dup_part = {s:False for s in side_labels}
                prev_partitionings = {s:set() for s in side_labels}
                #Initialize stats
                for side in side_labels:
                    edge_below_cov[side] = rr.init_side_stats(
                                        rep, side, repeat_edges, self.args, 
                                        self.buffer_count, self.max_iter, 
                                        self.min_edge_cov, position_path,
                                        partitioning.format(zero_it, side), 
                                        prev_partitionings[side], 
                                        template_len, 
                                        side_stats.format(side))
                rr.init_int_stats(rep, repeat_edges, zero_it, position_path, 
                                  partitioning, repeat_reads, template_len, 
                                  integrated_stats)
                logger.info("c.iterative procedure")
                for it in range(1, self.max_iter+1):
                    both_break = True
                    for side in side_labels:
                        if edge_below_cov[side] or dup_part[side]:
                            continue
                        else:
                            logger.info("iteration {0}, '{1}'".format(it, side))
                            both_break = False
                        #4a. Call consensus on partitioned reads
                        for edge_id in sorted(repeat_edges[rep][side]):
                            pol_con_dir = polishing_dir.format(
                                        it, side, edge_id)
                            curr_reads = edge_reads.format(it, side, edge_id)
                            rr.write_edge_reads(
                                        it, side, edge_id,
                                        repeat_reads, 
                                        partitioning.format(it-1, side), 
                                        curr_reads)
                            curr_extended = polished_extended[(side, edge_id)]
                            logger.info("polishing '{0} {1}' reads".format(side, edge_id))
                            pol_con_job = JobPolishing(self.args, 
                                       orient_dir,
                                       self.log_file, 
                                       [curr_reads], 
                                       curr_extended,
                                       pol_con_dir)
                            pol_con_job.run()
                            pol_con_out = pol_con_job.out_files["contigs"]
                            polished_consensus[(it, side, edge_id)] = pol_con_out
                        #4b. Partition reads using divergent positions
                        for edge_id in sorted(repeat_edges[rep][side]):
                            cons_al_file = cons_align.format(it, side, edge_id)
                            aln.make_alignment(polished_template, 
                                               [polished_consensus[(it, side, edge_id)]], 
                                                self.args.threads, 
                                                orient_dir, 
                                                self.args.platform, 
                                                cons_al_file)
                            read_al_file = read_align.format(it, side, edge_id)
                            aln.make_alignment(polished_consensus[(it, side, edge_id)], 
                                               [repeat_reads], 
                                               self.args.threads, 
                                               orient_dir, 
                                               self.args.platform, 
                                               read_al_file)
                        logger.info("partitioning '{0}' reads".format(side))
                        rr.partition_reads(repeat_edges[rep][side], it, side, 
                                           position_path, cons_align, 
                                           polished_template, read_align, 
                                           polished_consensus, 
                                           confirmed_pos_path.format(it, side), 
                                           partitioning.format(it, side), 
                                           all_edge_headers[rep], 
                                           config.vals["min_aln_rate"],
                                           self.buffer_count,
                                           test_pos.format(it, side), num_test)
                        #4c. Write stats file for current iteration
                        edge_pairs = sorted(combinations(repeat_edges[rep][side], 2))
                        for edge_one, edge_two in edge_pairs:
                            cons_cons_file = cons_vs_cons.format(it, side, edge_one, 
                                                                 it, side, edge_two)
                            aln.make_alignment(polished_consensus[(it, side, edge_two)], 
                                               [polished_consensus[(it, side, edge_one)]], 
                                                self.args.threads, 
                                                orient_dir, 
                                                self.args.platform, 
                                                cons_cons_file)
                        edge_below_cov[side], dup_part[side] = rr.update_side_stats(
                                            repeat_edges[rep][side], it, side, 
                                            cons_align, polished_template, 
                                            config.vals["min_aln_rate"], 
                                            confirmed_pos_path.format(it, side), 
                                            partitioning.format(it, side), 
                                            self.min_edge_cov, 
                                            prev_partitionings[side], 
                                            side_stats.format(side))
                        side_it[side] = it
                    rr.update_int_stats(rep, repeat_edges, side_it, cons_align, 
                                        polished_template, 
                                        template_len,
                                        config.vals["min_aln_rate"], 
                                        confirmed_pos_path, int_confirmed_path, 
                                        partitioning, integrated_stats)
                    if both_break:
                        break
                logger.info("writing stats files")
                for side in side_labels:
                    rr.finalize_side_stats(repeat_edges[rep][side], side_it[side], 
                                           side, cons_align, polished_template, 
                                        config.vals["min_aln_rate"], 
                                        cons_vs_cons, polished_consensus, 
                                        confirmed_pos_path.format(side_it[side], side), 
                                        partitioning.format(side_it[side], side), 
                                        self.max_iter, edge_below_cov[side],
                                        dup_part[side], 
                                        side_stats.format(side))
                bridged, repeat_seqs = rr.finalize_int_stats(rep, repeat_edges, side_it, 
                                                cons_align, polished_template, 
                                                template_len, 
                                                config.vals["min_aln_rate"], 
                                                cons_vs_cons, 
                                                polished_consensus, 
                                                int_confirmed_path, 
                                                partitioning, 
                                                self.min_bridge_count, 
                                                self.min_bridge_diff, 
                                                integrated_stats, 
                                                resolved_rep_path)
                all_resolved_reps_dict.update(repeat_seqs)
        fp.write_fasta_dict(all_resolved_reps_dict, self.out_files["reps"])

def _create_job_list(args, work_dir, log_file):
    """
    Build pipeline as a list of consecutive jobs
    """
    jobs = []

    #Resolve Unbridged Repeats
    jobs.append(JobTrestle(args, work_dir, log_file))
        
    """
        
        
        #3. Iteration 0 starts
        #Define current partitioning as the set of starting reads
        #Define the initial template as the precursor sequences of the repeat copies
        #Align the reads to the template, polish them, and identify starting points
        #Collect output for this iteration and save it
        
        #4. Start Iterations, 1-max_iter
    

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
                            scaffold_links))"""

    return jobs


def _set_kmer_size(args):
    """
    Select k-mer size based on the target genome size
    """
    if args.genome_size.isdigit():
        args.genome_size = int(args.genome_size)
    else:
        args.genome_size = human2bytes(args.genome_size.upper())

    logger.debug("Genome size: {0}".format(args.genome_size))
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

    if args.max_mult < args.min_mult:
        raise ResumeException("Max-mult cannot be greater than min-mult")

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

    """Repeat Resolutions inputs
    -repeats_dump
    -graph_final.fasta
    -all_reads for whole genome
    """
    parser.add_argument("-r", "--read_files", dest="read_files",
                        metavar="read_files", required=True,
                        help="reads file for the entire assembly")
    parser.add_argument("-d", "--repeats-dump", dest="repeats_dump",
                        metavar="repeats", required=True,
                        help="repeats_dump file from Flye assembly")
    parser.add_argument("-g", "--graph-edges", dest="graph_edges",
                        metavar="graph", required=True,
                        help="graph_final.fasta file from Flye assembly")
    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("-u", "--min-mult", dest="min_mult",
                        type=lambda v: check_int_range(v, 2, 100),
                        default=2, metavar="int", help="minimum multiplicity "
                        "to attempt to resolve")
    parser.add_argument("-x", "--max-mult", dest="max_mult",
                        type=lambda v: check_int_range(v, 2, 100),
                        default=2, metavar="int", help="maximum multiplicity "
                        "to attempt to resolve")
    parser.add_argument("-e", "--extend-len", dest="extend_len",
                        type=lambda v: check_int_range(v, 10, 100000),
                        default=10000, metavar="int", help="length to extend "
                        "into repeat-adjacent edges")

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
    parser.add_argument("-s", "--sub_thresh", dest="sub_thresh", metavar="float",
                        default=0.1, help="threshold for substitution calls "
                        "(default: 0.1)")
    parser.add_argument("-l", "--del_thresh", dest="del_thresh", metavar="float",
                        default=0.2, help="threshold for deletion calls "
                        "(default: 0.2)")
    parser.add_argument("-n", "--ins_thresh", dest="ins_thresh", metavar="float",
                        default=0.3, help="threshold for insertion calls "
                        "(default: 0.3)")

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

    args.reads = [args.read_files]
    args.platform = "pacbio"
    args.read_type = "raw"
    
    """if args.pacbio_raw:
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
    """
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    args.out_dir = os.path.abspath(args.out_dir)

    args.log_file = os.path.join(args.out_dir, "flye.log")
    _enable_logging(args.log_file, args.debug,
                    overwrite=not args.resume and not args.resume_from)

    aln.check_binaries()
    pol.check_binaries()
    asm.check_binaries()

    #_set_kmer_size(args)
    #_set_read_attributes(args)

    try:
        _run(args)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException, repeat.RepeatException,
            rr.ProcessingException, ResumeException) as e:
        logger.error(e)
        return 1

    return 0
