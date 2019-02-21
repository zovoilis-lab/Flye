#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Created on Wed Jan  4 03:50:31 2017

@author: jeffrey_yuan
"""

import os
import logging
from itertools import combinations
import copy
import multiprocessing, signal
from sklearn.cluster import AgglomerativeClustering

import flye.polishing.alignment as flye_aln
from flye.polishing.alignment import Alignment
import flye.utils.fasta_parser as fp
import flye.config.py_cfg as config
import flye.polishing.polish as pol

import flye.trestle.divergence as div
import flye.trestle.trestle_config as trestle_config

logger = logging.getLogger()


def test_divergence(args, phasing_dir, uniques_dict):
    SUB_THRESH = trestle_config.vals["sub_thresh"]
    DEL_THRESH = trestle_config.vals["del_thresh"]
    INS_THRESH = trestle_config.vals["ins_thresh"]
    MIN_ALN_RATE = trestle_config.vals["min_aln_rate"]
    NUM_POL_ITERS = trestle_config.vals["num_pol_iters"]
    ORIENT_CONFIG = trestle_config.vals["orientations_to_run"]
    PHASE_MIN_LONGEST_DIV = trestle_config.vals["phase_min_longest_div"]
    all_file_names = define_file_names()
    (pol_dir_names, initial_file_names,
     pre_file_names, div_file_names, aln_names,
     middle_file_names, output_file_names) = all_file_names
    (edge_label, template_name, initial_reads_name, 
     left_reads_name, right_reads_name) = initial_file_names
    
    reads_dict = {}
    for read_file in args.reads:
        reads_dict.update(fp.read_sequence_dict(read_file))
    if not reads_dict:
        raise ProcessingException("No reads found from {0}".format(args.reads))
    
    longest_edge = None
    longest_edge_len = 0
    for edge in sorted(uniques_dict, reverse=True):
        #Checks presence of reverse strand
        #One run processes both forward and reverse strand of edge
                
        if edge <= 0:
            continue
        
        valid_edge = True
        if -edge not in uniques_dict:
            logger.debug("Edge {0} missing reverse strand".format(edge))
            valid_edge = False
        if not valid_edge:
            continue
        
        template_seq = uniques_dict[edge].sequences["template"]
        if len(template_seq) > longest_edge_len:
            longest_edge_len = len(template_seq)
            longest_edge = edge
    
    haploid = True
    if not longest_edge:
        return True
    
    logger.debug("Testing divergence on longest edge {0}".format(longest_edge))
    
    edge = longest_edge
    #Make edge dir
    edge_dir = os.path.join(phasing_dir, edge_label.format(edge))
    if not os.path.isdir(edge_dir):
        os.mkdir(edge_dir)
    
    #Only run the forward strand for the test if ORIENT_CONFIG = "both"
    curr_label, curr_edge = None, None
    if ORIENT_CONFIG == "forward":
        curr_label, curr_edge = ("forward", edge)
    elif ORIENT_CONFIG == "reverse":
        curr_label, curr_edge = ("reverse", -edge)
    elif ORIENT_CONFIG == "both":
        curr_label, curr_edge = ("forward", edge)
    orient_path = os.path.join(edge_dir, curr_label)
    if not os.path.isdir(orient_path):
        os.mkdir(orient_path)
    template_path = os.path.join(orient_path, template_name)
    initial_reads_path = os.path.join(orient_path, initial_reads_name)
    all_reads_list = uniques_dict[curr_edge].all_reads
    
    template_dict = {}
    initial_reads_dict = {}
    
    template_seq = uniques_dict[curr_edge].sequences["template"]
    #if curr_label == "reverse":
    #    template_seq = fp.reverse_complement(graph_dict[edge])
    template_dict[curr_edge] = template_seq
    
    #rev_comp of read will be written if the header is -h
    
    for header in all_reads_list:
        if (not header) or (header[0] != '+' and header[0] != '-'):
            raise ProcessingException(
                "All reads format not recognized: {0}".format(header))
        if header[1:] not in reads_dict:
            raise ProcessingException(
                "Read header {0} not in any of {1}".format(
                    header[1:], args.reads))
        
        seq = reads_dict[header[1:]]
        if header[0] == '-':
            seq = fp.reverse_complement(seq)
        initial_reads_dict[header[1:]] = seq
    
    if template_dict and template_dict.values()[0]:
        fp.write_fasta_dict(template_dict, template_path)
    if initial_reads_dict and initial_reads_dict.values()[0]:
        fp.write_fasta_dict(initial_reads_dict, initial_reads_path)
    
    if not template_dict:
        raise ProcessingException("No template {0} found".format(
                                        curr_edge))
    if not initial_reads_dict:
        raise ProcessingException("No repeat reads {0} found".format(
                                        curr_edge))
    

    pol_temp_name, pol_cons_name = pol_dir_names
    div_freq_name, div_pos_name, div_summ_name = div_file_names
    (reads_template_aln_name, cons_temp_aln_name,
     cut_cons_temp_aln_name, reads_cons_aln_name) = aln_names
        
    #Polish template
    logger.debug("Polishing templates")
    pol_temp_dir = os.path.join(orient_path, pol_temp_name)
    if not os.path.isdir(pol_temp_dir):
        os.mkdir(pol_temp_dir)
    polished_template, _stats = \
        pol.polish(template_path, [initial_reads_path], pol_temp_dir, NUM_POL_ITERS,
                   args.threads, args.platform, output_progress=False)

    if not os.path.getsize(polished_template):
        return True
    
    #3. Find divergent positions
    logger.debug("Estimating divergence")
    frequency_path = os.path.join(orient_path, div_freq_name)
    position_path = os.path.join(orient_path, div_pos_name)
    summary_path = os.path.join(orient_path, div_summ_name)
    
    #logger.info("running Minimap2")
    alignment_file = os.path.join(orient_path, reads_template_aln_name)
    template_len = 0.0
    if os.path.getsize(polished_template):
        flye_aln.make_alignment(polished_template, [initial_reads_path], 
                   args.threads, orient_path, args.platform, 
                   alignment_file, reference_mode=True, sam_output=True)
        template_info = flye_aln.get_contigs_info(polished_template)
        template_len = template_info[str(edge)].length

    logger.debug("Finding tentative divergence positions")
    div.find_divergence(alignment_file, polished_template, 
                        template_info, frequency_path, position_path, 
                        summary_path, MIN_ALN_RATE, 
                        args.platform, args.threads,
                        SUB_THRESH, DEL_THRESH, INS_THRESH)

    pos_headers, pos = div.read_positions(position_path)
    num_pos = len(pos["total"])
    div_rate = num_pos / float(template_len)
    logger.debug("Edge {0}: divergence rate = {1} pos / {2} bp = {3}".format(
                        curr_edge, num_pos, template_len, div_rate))
    if div_rate >= PHASE_MIN_LONGEST_DIV:
        haploid = False
    return haploid

def phase_uniques(args, phasing_dir, uniques_info, summ_file,
                    phased_seqs):
    print "Here 1"
    all_file_names = define_file_names()
    (pol_dir_names, initial_file_names,
     pre_file_names, div_file_names, aln_names,
     middle_file_names, output_file_names) = all_file_names
    all_phased_dict = {}
    all_summaries = []
    init_summary(summ_file)
    #1. Process repeats from graph - generates a folder for each repeat
    logger.debug("Finding unique edges")
    process_outputs = process_repeats(args.reads, uniques_info,
                                      phasing_dir,
                                      initial_file_names)
    edge_list, all_edge_headers = process_outputs
    logger.info("Unique edges: {0}".format(len(edge_list)))
    #if not repeat_list:
    #    return
    print "Here 2"
    #Resolve every repeat in a separate thread
    def _thread_worker(func_args, log_file, results_queue, error_queue):
        try:
            #each thread logs to a separate file
            log_formatter = \
                logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                  "%(message)s", "%Y-%m-%d %H:%M:%S")
            file_handler = logging.FileHandler(log_file, mode="a")
            file_handler.setFormatter(log_formatter)
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
            logger.addHandler(file_handler)

            result = phase_each_edge(*func_args)
            results_queue.put(result)

        except Exception as e:
            error_queue.put(e)
    print "Here 3"
    job_chunks = [edge_list[i:i + args.threads]
              for i in xrange(0, len(edge_list), args.threads)]

    for job_chunk in job_chunks:
        manager = multiprocessing.Manager()
        results_queue = manager.Queue()
        error_queue = manager.Queue()
        print "Here 4"
        
        phase_threads = max(1, args.threads / len(job_chunk))
        orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
        threads = []
        for edge_id in sorted(job_chunk):
            func_args = (edge_id, all_edge_headers, args, phasing_dir,
                         uniques_info, all_file_names, phase_threads)
            log_file = os.path.join(phasing_dir,
                                    "edge_{0}".format(edge_id), "log.txt")
            threads.append(multiprocessing.Process(target=_thread_worker,
                                        args=(func_args, log_file,
                                              results_queue, error_queue)))
        signal.signal(signal.SIGINT, orig_sigint)
        print "Here 5"

        for t in threads:
            t.start()
        try:
            for t in threads:
                t.join()
        except KeyboardInterrupt:
            for t in threads:
                t.terminate()
            if t.exitcode == -9:
                logger.error("Looks like the system ran out of memory")
            if t.exitcode != 0:
                raise Exception("One of the processes exited with code: {0}"
                                .format(t.exitcode))
        if not error_queue.empty():
            raise error_queue.get()
        print "Here 6"
        while not results_queue.empty():
            phased_dict, summary_list = results_queue.get()
            all_phased_dict.update(phased_dict)
            all_summaries.extend(summary_list)
    print "Here 7"
    fp.write_fasta_dict(all_phased_dict, phased_seqs)
    num_phased = 0
    for summ_items in sorted(all_summaries):
        if len(summ_items) > 5 and summ_items[5]:
            num_phased += 1
        update_summary(summ_items, summ_file)
    logger.info("Phased: {0}".format(num_phased))


def phase_each_edge(edge_id, all_edge_headers, args,
                    phasing_dir, uniques_info, all_file_names,
                    num_threads):
    SUB_THRESH = trestle_config.vals["sub_thresh"]
    DEL_THRESH = trestle_config.vals["del_thresh"]
    INS_THRESH = trestle_config.vals["ins_thresh"]
    MAX_ITER = trestle_config.vals["max_iter"]
    MIN_ALN_RATE = trestle_config.vals["min_aln_rate"]
    NUM_POL_ITERS = trestle_config.vals["num_pol_iters"]
    ORIENT_CONFIG = trestle_config.vals["orientations_to_run"]
    PHASE_MIN_ANY_DIV = trestle_config.vals["phase_min_any_div"]
    zero_it = 0
    side_labels = ["left", "right"]
    phase_labels = [1, 2]
    print "Here 8"
    
    (pol_dir_names, initial_file_names,
     pre_file_names, div_file_names, aln_names,
     middle_file_names, output_file_names) = all_file_names

    pol_temp_name, pol_cons_name = pol_dir_names
    (edge_label, template_name, initial_reads_name, 
     left_reads_name, right_reads_name) = initial_file_names
    (pre_phased_reads_name, pre_read_aln_name, pre_partitioning_name,
     partitioning_name) = pre_file_names
    div_freq_name, div_pos_name, div_summ_name = div_file_names
    (reads_template_aln_name, cons_temp_aln_name,
     cut_cons_temp_aln_name, reads_cons_aln_name) = aln_names
    (confirmed_pos_name, phased_reads_name,
     cut_cons_name, cons_vs_cons_name) = middle_file_names
    (side_stats_name, int_stats_name, int_confirmed_pos_name,
     phased_edge_name, res_vs_res_name) = output_file_names
    print "Here 9"
    logger.info("Phasing unique edge {0}: {1}" \
        .format(edge_id, uniques_info[edge_id].path))
    edge_dir = os.path.join(phasing_dir,
                              edge_label.format(edge_id))

    run_orientations = []
    if ORIENT_CONFIG == "forward":
        run_orientations = [("forward", edge_id)]
    elif ORIENT_CONFIG == "reverse":
        run_orientations = [("reverse", -edge_id)]
    elif ORIENT_CONFIG == "both":
        run_orientations = [("forward", edge_id), ("reverse", -edge_id)]
    edge_phased = False
    phased_dict = {}
    summary_list = []
    print "Here 10"
    for orientation, edge in run_orientations:
        logger.debug("Orientation: " + orientation)
        orient_dir = os.path.join(edge_dir, orientation)
        if not os.path.isdir(orient_dir):
            os.mkdir(orient_dir)
        template = os.path.join(orient_dir, template_name)
        initial_reads = os.path.join(orient_dir, initial_reads_name)
        term_bool = {s:False for s in side_labels}
        
        #2. Polish template
        logger.debug("Polishing templates")
        pol_temp_dir = os.path.join(orient_dir, pol_temp_name)
        if not os.path.isdir(pol_temp_dir):
            os.mkdir(pol_temp_dir)
        polished_template, _stats = \
            pol.polish(template, [initial_reads], pol_temp_dir, NUM_POL_ITERS,
                       num_threads, args.platform, output_progress=False)
        print "Here 11"
        if not os.path.getsize(polished_template):
            for side in side_labels:
                term_bool[side] = True
        
        #3. Find divergent positions
        logger.debug("Estimating divergence")
        frequency_path = os.path.join(orient_dir, div_freq_name)
        position_path = os.path.join(orient_dir, div_pos_name)
        summary_path = os.path.join(orient_dir, div_summ_name)
        
        #logger.info("running Minimap2")
        alignment_file = os.path.join(orient_dir, reads_template_aln_name)
        template_len = 0.0
        if os.path.getsize(polished_template):
            flye_aln.make_alignment(polished_template, [initial_reads], 
                       num_threads, orient_dir, args.platform, 
                       alignment_file, reference_mode=True, sam_output=True)
            template_info = flye_aln.get_contigs_info(polished_template)
            template_len = template_info[str(edge)].length
        print "Here 12"
        logger.debug("Finding tentative divergence positions")
        div.find_divergence(alignment_file, polished_template, 
                            template_info, frequency_path, position_path, 
                            summary_path, MIN_ALN_RATE, 
                            args.platform, num_threads,
                            SUB_THRESH, DEL_THRESH, INS_THRESH) 
        pos_headers, pos = div.read_positions(position_path)
        div_rate = len(pos["total"]) / float(template_len)
        if div_rate < PHASE_MIN_ANY_DIV:
            logger.debug("Edge {0} div rate {1} below minimum".format(edge, 
                                                                     div_rate))
            for side in side_labels:
                term_bool[side] = True
        read_endpoints = find_read_endpoints(alignment_file, 
                                             polished_template)
        avg_cov = find_coverage(frequency_path)
        print "Here 13"
        #4. Initialize paths, variables, and stats
        #pre_partitioning = os.path.join(orient_dir, pre_partitioning_name)
        #pre_phased_reads = os.path.join(orient_dir, pre_phased_reads_name)
        #pre_read_align = os.path.join(orient_dir, pre_read_aln_name)
        partitioning = os.path.join(orient_dir, partitioning_name)
        cons_align = os.path.join(orient_dir, cons_temp_aln_name)
        cut_cons_align = os.path.join(orient_dir, cut_cons_temp_aln_name)
        read_align = os.path.join(orient_dir, reads_cons_aln_name)
        confirmed_pos_path = os.path.join(orient_dir, confirmed_pos_name)
        phased_reads = os.path.join(orient_dir, phased_reads_name)
        cut_cons = os.path.join(orient_dir, cut_cons_name)
        polishing_dir = os.path.join(orient_dir, pol_cons_name)
        cons_vs_cons = os.path.join(orient_dir, cons_vs_cons_name)
        side_stats = os.path.join(orient_dir, side_stats_name)
        integrated_stats = os.path.join(orient_dir, int_stats_name)
        int_confirmed_path = os.path.join(orient_dir, 
                                          int_confirmed_pos_name)
        phased_edge_path = os.path.join(orient_dir, phased_edge_name)
        res_vs_res = os.path.join(orient_dir, res_vs_res_name)
        
        '''#5. Re-align reads to extended and initialize partitioning 0
        logger.debug("Checking initial set of phased reads")
        for side in side_labels:
            for edge_id in repeat_edges[rep][side]:
                write_edge_reads(zero_it, side, edge_id,
                                 repeat_reads, 
                                 pre_partitioning.format(side), 
                                 pre_edge_reads.format(side, edge_id))
                flye_aln.make_alignment(polished_extended[(side, edge_id)], 
                                        [pre_edge_reads.format(side, edge_id)], 
                                        num_threads, orient_dir, args.platform, 
                                        pre_read_align.format(side, edge_id),
                                        reference_mode=True, sam_output=True)
            init_partitioning(repeat_edges[rep][side], 
                              side, pre_partitioning.format(side), 
                              pre_read_align, polished_extended, 
                              partitioning.format(zero_it, side))'''
        
        #5. Get initial set of read partitions
        #given the alignment of all reads against the genome,
        #the divergent positions, etc.
        print "Edge",edge
        window_size = 1000
        left_reads_file = os.path.join(orient_dir, left_reads_name)
        right_reads_file = os.path.join(orient_dir, right_reads_name)
        left_part = partitioning.format(zero_it, "left")
        right_part = partitioning.format(zero_it, "right")
        win_inds = _find_div_region(position_path, polished_template, 
                                    alignment_file, left_part, right_part, 
                                    window_size, phase_labels, 
                                    left_reads_file, right_reads_file, 
                                    initial_reads, all_edge_headers[edge])
        win_start, win_end = win_inds
        print "win_st, win_end", win_start, win_end
        cut_consensus = {}
        side_it = {s:0 for s in side_labels}
        iter_pairs = []
        phase_below_cov = {s:False for s in side_labels}
        dup_part = {s:False for s in side_labels}
        prev_partitionings = {s:set() for s in side_labels}
        print "Here 14"
        #6. Initialize stats
        for side in side_labels:
            phase_below_cov[side] = init_side_stats(
                                edge, side, phase_labels, args.min_overlap, 
                                position_path, 
                                partitioning.format(zero_it, side), 
                                prev_partitionings[side], 
                                template_len, 
                                side_stats.format(side))
        init_int_stats(edge, side_labels, phase_labels, zero_it, position_path, 
                       partitioning, initial_reads, template_len, 
                       avg_cov, integrated_stats)
        print "Here 15"
        #7. Start iterations
        logger.debug("Iterative procedure")
        for it in range(1, MAX_ITER + 1):
            both_break = True
            for side in side_labels:
                if (phase_below_cov[side] or dup_part[side] or 
                    term_bool[side]):
                    continue
                else:
                    logger.debug("Iteration {0}, '{1}'".format(it, side))
                    both_break = False
                print "Here 16"
                if side == "left":
                    side_reads = left_reads_file
                elif side == "right":
                    side_reads = right_reads_file
                for phase_id in phase_labels:
                    #7a. Call consensus on partitioned reads
                    pol_con_dir = polishing_dir.format(
                                it, side, phase_id)
                    curr_reads = phased_reads.format(it, side, phase_id)
                    write_phased_reads(
                                it, side, phase_id,
                                side_reads, 
                                partitioning.format(it - 1, side), 
                                curr_reads)
                    pol_reads_str = "\tPolishing '{0} {1}' reads"
                    logger.debug(pol_reads_str.format(side, phase_id))
                    print "Here 17"
                    
                    if not os.path.isdir(pol_con_dir):
                        os.mkdir(pol_con_dir)
                    pol_con_out, _stats = \
                        pol.polish(polished_template, [curr_reads], pol_con_dir, 
                                   NUM_POL_ITERS, num_threads, args.platform,
                                   output_progress=False)
                    print "Here 18"
                    #7b. Cut consensus where coverage drops
                    if side == "left":
                        win_cutpoint = win_end
                    elif side == "right":
                        win_cutpoint = win_start
                    cutpoint = locate_consensus_cutpoint(
                                    side, read_endpoints,
                                    curr_reads)
                    if os.path.getsize(pol_con_out):
                        cons_al_file = cons_align.format(it, side, phase_id)
                        flye_aln.make_alignment(polished_template, [pol_con_out], 
                                                num_threads, orient_dir, 
                                                args.platform, cons_al_file,
                                                reference_mode=True, sam_output=True)
                    else:
                        term_bool[side] = True
                    curr_cut_cons = cut_cons.format(it, side, phase_id)
                    cut_consensus[(it, side, phase_id)] = curr_cut_cons
                    if os.path.getsize(pol_con_out) and os.path.isfile(cons_al_file):
                        truncate_consensus(side, cutpoint, cons_al_file, 
                                           polished_template,
                                           pol_con_out, curr_cut_cons, 
                                           win_cutpoint)
                    else:
                        term_bool[side] = True
                    print "Here 19"
                    #7c. Align consensuses to template 
                    #    and reads to consensuses
                    if os.path.isfile(curr_cut_cons):
                        cut_cons_al_file = cut_cons_align.format(it, side, 
                                                                 phase_id)
                        flye_aln.make_alignment(polished_template, [curr_cut_cons], 
                                                num_threads, orient_dir, 
                                                args.platform, cut_cons_al_file,
                                                reference_mode=True, sam_output=True)
                        read_al_file = read_align.format(it, side, phase_id)
                        flye_aln.make_alignment(curr_cut_cons, [side_reads], 
                                                num_threads, orient_dir, 
                                                args.platform, read_al_file,
                                                reference_mode=True, sam_output=True)
                    else:
                        term_bool[side] = True
                print "Here 20"
                #7d. Partition reads using divergent positions
                logger.debug("\tPartitioning '{0}' reads".format(side))
                partition_reads(phase_labels, it, side, 
                                   position_path, cut_cons_align, 
                                   polished_template, read_align, 
                                   cut_consensus, confirmed_pos_path, 
                                   partitioning, all_edge_headers[edge])
                print "Here 21"
                #7e. Write stats file for current iteration
                phase_pairs = sorted(combinations(phase_labels, 
                                                 2))
                for phase_one, phase_two in phase_pairs:
                    cons_one = cut_consensus[(it, side, phase_one)]
                    cons_two = cut_consensus[(it, side, phase_two)]
                    if (not os.path.isfile(cons_one) or 
                        not os.path.isfile(cons_two)):
                        continue
                    cons_cons_file = cons_vs_cons.format(
                                            it, side, phase_one, 
                                            it, side, phase_two)
                    flye_aln.make_alignment(cons_two, [cons_one], 
                                            num_threads, orient_dir, 
                                            args.platform, cons_cons_file,
                                            reference_mode=True, sam_output=True)
                print "Here 22"
                side_stat_outputs = update_side_stats(
                                    phase_labels, it, side, 
                                    cut_cons_align, polished_template, 
                                    confirmed_pos_path.format(it, side), 
                                    partitioning.format(it, side), 
                                    prev_partitionings[side], 
                                    side_stats.format(side))
                phase_below_cov[side], dup_part[side] = side_stat_outputs
                side_it[side] = it
                print "Here 23"
            iter_pairs.append((side_it[side_labels[0]], 
                               side_it[side_labels[1]]))
            update_int_stats(edge, side_labels, phase_labels, side_it, cut_cons_align, 
                                polished_template, 
                                template_len,
                                confirmed_pos_path, int_confirmed_path, 
                                partitioning, integrated_stats)
            print "Here 24"
            if both_break:
                break
        print "Here 25"
        #8. Finalize stats files
        logger.debug("Writing stats files")
        for side in side_labels:
            finalize_side_stats(phase_labels, side_it[side], 
                                side, cut_cons_align, polished_template, 
                                cons_vs_cons, cut_consensus, 
                                confirmed_pos_path.format(side_it[side], 
                                                          side), 
                                partitioning.format(side_it[side], side), 
                                phase_below_cov[side],
                                dup_part[side], term_bool[side], 
                                side_stats.format(side))
        final_int_outputs = finalize_int_stats(edge, side_labels, phase_labels, 
                                               side_it, cut_cons_align, 
                                               polished_template, 
                                               template_len, cons_vs_cons, 
                                               cut_consensus, 
                                               int_confirmed_path, 
                                               partitioning, 
                                               integrated_stats, 
                                               phased_edge_path)
        phased, phased_edges, summ_vals = final_int_outputs
        print "Here 26"
        #9. Generate summary and phased edge file
        logger.debug("Generating summary and phased edge file")
        avg_div = 0.0
        both_phased_present = False
        if phased:
            res_inds = phase_labels
            for res_one, res_two in sorted(combinations(res_inds, 2)):
                res_one_path = phased_edge_path.format(edge, res_one)
                res_two_path = phased_edge_path.format(edge, res_two)
                if (os.path.isfile(res_one_path) and
                    os.path.isfile(res_two_path)):
                    both_phased_present = True
                    edge_phased = True
                    flye_aln.make_alignment(res_two_path, [res_one_path], 
                                        num_threads, orient_dir, 
                                        args.platform, 
                                        res_vs_res.format(edge, res_one, res_two),
                                        reference_mode=True, sam_output=True)
            if both_phased_present:
                avg_div = int_stats_postscript(edge, side_labels, phase_labels,  
                                               integrated_stats, 
                                               phased_edge_path, 
                                               res_vs_res)
        print "Here 27"
        if both_phased_present:
            phased_dict.update(phased_edges)
        summary_list.append((edge, template_len, 
                             avg_cov, summ_vals, avg_div, 
                             both_phased_present))
        """remove_unneeded_files(repeat_edges, rep, side_labels, side_it, 
                              orient_dir, template, extended, pol_temp_dir, 
                              pol_ext_dir, pre_phased_reads, 
                              pre_partitioning, pre_read_align, 
                              partitioning, cons_align, cut_cons_align, 
                              read_align, confirmed_pos_path, phased_reads, 
                              cut_cons, polishing_dir, cons_vs_cons, 
                              int_confirmed_path, repeat_reads, 
                              frequency_path, alignment_file, 
                              NUM_POL_ITERS, iter_pairs)"""
    if edge_phased:
        logger.info("Edge successfully phased")
    else:
        logger.info("Edge not phased")
    return phased_dict, summary_list

def _find_div_region(pos_file, pol_template_file, alignment_file, 
                     left_part_file, right_part_file, window_size, phase_labels, 
                     left_reads_file, right_reads_file, initial_reads_file, 
                     headers_to_id):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    PHASE_MIN_WIN_DIV = trestle_config.vals["phase_min_win_div"]
    
    template = ""
    if os.path.getsize(pol_template_file):
        template_dict = fp.read_sequence_dict(pol_template_file)
        template = template_dict.values()[0]
    
    pos_headers, positions = div.read_positions(pos_file)
    all_pos = positions["total"]
    num_windows = ((len(template) - 1) / window_size) + 1
    win_count = [0 for _ in range(num_windows)]
    curr_win = 0
    high_win = -1
    high_count = 0
    for p in sorted(all_pos):
        while p >= (curr_win + 1) * window_size:
            curr_win += 1
        if curr_win < len(win_count):
            win_count[curr_win] += 1
            if win_count[curr_win] >= high_count:
                high_win = curr_win
                high_count = win_count[curr_win]
    
    win_start = high_win * window_size
    win_end = (high_win + 1) * window_size
    if high_count < PHASE_MIN_WIN_DIV * window_size:
        _write_partitioning_file([], left_part_file)
        _write_partitioning_file([], right_part_file)
    else:
        aligns = _read_alignment(alignment_file, pol_template_file,
                                 CONS_ALN_RATE)
        if aligns and aligns[0]:
            win_output = _all_reads_in_win(aligns[0], high_win, window_size)
            win_headers, win_dists, left_reads, right_reads = win_output
            
            head_ids = []
            headers = []
            distances = []
            for i, (h, d) in enumerate(sorted(zip(win_headers, win_dists))):
                head_ids.append(i)
                headers.append(h)
                distances.append([d])
            
            labels = []
            if len(distances) >= 2:
                clustering = AgglomerativeClustering(affinity='euclidean', 
                                                n_clusters=2, 
                                                linkage='average')
                model = clustering.fit(distances)
                print model.labels_
                labels = model.labels_
            cluster_one = []
            cluster_two = []
            for i, lab in enumerate(labels):
                if lab == 0:
                    cluster_one.append(headers[i])
                elif lab == 1:
                    cluster_two.append(headers[i])
            
            left_parts = _make_part_list(cluster_one, cluster_two, left_reads, 
                                         phase_labels, headers_to_id)
            _write_partitioning_file(left_parts, left_part_file)
            write_side_reads(initial_reads_file, left_part_file, left_reads_file)
            right_parts = _make_part_list(cluster_one, cluster_two, right_reads, 
                                          phase_labels, headers_to_id)
            _write_partitioning_file(right_parts, right_part_file)
            write_side_reads(initial_reads_file, right_part_file, right_reads_file)
        else:
            raise Exception("Unreadable alignment: {0}".format(alignment_file))
    return win_start, win_end
    
    
def _all_reads_in_win(alns, high_win, win_size):
    #require that the read completely overlaps the window
    win_st = high_win * win_size
    win_end = (high_win + 1) * win_size
    """
    Reads and returns all lines in a sam file
    0    1    2      3    4      5     6    7      8      9   10    11   12
    qID tID qStart qEnd qStrand qLen tStart tEnd tStrand tLen qSeq tSeq errRate
    """
    win_headers = []
    win_dists = []
    left_reads = []
    right_reads = []
    all_reads = set()
    #Assumes all strands are forward strands
    for i, aln in enumerate(alns):
        q_id = aln.qry_id
        t_st = aln.trg_start
        t_end = aln.trg_end
        #t_strand = aln.trg_sign
        q_align = aln.qry_seq
        t_align = aln.trg_seq
        if q_id in all_reads:
            continue
        if win_st >= t_st and t_end >= win_end:
            #Get read segments for each window
            trg_aln, aln_trg = _index_mapping(aln.trg_seq)
            al_st = trg_aln[win_st - t_st]
            al_end = trg_aln[win_end - t_st]
            matches = 0
            for q, t in zip(q_align[al_st:al_end], t_align[al_st:al_end]):
                if q == t:
                    matches += 1
            win_dists.append(1 - matches/float(al_end - al_st))
            win_headers.append(q_id)
        elif t_end < win_end:
            left_reads.append(q_id)
        elif win_st < t_st:
            right_reads.append(q_id)
        all_reads.add(q_id)
        
    return win_headers, win_dists, left_reads, right_reads

def _make_part_list(cluster_one, cluster_two, inner_reads, phase_labels, 
                    headers_to_id):
    part_list = []
    for h in cluster_one:
        read_id = headers_to_id[h]
        part_list.append((read_id, "Partitioned", phase_labels[0], 1, 0, h))
    for h in cluster_two:
        read_id = headers_to_id[h]
        part_list.append((read_id, "Partitioned", phase_labels[1], 1, 0, h))
    #Note that the qid is saved for the inner reads but not the clusters
    for h in inner_reads:
        read_id = headers_to_id[h]
        part_list.append((read_id, "None", "NA", 0, 0, h))
    return part_list

def write_side_reads(all_reads, partitioning, out_file):
    all_reads_dict = fp.read_sequence_dict(all_reads)
    part_list = _read_partitioning_file(partitioning)
    side_reads = {}
    for part in part_list:
        read_id, status, phase, top_sc, total_sc, header = part
        seq = all_reads_dict[header]
        side_reads[header] = seq
    if side_reads and side_reads.values()[0]:
        fp.write_fasta_dict(side_reads, out_file)
        

def define_file_names():
    #Defining directory and file names for phasing output
    pol_temp_name = "Polishing.Template"
    pol_cons_name = "Polishing.Consensus.{0}.{1}.{2}"
    pol_dir_names = pol_temp_name, pol_cons_name
    
    edge_label = "edge_{0}"
    template_name = "template.fasta"
    initial_reads_name = "initial_reads.fasta"
    left_reads_name = "left_reads.fasta"
    right_reads_name = "right_reads.fasta"
    initial_file_names = (edge_label, template_name, initial_reads_name, 
                          left_reads_name, right_reads_name)
    
    pre_phased_reads_name = "pre_phased_reads.{0}.{1}.txt"
    pre_read_aln_name = "pre_phased_reads.{0}.{1}.vs.extended.minimap.sam"
    pre_partitioning_name = "pre_partitioning.{0}.{1}.txt"
    partitioning_name = "partitioning.{0}.{1}.txt"
    pre_file_names = (pre_phased_reads_name, pre_read_aln_name, 
                      pre_partitioning_name, partitioning_name)
    
    div_freq_name = "divergence_frequencies.txt"
    div_pos_name = "divergent_positions.txt"
    div_summ_name = "divergence_summary.txt"
    div_file_names = div_freq_name, div_pos_name, div_summ_name
    
    reads_template_aln_name = "reads.vs.template.minimap.sam"
    cons_temp_aln_name = "uncut_consensus.{0}.{1}.{2}.vs.template.minimap.sam"
    cut_cons_temp_aln_name = "consensus.{0}.{1}.{2}.vs.template.minimap.sam"
    reads_cons_aln_name = "reads.vs.consensus.{0}.{1}.{2}.minimap.sam"
    aln_names = (reads_template_aln_name, cons_temp_aln_name, 
                 cut_cons_temp_aln_name, reads_cons_aln_name)
    
    confirmed_pos_name = "confirmed_positions.{0}.{1}.txt"
    phased_reads_name = "phased_reads.{0}.{1}.{2}.fasta"
    cut_cons_name = "consensus.{0}.{1}.{2}.fasta"
    cons_vs_cons_name = "".join(["consensus.{0}.{1}.{2}.vs.",
                                 "consensus.{3}.{4}.{5}.minimap.sam"])
    middle_file_names = (confirmed_pos_name, phased_reads_name, 
                         cut_cons_name, cons_vs_cons_name)
                         
    side_stats_name = "stats_{0}.txt"
    int_stats_name = "stats_integrated.txt"
    int_confirmed_pos_name = "integrated_confirmed_positions.{0}.{1}.txt"
    phased_edge_name = "phased_edge_{0}.copy.{1}.fasta"
    res_vs_res_name = "phased_edge_{0}.copy.{1}.vs.{2}.minimap.sam"
    output_file_names = (side_stats_name, int_stats_name, 
                         int_confirmed_pos_name, phased_edge_name, 
                         res_vs_res_name)
                         
    all_file_names = (pol_dir_names, initial_file_names, 
                      pre_file_names, div_file_names, aln_names, 
                      middle_file_names, output_file_names)
    return all_file_names


#Process Repeats functions
class ProcessingException(Exception):
    pass


def process_repeats(reads, uniques_dict, work_dir, initial_file_names):
    """Generates unique dirs and files given reads, uniques_dict and
    graph_edges files."""
    ORIENT_CONFIG = trestle_config.vals["orientations_to_run"]
    
    (edge_label, template_name, initial_reads_name, 
     left_reads_name, right_reads_name) = initial_file_names
    
    #Reads input files
    #repeats_dict = _read_repeats_dump(repeats_dump)
    if not uniques_dict:
        #logger.debug("Empty repeats_dump file: {0}".format(repeats_dump))
        return [], {}, {}
        
    reads_dict = {}
    for read_file in reads:
        reads_dict.update(fp.read_sequence_dict(read_file))
    #orig_graph = fp.read_sequence_dict(graph_edges)
    #graph_dict = {int(h.split('_')[1]):orig_graph[h] for h in orig_graph}
    
    if not reads_dict:
        raise ProcessingException("No reads found from {0}".format(reads))
    #if not graph_dict:
    #    raise ProcessingException("No edges found from {0}".format(
    #        graph_edges))
        
    edge_list = []
    all_edge_headers = {}
    for edge in sorted(uniques_dict, reverse=True):
        #Checks presence of reverse strand
        #One run processes both forward and reverse strand of edge
                
        if edge <= 0:
            continue
        
        valid_edge = True
        if -edge not in uniques_dict:
            logger.debug("Edge {0} missing reverse strand".format(edge))
            valid_edge = False
        if not valid_edge:
            continue
        
        #Makes edge dirs
        edge_dir = os.path.join(work_dir, edge_label.format(edge))
        if not os.path.isdir(edge_dir):
            os.mkdir(edge_dir)
        edge_list.append(edge)
        
        run_orientations = []
        if ORIENT_CONFIG == "forward":
            run_orientations = [("forward", edge)]
        elif ORIENT_CONFIG == "reverse":
            run_orientations = [("reverse", -edge)]
        elif ORIENT_CONFIG == "both":
            run_orientations = [("forward", edge), ("reverse", -edge)]
        for curr_label, curr_edge in run_orientations:
            orient_path = os.path.join(edge_dir, curr_label)
            if not os.path.isdir(orient_path):
                os.mkdir(orient_path)
            template_path = os.path.join(orient_path, template_name)
            initial_reads_path = os.path.join(orient_path, initial_reads_name)
                        
            #(mult, all_reads_list, inputs_dict,
            # outputs_dict) = repeats_dict[curr_rep]
            #mult = repeats_dict[curr_rep].multiplicity
            all_reads_list = uniques_dict[curr_edge].all_reads
            
            template_dict = {}
            initial_reads_dict = {}
            read_id = 0
            
            template_seq = uniques_dict[curr_edge].sequences["template"]
            #if curr_label == "reverse":
            #    template_seq = fp.reverse_complement(graph_dict[edge])
            template_dict[curr_edge] = template_seq
            
            all_edge_headers[curr_edge] = {}
            #rev_comp of read will be written if the header is -h
            
            for header in all_reads_list:
                if (not header) or (header[0] != '+' and header[0] != '-'):
                    raise ProcessingException(
                        "All reads format not recognized: {0}".format(header))
                if header[1:] not in reads_dict:
                    raise ProcessingException(
                        "Read header {0} not in any of {1}".format(
                            header[1:], reads))
                
                seq = reads_dict[header[1:]]
                if header[0] == '-':
                    seq = fp.reverse_complement(seq)
                initial_reads_dict[header[1:]] = seq
                
                if header[1:] not in all_edge_headers[curr_edge]:
                    all_edge_headers[curr_edge][header[1:]] = read_id
                    read_id += 1
            
            if template_dict and template_dict.values()[0]:
                fp.write_fasta_dict(template_dict, template_path)
            if initial_reads_dict and initial_reads_dict.values()[0]:
                fp.write_fasta_dict(initial_reads_dict, initial_reads_path)
            
            if not template_dict:
                raise ProcessingException("No template {0} found".format(
                                                curr_edge))
            if not initial_reads_dict:
                raise ProcessingException("No repeat reads {0} found".format(
                                                curr_edge))
    return edge_list, all_edge_headers


def _read_uniques_dump(uniques_dump):
    """Read uniques_dump.txt file following format guidelines"""
    uniques_dict = {}
    
    curr_repeat = 0
    all_reads_list = []
    all_bool = False
    with open(uniques_dump, 'r') as uf:
        for i, line in enumerate(uf):
            line = line.strip()
            if line:
                if line[0] == '#':
                    all_bool = False
                    parts = line.split('\t')
                    header = parts[0]
                    head_parts = header.split(' ')
                    if head_parts[0] == '#Unique':
                        if curr_repeat != 0:
                            uniques_dict[curr_repeat] = ("", all_reads_list)
                            all_reads_list = []
                        
                        curr_repeat = int(head_parts[1])
                    elif head_parts[0] == '#All':
                        all_bool = True
                    else:
                        raise ProcessingException(
                        "Unexpected repeats dump file format")
                else:
                    if all_bool:
                        all_reads_list.append(line)
                    else:
                        raise ProcessingException(
                        "Unexpected repeats dump file format")
        if curr_repeat != 0:
            uniques_dict[curr_repeat] = ("", all_reads_list)
    return uniques_dict
        

def _write_partitioning_file(part_list, part_path):
    with open(part_path, "w") as f:
        header_labels = ["Read_ID", "Status", "Phase", "Top Score", 
                         "Total Score", "Header"]
        spaced_header = map("{:11}".format, header_labels)
        f.write("\t".join(spaced_header))
        f.write("\n")
        for read_label in sorted(part_list):
            spaced_label = map("{:11}".format, read_label)
            f.write("\t".join(spaced_label))
            f.write("\n")


def _read_partitioning_file(partitioning_file):
    part_list = []
    with open(partitioning_file, "r") as f:
        for i, line in enumerate(f):
            if i > 0:
                line = line.strip()
                tokens = [t.strip() for t in line.split("\t")]
                for int_ind in [0, 3, 4]:
                    tokens[int_ind] = int(tokens[int_ind])
                part_list.append(tuple(tokens))
    return part_list
    

def find_coverage(frequency_file):
    coverage = 0.0
    if os.path.isfile(frequency_file):
        header, freqs = div.read_frequency_path(frequency_file)
        cov_ind = header.index("Cov")
        all_covs = [f[cov_ind] for f in freqs]
        coverage = _mean(all_covs)
        #print min(all_covs), _mean(all_covs), max(all_covs)
    return coverage


def write_phased_reads(it, side, phase_id, all_reads, partitioning, out_file):
    all_reads_dict = fp.read_sequence_dict(all_reads)
    part_list = _read_partitioning_file(partitioning)
    phase_header_name = "Read_{0}|Iter_{1}|Side_{2}|Phase_{3}|{4}"
    phased_reads = {}
    for part in part_list:
        read_id, status, phase, top_sc, total_sc, header = part
        if status == "Partitioned" and phase != "NA" and int(phase) == phase_id:
            phase_seq = all_reads_dict[header]
            phase_header = phase_header_name.format(read_id, it, 
                                                  side, phase_id, header)
            phased_reads[phase_header] = phase_seq
    if phased_reads and phased_reads.values()[0]:
        fp.write_fasta_dict(phased_reads, out_file)


def init_partitioning(phases, side, pre_partitioning, pre_read_align, extended, 
                      partitioning):
    FLANKING_LEN = trestle_config.vals["flanking_len"]
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    #dict from read_header to phase
    extend_overlap_reads = {}
    for phase in phases:
        non_overlap_reads = 0
        aligns = _read_alignment(pre_read_align.format(side, phase), 
                                 extended[(side, phase)], CONS_ALN_RATE)
        if aligns and aligns[0]:
            for aln in aligns[0]:
                phase_header = aln.qry_id
                read_header = phase_header.split("|")[-1]
                if ((side == "in" and 
                        aln.trg_start < FLANKING_LEN) or 
                    (side == "out" and 
                        aln.trg_end >= aln.trg_len - FLANKING_LEN)):
                    extend_overlap_reads[read_header] = str(phase)
                else:
                    non_overlap_reads += 1
        logger.debug("Side {0}, phase {1}, non-overlap reads = {2}".format(
                        side, phase, non_overlap_reads))
    partitioned_reads = []
    part_list = _read_partitioning_file(pre_partitioning)
    for part in part_list:
        read_id, status, phase, top_sc, total_sc, header = part
        if header in extend_overlap_reads:
            partitioned_reads.append((read_id, "Partitioned", 
                                      extend_overlap_reads[header],
                                      1, 0, header))
        else:
            partitioned_reads.append((read_id, "None", "NA", 0, 0, header))
    _write_partitioning_file(partitioned_reads, partitioning)


#Cut Consensus Functions
def find_read_endpoints(alignment_file, template):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    read_endpoints = {}
    aligns = _read_alignment(alignment_file, template, CONS_ALN_RATE)
    if aligns and aligns[0]:
        for aln in aligns[0]:
            read_header = aln.qry_id
            start = aln.trg_start
            end = aln.trg_end
            if read_header not in read_endpoints:
                read_endpoints[read_header] = (start, end)
    else:
        logger.debug("No read alignment to template, no read_endpoints")
    return read_endpoints


def locate_consensus_cutpoint(side, read_endpoints, phased_read_file):
    MIN_EDGE_COV = trestle_config.vals["min_edge_cov"]
    all_endpoints = []
    max_endpoint = 0
    phased_reads = fp.read_sequence_dict(phased_read_file)
    for phase_header in phased_reads:
        parts = phase_header.split("|")
        read_header = parts[-1]
        if read_header in read_endpoints:
            endpoint = read_endpoints[read_header]
            if max(endpoint) > max_endpoint:
                max_endpoint = max(endpoint)
            all_endpoints.append(endpoint)
    coverage = [0 for _ in range(max_endpoint + 1)]
    for start, end in all_endpoints:
        for x in range(start, end):
            coverage[x] += 1
    window_len = 100
    cutpoint = -1
    for i in range(len(coverage) - window_len):
        if side == "right":
            window_start = (len(coverage) - window_len) - i
            window_end = len(coverage) - i
            if window_start < 0:
                window_start = 0
            if window_end > len(coverage):
                window_end = len(coverage)
            avg_cov = _mean(coverage[window_start:window_end])
            if avg_cov >= MIN_EDGE_COV:
                cutpoint = window_end
                break
        elif side == "left":
            window_start = i
            window_end = i + window_len
            if window_start < 0:
                window_start = 0
            if window_end > len(coverage):
                window_end = len(coverage)
            avg_cov = _mean(coverage[window_start:window_end])
            if avg_cov >= MIN_EDGE_COV:
                cutpoint = window_start
                break
    return cutpoint


def truncate_consensus(side, cutpoint, cons_al_file, template,
                       polished_consensus, cut_cons_file, win_cutpoint):
    if cutpoint == -1:
        logger.debug("No cutpoint for consensus file")
        return
        
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    cons_al = _read_alignment(cons_al_file, template, CONS_ALN_RATE)
    consensus_endpoint = -1
    if cons_al and cons_al[0]:
        consensus_endpoint = _find_consensus_endpoint(cutpoint, cons_al, side,
                                                      False)
        win_endpoint = _find_consensus_endpoint(win_cutpoint, cons_al, side,
                                                      True)
    else:
        logger.debug("No cons alignment to template, no cut consensus")
        return
    
    if consensus_endpoint != -1 and win_endpoint != -1:
        cons_seqs = fp.read_sequence_dict(polished_consensus)
        cons_head = cons_seqs.keys()[0]
        consensus = cons_seqs.values()[0]
        if side == "right":
            start = 0
            if win_endpoint >= 0:
                start = win_endpoint
            end = len(consensus)
            if consensus_endpoint < len(consensus):
                end = consensus_endpoint
        elif side == "left":
            start = 0
            if consensus_endpoint >= 0:
                start = consensus_endpoint
            end = len(consensus)
            if win_endpoint < len(consensus):
                end = win_endpoint
        cut_head = "".join([cons_head, "|{0}_{1}".format(start, end)])
        cut_dict = {cut_head:consensus[start:end]}
        fp.write_fasta_dict(cut_dict, cut_cons_file)


def _find_consensus_endpoint(cutpoint, aligns, side, win_bool):
    consensus_endpoint = -1
    #first try collapsing
    coll_aln = _collapse_cons_aln(aligns)
    print "Here z", cutpoint, win_bool, coll_aln.trg_end
    if cutpoint >= coll_aln.trg_start and cutpoint < coll_aln.trg_end:
        trg_aln, aln_trg = _index_mapping(coll_aln.trg_seq)
        qry_aln, aln_qry = _index_mapping(coll_aln.qry_seq)
        cutpoint_minus_start = cutpoint - coll_aln.trg_start
        print "Here a", cutpoint_minus_start, len(trg_aln)
        aln_ind = trg_aln[cutpoint_minus_start]
        print "Here a", aln_ind, len(aln_qry)
        qry_ind = aln_qry[aln_ind]
        consensus_endpoint = qry_ind + coll_aln.qry_start
    elif cutpoint == coll_aln.trg_end:
        consensus_endpoint = coll_aln.qry_end
    else:
        #otherwise try each alignment
        MIN_SUPP_ALN_LEN = trestle_config.vals["min_supp_align_len"]
        #save tuples of cutpoint distance, cutpoint
        aln_endpoints = []
        for i, aln in enumerate(aligns[0]):
            if i == 0 or len(aln.trg_seq) >= MIN_SUPP_ALN_LEN:
                if cutpoint >= aln.trg_start and cutpoint < aln.trg_end:
                    trg_aln, aln_trg = _index_mapping(aln.trg_seq)
                    qry_aln, aln_qry = _index_mapping(aln.qry_seq)
                    cutpoint_minus_start = cutpoint - aln.trg_start
                    if cutpoint_minus_start < 0:
                        #print aln.qry_id, aln.trg_id, side, cutpoint, cutpoint_minus_start
                        aln_ind = trg_aln[0]
                    elif cutpoint_minus_start >= len(trg_aln):
                        #print aln.qry_id, aln.trg_id, side, cutpoint, cutpoint_minus_start
                        aln_ind = trg_aln[-1]
                    else:
                        aln_ind = trg_aln[cutpoint_minus_start]
                    qry_ind = aln_qry[aln_ind]
                    endpoint = qry_ind + coll_aln.qry_start
                    aln_endpoints.append((0, endpoint))
                elif side == "right" and cutpoint >= aln.trg_end and not win_bool:
                    endpoint = aln.qry_end
                    distance = cutpoint - aln.trg_end
                    aln_endpoints.append((distance, endpoint))
                elif side == "left" and cutpoint < aln.trg_start and not win_bool:
                    endpoint = aln.qry_start
                    distance = aln.trg_start - cutpoint
                    aln_endpoints.append((distance, endpoint))
                elif side == "right" and cutpoint < aln.trg_start and win_bool:
                    endpoint = aln.qry_start
                    distance = aln.trg_start - cutpoint
                    aln_endpoints.append((distance, endpoint))
                elif side == "left" and cutpoint >= aln.trg_end and win_bool:
                    endpoint = aln.qry_end
                    distance = cutpoint - aln.trg_end
                    aln_endpoints.append((distance, endpoint))
        if aln_endpoints:
            consensus_endpoint = sorted(aln_endpoints)[0][1]
    return consensus_endpoint

#Partition Reads Functions


def partition_reads(phases, it, side, position_path, cons_align_path, 
                    template, read_align_path, consensuses, 
                    confirmed_pos_path, part_file, 
                    headers_to_id):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    BUFFER_COUNT = trestle_config.vals["buffer_count"]
    
    skip_bool = False
    pos_headers, pos = div.read_positions(position_path)
        
    cons_aligns = {}
    for phase_id in phases:
        if not os.path.isfile(cons_align_path.format(it, side, phase_id)):
            skip_bool = True
        else:
            cons_aligns[phase_id] = _read_alignment(cons_align_path.format(it, 
                                                        side, 
                                                        phase_id), 
                                                   template, 
                                                   CONS_ALN_RATE)
        if (skip_bool or 
                not cons_aligns or
                not cons_aligns[phase_id] or 
                not cons_aligns[phase_id][0]):
            logger.debug("No cons alignment found for phase {0}".format(
                            phase_id))
            skip_bool = True
    if skip_bool:
        if it <= 1:
            confirmed_pos = {"total":[], "sub":[], "ins":[], "del":[]}
            rejected_pos = {"total":[], "sub":[], "ins":[], "del":[]}
            consensus_pos = pos
        else:
            previous_pos = _read_confirmed_positions(
                                confirmed_pos_path.format(it - 1, side))
            confirmed_pos, rejected_pos, consensus_pos = previous_pos
    else:
        #collapse multiple consensus alignments to the template
        coll_cons_aligns = {}
        for phase_id in cons_aligns:
            aln = cons_aligns[phase_id]
            coll_cons_aligns[phase_id] = _collapse_cons_aln(aln)
        curr_pos = _evaluate_positions(pos, coll_cons_aligns, side)
        confirmed_pos, rejected_pos, consensus_pos = curr_pos
    _write_confirmed_positions(confirmed_pos, rejected_pos, pos, 
                               confirmed_pos_path.format(it, side))
    read_aligns = {}
    for phase_id in phases:
        if (not os.path.isfile(read_align_path.format(it, side, phase_id)) or
            not os.path.isfile(consensuses[(it, side, phase_id)])):
            skip_bool = True
        elif not skip_bool:
            read_aligns[phase_id] = _read_alignment(
                                        read_align_path.format(it, side, 
                                                               phase_id), 
                                        consensuses[(it, side, phase_id)], 
                                        CONS_ALN_RATE)
        if (skip_bool or 
                not read_aligns or 
                not read_aligns[phase_id] or 
                not read_aligns[phase_id][0]):
            read_aln_str = "No read alignment found for phase {0}"
            logger.debug(read_aln_str.format(phase_id))
            skip_bool = True
    if skip_bool:
        partitioning = _read_partitioning_file(part_file.format(it - 1, side))
    else:
        partitioning = _classify_reads(read_aligns, consensus_pos, 
                        headers_to_id, BUFFER_COUNT)
    _write_partitioning_file(partitioning, part_file.format(it, side))
    

def _read_alignment(alignment, target_path, min_aln_rate):
    alignments = []
    aln_reader = flye_aln.SynchronizedSamReader(alignment,
                                       fp.read_sequence_dict(target_path),
                                       config.vals["max_read_coverage"])
    aln_reader.init_reading()
    while not aln_reader.is_eof():
        ctg_id, ctg_aln = aln_reader.get_chunk()
        if ctg_id is None:
            break
        alignments.append(ctg_aln)
    
    return alignments


def _collapse_cons_aln(cons_aligns):
    MAX_SUPP_ALIGN_OVERLAP = trestle_config.vals["max_supp_align_overlap"]
    coll_aln = None
    for aln in cons_aligns[0]:
        if coll_aln is None:
            coll_aln = aln
        elif _overlap(coll_aln, aln) <= MAX_SUPP_ALIGN_OVERLAP:
            coll_aln = _collapse(coll_aln, aln)
    return coll_aln


def _overlap(aln_one, aln_two):
    qry_overlap_lens = []
    if (aln_one.qry_start >= aln_two.qry_start and 
            aln_one.qry_start < aln_two.qry_end):
        if aln_one.qry_end >= aln_two.qry_end:
            qry_overlap_lens.append(aln_two.qry_end - aln_one.qry_start)
        else:
            qry_overlap_lens.append(aln_one.qry_end - aln_one.qry_start)
    if (aln_one.qry_end > aln_two.qry_start and
            aln_one.qry_end <= aln_two.qry_end):
        if aln_one.qry_start <= aln_two.qry_start:
            qry_overlap_lens.append(aln_one.qry_end - aln_two.qry_start)
        else:
            qry_overlap_lens.append(aln_one.qry_end - aln_one.qry_start)
    if (aln_two.qry_start >= aln_one.qry_start and
            aln_two.qry_start < aln_one.qry_end):
        if aln_two.qry_end >= aln_one.qry_end:
            qry_overlap_lens.append(aln_one.qry_end - aln_two.qry_start)
        else:
            qry_overlap_lens.append(aln_two.qry_end - aln_two.qry_start)
    if (aln_two.qry_end > aln_one.qry_start and
            aln_two.qry_end <= aln_one.qry_end):
        if aln_two.qry_start <= aln_one.qry_start:
            qry_overlap_lens.append(aln_two.qry_end - aln_one.qry_start)
        else:
            qry_overlap_lens.append(aln_two.qry_end - aln_two.qry_start)
    qry_len = 0
    if qry_overlap_lens:
        qry_len = min(qry_overlap_lens)
    trg_overlap_lens = []
    if (aln_one.trg_start >= aln_two.trg_start and 
            aln_one.trg_start < aln_two.trg_end):
        if aln_one.trg_end >= aln_two.trg_end:
            trg_overlap_lens.append(aln_two.trg_end - aln_one.trg_start)
        else:
            trg_overlap_lens.append(aln_one.trg_end - aln_one.trg_start)
    if (aln_one.trg_end > aln_two.trg_start and
            aln_one.trg_end <= aln_two.trg_end):
        if aln_one.trg_start <= aln_two.trg_start:
            trg_overlap_lens.append(aln_one.trg_end - aln_two.trg_start)
        else:
            trg_overlap_lens.append(aln_one.trg_end - aln_one.trg_start)
    if (aln_two.trg_start >= aln_one.trg_start and
            aln_two.trg_start < aln_one.trg_end):
        if aln_two.trg_end >= aln_one.trg_end:
            trg_overlap_lens.append(aln_one.trg_end - aln_two.trg_start)
        else:
            trg_overlap_lens.append(aln_two.trg_end - aln_two.trg_start)
    if (aln_two.trg_end > aln_one.trg_start and
            aln_two.trg_end <= aln_one.trg_end):
        if aln_two.trg_start <= aln_one.trg_start:
            trg_overlap_lens.append(aln_two.trg_end - aln_one.trg_start)
        else:
            trg_overlap_lens.append(aln_two.trg_end - aln_two.trg_start)
    trg_len = 0
    if trg_overlap_lens:
        trg_len = min(trg_overlap_lens)
    return max([qry_len, trg_len])


def _collapse(aln_one, aln_two):
    MAX_SUPP_ALIGN_OVERLAP = trestle_config.vals["max_supp_align_overlap"]
    out_aln = copy.deepcopy(aln_one)
    if (aln_one.qry_sign == "-" or aln_two.qry_sign == "-" or 
            _overlap(aln_one, aln_two) > MAX_SUPP_ALIGN_OVERLAP):
        return out_aln
    if (aln_one.qry_start <= aln_two.qry_start and 
            aln_one.trg_start <= aln_two.trg_start):
        qry_merge_outs = _merge_alns(aln_one.qry_start, aln_one.qry_end, 
                                     aln_one.qry_seq, aln_two.qry_start, 
                                     aln_two.qry_end, aln_two.qry_seq)
        one_qry_seq, two_qry_seq, out_qry_end = qry_merge_outs
        trg_merge_outs = _merge_alns(aln_one.trg_start, aln_one.trg_end, 
                                     aln_one.trg_seq, aln_two.trg_start, 
                                     aln_two.trg_end, aln_two.trg_seq)
        one_trg_seq, two_trg_seq, out_trg_end = trg_merge_outs
        fill_qry = ""
        fill_trg = ""
        qry_lens = len(one_qry_seq) + len(two_qry_seq)
        trg_lens = len(one_trg_seq) + len(two_trg_seq)
        if qry_lens > trg_lens:
            diff = qry_lens - trg_lens
            fill_trg = "-" * diff
        elif trg_lens > qry_lens:
            diff = trg_lens - qry_lens
            fill_qry = "-" * diff
        out_qry_seq = "".join([one_qry_seq, fill_qry, two_qry_seq])
        out_trg_seq = "".join([one_trg_seq, fill_trg, two_trg_seq])
        out_err_rate = ((aln_one.err_rate * len(aln_one.trg_seq) + 
                         aln_two.err_rate * len(aln_two.trg_seq)) /
                         (len(aln_one.trg_seq) + len(aln_two.trg_seq)))
        out_aln = Alignment(aln_one.qry_id, aln_one.trg_id, aln_one.qry_start, 
                            out_qry_end, aln_one.qry_sign, aln_one.qry_len, 
                            aln_one.trg_start, out_trg_end, aln_one.trg_sign, 
                            aln_one.trg_len, out_qry_seq, out_trg_seq, 
                            out_err_rate)
        return out_aln
    elif (aln_two.qry_start <= aln_one.qry_start and 
            aln_two.trg_start <= aln_one.trg_start):
        qry_merge_outs = _merge_alns(aln_two.qry_start, aln_two.qry_end, 
                                     aln_two.qry_seq, aln_one.qry_start, 
                                     aln_one.qry_end, aln_one.qry_seq)
        two_qry_seq, one_qry_seq, out_qry_end = qry_merge_outs
        trg_merge_outs = _merge_alns(aln_two.trg_start, aln_two.trg_end, 
                                     aln_two.trg_seq, aln_one.trg_start, 
                                     aln_one.trg_end, aln_one.trg_seq)
        two_trg_seq, one_trg_seq, out_trg_end = trg_merge_outs
        fill_qry = ""
        fill_trg = ""
        qry_lens = len(two_qry_seq) + len(one_qry_seq)
        trg_lens = len(two_trg_seq) + len(one_trg_seq)
        if qry_lens > trg_lens:
            diff = qry_lens - trg_lens
            fill_trg = "-" * diff
        elif trg_lens > qry_lens:
            diff = trg_lens - qry_lens
            fill_qry = "-" * diff
        out_qry_seq = "".join([two_qry_seq, fill_qry, one_qry_seq])
        out_trg_seq = "".join([two_trg_seq, fill_trg, one_trg_seq])
        out_err_rate = ((aln_one.err_rate * len(aln_one.trg_seq) + 
                         aln_two.err_rate * len(aln_two.trg_seq)) /
                         (len(aln_one.trg_seq) + len(aln_two.trg_seq)))
        out_aln = Alignment(aln_one.qry_id, aln_one.trg_id, aln_two.qry_start, 
                            out_qry_end, aln_one.qry_sign, aln_one.qry_len, 
                            aln_two.trg_start, out_trg_end, aln_one.trg_sign, 
                            aln_one.trg_len, out_qry_seq, out_trg_seq, 
                            out_err_rate)
        return out_aln
    return out_aln


def _merge_alns(first_start, first_end, first_seq, 
                second_start, second_end, second_seq):
    first_out_seq = first_seq
    second_out_seq = second_seq
    out_end = second_end
    if first_end <= second_start:
        fill_qry_seq = "N" * (second_start - first_end)
        first_out_seq = "".join([first_seq, fill_qry_seq])
        second_out_seq = second_seq
    else:
        if first_end < second_end:
            overlap = first_end - second_start
            two_cut_ind = _overlap_to_aln_ind(overlap, second_seq)
            first_out_seq = first_seq
            second_out_seq = second_seq[two_cut_ind:]
        else:
            first_out_seq = first_seq
            second_out_seq = ""
            out_end = first_end
    return first_out_seq, second_out_seq, out_end


def _overlap_to_aln_ind(overlap, aln):
    num_bases = 0
    for i, base in enumerate(aln):
        if base != "-":
            num_bases += 1
        if num_bases == overlap:
            return i + 1
    return len(aln)


class PhaseAlignment:
    __slots__ = ("phase_id", "qry_seq", "trg_seq", "qry_start", "trg_start", 
                 "trg_end", "in_alignment", "curr_aln_ind", "curr_qry_nuc", 
                 "curr_trg_nuc", "curr_ins_nuc")

    def __init__(self, phase_id, qry_seq, trg_seq, qry_start, trg_start, 
                 trg_end):
        self.phase_id = phase_id
        self.qry_seq = flye_aln.shift_gaps(trg_seq, qry_seq)
        self.trg_seq = flye_aln.shift_gaps(self.qry_seq, trg_seq)
        self.qry_start = qry_start
        self.trg_start = trg_start
        self.trg_end = trg_end
        self.in_alignment = False
        self.curr_aln_ind = -1
        self.curr_qry_ind = -1
        self.curr_qry_nuc = ""
        self.curr_trg_nuc = ""
        self.curr_ins_nuc = ""       
    
    def reset_nucs(self):
        self.curr_qry_nuc = ""
        self.curr_trg_nuc = ""
        self.curr_ins_nuc = ""


def _evaluate_positions(pos, cons_aligns, side):
    #Includes insertions!
    confirmed_pos = {"total":[], "sub":[], "ins":[], "del":[]}
    rejected_pos = {"total":[], "sub":[], "ins":[], "del":[]}
    consensus_pos = {e:[] for e in cons_aligns}    
    
    alns = {}
    for phase_id in cons_aligns:
        orig_aln = cons_aligns[phase_id]
        alns[phase_id] = PhaseAlignment(phase_id, orig_aln.qry_seq, 
                                      orig_aln.trg_seq, orig_aln.qry_start, 
                                      orig_aln.trg_start, orig_aln.trg_end)
    
    min_start_phase = min([alns[e].trg_start for e in alns])
    max_end_phase = max([alns[e].trg_end for e in alns])
    #end indices for conservatively defining confirmed positions
    min_end_phase = min([alns[e].trg_end for e in alns])
    max_start_phase = max([alns[e].trg_start for e in alns])
    
    for trg_ind in range(min_start_phase, max_end_phase):
        for phase_id in alns:
            aln = alns[phase_id]
            if aln.trg_start == trg_ind:
                aln.curr_aln_ind = 0
                aln.curr_qry_ind = aln.qry_start
                aln.in_alignment = True
            
            if aln.trg_start > trg_ind or aln.trg_end <= trg_ind:
                aln.in_alignment = False
            
            if aln.in_alignment:
                while aln.trg_seq[aln.curr_aln_ind] == "-":                    
                    if aln.qry_seq[aln.curr_aln_ind] != "-":
                        aln.curr_ins_nuc += aln.qry_seq[aln.curr_aln_ind]
                        aln.curr_qry_ind += 1
                    aln.curr_aln_ind += 1
                aln.curr_qry_nuc = aln.qry_seq[aln.curr_aln_ind]
                aln.curr_trg_nuc = aln.trg_seq[aln.curr_aln_ind]
                
        
        if trg_ind in pos["total"]:
            if ((side == "right" and trg_ind < min_end_phase) or
                (side == "left" and trg_ind >= max_start_phase)):
                ins_confirmed = False
                del_confirmed = False
                sub_confirmed = False
                qry_nuc = ""
                trg_nuc = ""
                for phase_id in alns:
                    aln = alns[phase_id]
                    if aln.in_alignment:
                        #Directly add positions only to consensuses 
                        # where insertions occur
                        #Add the position prior to curr_qry_ind to 
                        # account for insertion
                        if aln.curr_ins_nuc:
                            ins_confirmed = True
                            consensus_pos[phase_id].append(aln.curr_qry_ind - 1)
                        
                        if qry_nuc and qry_nuc != aln.curr_qry_nuc:
                            if qry_nuc != "N" and aln.curr_qry_nuc != "N":
                                if qry_nuc == "-":
                                    del_confirmed = True
                                else:
                                    sub_confirmed = True
                        else:
                            qry_nuc = aln.curr_qry_nuc
                        if (trg_nuc and trg_nuc != aln.curr_trg_nuc and
                                trg_nuc != "N" and aln.curr_trg_nuc != "N"):
                            incon_str = "".join(["Inconsistent trg_nuc,",
                                                 " {0} {1} {2} {3}"])
                            logger.debug(incon_str.format(phase_id, 
                                                      trg_ind, trg_nuc, 
                                                      aln.curr_trg_nuc))
                        trg_nuc = aln.curr_trg_nuc
                if ins_confirmed or del_confirmed or sub_confirmed:
                    confirmed_pos["total"].append(trg_ind)
                    #Add positions to consensuses for only subs/deletions
                    if del_confirmed or sub_confirmed:
                        for phase_id in alns:
                            aln = alns[phase_id]
                            if aln.in_alignment:
                                consensus_pos[phase_id].append(aln.curr_qry_ind)
                    if trg_ind in pos["ins"]:
                        if ins_confirmed:
                            confirmed_pos["ins"].append(trg_ind)
                        else:
                            rejected_pos["ins"].append(trg_ind)
                    if trg_ind in pos["del"]:
                        if del_confirmed:
                            confirmed_pos["del"].append(trg_ind)
                        else:
                            rejected_pos["del"].append(trg_ind)
                    if trg_ind in pos["sub"]:
                        if sub_confirmed:
                            confirmed_pos["sub"].append(trg_ind)
                        else:
                            rejected_pos["sub"].append(trg_ind)
                else:
                    rejected_pos["total"].append(trg_ind)
                    if trg_ind in pos["ins"]:
                        rejected_pos["ins"].append(trg_ind)
                    if trg_ind in pos["del"]:
                        rejected_pos["del"].append(trg_ind)
                    if trg_ind in pos["sub"]:
                        rejected_pos["sub"].append(trg_ind)
        
        for phase_id in alns:
            aln = alns[phase_id]
            if aln.in_alignment:
                if aln.qry_seq[aln.curr_aln_ind] != "-":
                    aln.curr_qry_ind += 1
                aln.curr_aln_ind += 1
                
                aln.reset_nucs()
    
    return confirmed_pos, rejected_pos, consensus_pos


def _write_confirmed_positions(confirmed, rejected, pos, out_file):
    with open(out_file, 'w') as f:
        f.write(">Confirmed_total_positions_{0}\n".format(
                        len(confirmed["total"])))
        f.write(",".join(map(str, sorted(confirmed["total"]))))
        f.write("\n")
        f.write(">Confirmed_sub_positions_{0}\n".format(len(confirmed["sub"])))
        f.write(",".join(map(str, sorted(confirmed["sub"]))))
        f.write("\n")
        f.write(">Confirmed_del_positions_{0}\n".format(len(confirmed["del"])))
        f.write(",".join(map(str, sorted(confirmed["del"]))))
        f.write("\n")
        f.write(">Confirmed_ins_positions_{0}\n".format(len(confirmed["ins"])))
        f.write(",".join(map(str, sorted(confirmed["ins"]))))
        f.write("\n")
        f.write(">Rejected_total_positions_{0}\n".format(
                        len(rejected["total"])))
        f.write(",".join(map(str, sorted(rejected["total"]))))
        f.write("\n")
        f.write(">Rejected_sub_positions_{0}\n".format(len(rejected["sub"])))
        f.write(",".join(map(str, sorted(rejected["sub"]))))
        f.write("\n")
        f.write(">Rejected_del_positions_{0}\n".format(len(rejected["del"])))
        f.write(",".join(map(str, sorted(rejected["del"]))))
        f.write("\n")
        f.write(">Rejected_ins_positions_{0}\n".format(len(rejected["ins"])))
        f.write(",".join(map(str, sorted(rejected["ins"]))))
        f.write("\n")
        f.write(">Tentative_total_positions_{0}\n".format(len(pos["total"])))
        f.write(",".join(map(str, sorted(pos["total"]))))
        f.write("\n")
        f.write(">Tentative_sub_positions_{0}\n".format(len(pos["sub"])))
        f.write(",".join(map(str, sorted(pos["sub"]))))
        f.write("\n")
        f.write(">Tentative_del_positions_{0}\n".format(len(pos["del"])))
        f.write(",".join(map(str, sorted(pos["del"]))))
        f.write("\n")
        f.write(">Tentative_ins_positions_{0}\n".format(len(pos["ins"])))
        f.write(",".join(map(str, sorted(pos["ins"]))))
        f.write("\n")


def _read_confirmed_positions(confirmed_file):
    confirmed = {"total":[], "sub":[], "ins":[], "del":[]}
    rejected = {"total":[], "sub":[], "ins":[], "del":[]}
    pos = {"total":[], "sub":[], "ins":[], "del":[]}
    with open(confirmed_file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i == 1 and line:
                confirmed["total"] = map(int, line.split(","))
            elif i == 3 and line:
                confirmed["sub"] = map(int, line.split(","))
            elif i == 5 and line:
                confirmed["del"] = map(int, line.split(","))
            elif i == 7 and line:
                confirmed["ins"] = map(int, line.split(","))
            elif i == 9 and line:
                rejected["total"] = map(int, line.split(","))
            elif i == 11 and line:
                rejected["sub"] = map(int, line.split(","))
            elif i == 13 and line:
                rejected["del"] = map(int, line.split(","))
            elif i == 15 and line:
                rejected["ins"] = map(int, line.split(","))
            elif i == 17 and line:
                pos["total"] = map(int, line.split(","))
            elif i == 19 and line:
                pos["sub"] = map(int, line.split(","))
            elif i == 21 and line:
                pos["del"] = map(int, line.split(","))
            elif i == 23 and line:
                pos["ins"] = map(int, line.split(","))
    return confirmed, rejected, pos
                

def _classify_reads(read_aligns, consensus_pos, 
                    headers_to_id, buffer_count):
    #Includes insertion positions where an insertion occurs right before the 
    #position for the read.
    #partitioning format same as above:
    #list of (read_id, status, phase_id, top_score, total_score, header)
    partitioning = []
    
    read_scores = {}    
    for phase_id in read_aligns:
        read_counts = {}
        for aln in read_aligns[phase_id][0]:
            read_header = aln.qry_id
            cons_header = aln.trg_id
            #Unmapped segments will not be scored
            if cons_header == "*":
                continue
            if read_header not in read_scores:
                read_scores[read_header] = {}
            read_scores[read_header][phase_id] = 0
            if read_header not in read_counts:
                read_counts[read_header] = 1
            else:
                read_counts[read_header] += 1
            #Any alignments after the first supplementary will not be scored
            if read_counts[read_header] > 2:
                continue
            positions = consensus_pos[phase_id]
            trg_aln, aln_trg = _index_mapping(aln.trg_seq)
            for pos in positions:
                if pos >= aln.trg_start and pos < aln.trg_end:
                    pos_minus_start = pos - aln.trg_start
                    aln_ind = trg_aln[pos_minus_start]
                    if aln.qry_seq[aln_ind] == aln.trg_seq[aln_ind]:
                        read_scores[read_header][phase_id] += 1
    #Iterate through all read_headers so partitioning will be a complete set
    for read_header in headers_to_id:
        read_id = headers_to_id[read_header]
        if read_header in read_scores:
            tie_bool = False
            top_phase = 0
            top_score = 0
            total_score = 0
            for phase_id in read_scores[read_header]:
                phase_score = read_scores[read_header][phase_id]
                #print phase_id, phase_score, top_score
                if phase_score - buffer_count > top_score:
                    top_phase = phase_id
                    top_score = phase_score
                    tie_bool = False
                elif (phase_score - buffer_count <= top_score and 
                      phase_score >= top_score):
                    top_score = phase_score
                    tie_bool = True
                elif (phase_score >= top_score - buffer_count and
                      phase_score < top_score):
                    tie_bool = True
                total_score += phase_score
            
            if total_score == 0:
                status_label = "None"
                phase_label = "NA"
            elif tie_bool:
                status_label = "Tied"
                phase_label = "NA"
            else:
                status_label = "Partitioned"
                phase_label = str(top_phase)
            partitioning.append((read_id, status_label, phase_label, 
                                 top_score, total_score, read_header))
        else:
            status_label = "None"
            phase_label = "NA"
            top_score = 0
            total_score = 0
            partitioning.append((read_id, status_label, phase_label, 
                                top_score, total_score, read_header))    
    return partitioning


def _index_mapping(aln):
    #Given a genomic index, return the alignment index of the alignment
    al_inds = []
    #Given an alignment index, return the genomic index at that position
    gen_inds = []
    for i,b in enumerate(aln):
        gen_inds.append(len(al_inds))
        if b != '-':
            al_inds.append(i)
    return al_inds, gen_inds


def init_side_stats(edge, side, phase_labels, min_overlap, position_path, 
                    partitioning, prev_parts, template_len, stats_file):
    SUB_THRESH = trestle_config.vals["sub_thresh"]
    DEL_THRESH = trestle_config.vals["del_thresh"]
    INS_THRESH = trestle_config.vals["ins_thresh"]
    FLANKING_LEN = trestle_config.vals["flanking_len"]
    BUFFER_COUNT = trestle_config.vals["buffer_count"]
    MAX_ITER = trestle_config.vals["max_iter"]
    MIN_EDGE_COV = trestle_config.vals["min_edge_cov"]
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    
    pos_headers, pos = div.read_positions(position_path)
    #Count partitioned reads
    phase_below_cov = False
    part_list = _read_partitioning_file(partitioning)
    phased_reads, tied_reads, unassigned_reads = _get_partitioning_info(
                                                    part_list, 
                                                    phase_labels)
    #Check break condition for iteration loop
    for phase in phase_labels:
        if phased_reads[phase] < MIN_EDGE_COV:
            phase_below_cov = True
    prev_parts.add(tuple(part_list))
    #Prepare header for iteration stats
    #Iter,Rep Lens,Confirmed/Rejected Pos,Partitioned Reads
    header_labels = ["Iter"]
    for phase in phase_labels:
        header_labels.extend(["Len {0}".format(phase)])
    header_labels.extend(["Confirmed Pos", "Rejected Pos"])
    for phase in phase_labels:
        header_labels.extend(["#Reads {0}".format(phase)])
    header_labels.extend(["#Tied", "#Unassigned"])
    spaced_header = map("{:11}".format, header_labels)
    #Write stats output
    with open(stats_file, 'w') as f:
        f.write("{0:25}\t{1}\n".format("Edge:", edge))
        f.write("{0:25}\t'{1}'\n".format("Side:", side))
        f.write("{0:25}\t".format("Phases:"))
        f.write(", ".join(map(str,phase_labels)))
        f.write("\n")
        f.write("{0:25}\t{1}\n\n".format("Template Length:", template_len))
        f.write("Initial Option Values\n")
        f.write("{0:25}\t{1}\n".format("min_overlap:", min_overlap))
        f.write("{0:25}\t{1}\n".format("sub_thresh:", SUB_THRESH))
        f.write("{0:25}\t{1}\n".format("del_thresh:", DEL_THRESH))
        f.write("{0:25}\t{1}\n".format("ins_thresh:", INS_THRESH))
        f.write("{0:25}\t{1}\n".format("flanking_len:", FLANKING_LEN))
        f.write("{0:25}\t{1}\n".format("buffer_count:", BUFFER_COUNT))
        f.write("{0:25}\t{1}\n".format("max_iter:", MAX_ITER))
        f.write("{0:25}\t{1}\n".format("min_edge_cov:", MIN_EDGE_COV))
        f.write("{0:25}\t{1}\n".format("cons_aln_rate:", CONS_ALN_RATE))
        f.write("\n")
        f.write("The following numbers are calculated based on moving ")
        f.write("'{0}' from the highly divergent window\n\n".format(side))
        f.write("{0}\n".format("Divergent Positions:"))
        f.write("{0:25}\t{1}\n".format("Total", len(pos["total"])))
        f.write("{0:25}\t{1}\n".format("Substitutions", len(pos["sub"])))
        f.write("{0:25}\t{1}\n".format("Deletions", len(pos["del"])))
        f.write("{0:25}\t{1}\n".format("Insertions", len(pos["ins"])))
        f.write("\n")
        f.write("{0:25}\t{1}\n".format("Total Starting Reads:", 
                                       sum(phased_reads.values())))
        for phase in phase_labels:
            f.write("{0}{1}{2:18}\t{3}\n".format("Phase ", phase, 
                                                 " starting reads:", 
                                                 phased_reads[phase]))
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")
    
    return phase_below_cov


def update_side_stats(phases, it, side, cons_align_path, template, 
                      confirmed_pos_path, partitioning, prev_parts, 
                      stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MIN_EDGE_COV = trestle_config.vals["min_edge_cov"]
    #Write stats for each iteration
    #Iter,Rep Lens,Confirmed/Rejected Pos,Partitioned Reads
    stats_out = [str(it)]
    for phase_id in phases:
        hap_len = 0
        if os.path.isfile(cons_align_path.format(it, side, phase_id)):
            cons_align = _read_alignment(cons_align_path.format(it, side, 
                                                                phase_id), 
                                         template, 
                                         CONS_ALN_RATE)
            if cons_align and cons_align[0]:
                if side == "left":
                    hap_len = cons_align[0][0].qry_len
                elif side == "right":
                    hap_len = cons_align[0][0].qry_len
        stats_out.extend([str(hap_len)])
    confirmed_total = 0
    rejected_total = 0
    if it > 0:
        confirmed, rejected, pos = _read_confirmed_positions(
                                            confirmed_pos_path)
        confirmed_total = len(confirmed["total"])
        rejected_total = len(rejected["total"])
    stats_out.extend([str(confirmed_total), 
                      str(rejected_total)])
    phase_below_cov = False
    dup_part = False
    part_list = _read_partitioning_file(partitioning)
    phased_reads, tied_reads, unassigned_reads = _get_partitioning_info(
                                                        part_list, phases)
    for phase_id in sorted(phases):
        stats_out.extend([str(phased_reads[phase_id])])
    stats_out.extend([str(tied_reads), str(unassigned_reads)])
    #Check break conditions for iteration loop
    for phase in phases:
        if phased_reads[phase] < MIN_EDGE_COV:
            phase_below_cov = True
    if tuple(part_list) in prev_parts:
        dup_part = True
    else:
        prev_parts.add(tuple(part_list))
    spaced_header = map("{:11}".format, stats_out)
    with open(stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")
        
    return phase_below_cov, dup_part


def finalize_side_stats(phases, it, side, cons_align_path, template, 
                 cons_vs_cons_path, consensuses, confirmed_pos_path, 
                 partitioning, phase_below_cov, dup_part, term_bool, 
                 stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MAX_ITER = trestle_config.vals["max_iter"]

    with open(stats_file, "a") as f:
        f.write("\n\n")
        f.write("{0:26}\t{1}\n\n".format("Final Iter:", it))
        f.write("Iteration terminated because:\n")
        if it == MAX_ITER:
            f.write("Max iter reached\n")
        if phase_below_cov:
            f.write("Phase coverage fell below min_edge_cov\n")
        if dup_part:
            f.write("Partitioning was identical to a previous iteration\n")
        if term_bool:
            f.write("Encountered empty consensus sequence or alignment\n")
        f.write("\n")
        #Write out alignment indices for phases vs template
        limit_ind = None
        limit_label = ""
        if side == "right":
            limit_label = "Min Template End"
        elif side == "left":
            limit_label = "Max Template Start"
        for phase_id in sorted(phases):
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            curr_cons_path = cons_align_path.format(it, side, phase_id)
            if os.path.isfile(curr_cons_path):
                cons_align = _read_alignment(curr_cons_path, 
                                             template, 
                                             CONS_ALN_RATE)
                if cons_align and cons_align[0]:
                    #collapse multiple consensus alignments
                    coll_cons = _collapse_cons_aln(cons_align)
                    qry_start = coll_cons.qry_start
                    qry_end = coll_cons.qry_end
                    qry_len = coll_cons.qry_len
                    trg_start = coll_cons.trg_start
                    trg_end = coll_cons.trg_end
                    trg_len = coll_cons.trg_len
                    if limit_ind is None or (
                            (side == "right" and trg_end < limit_ind) or
                            (side == "left" and trg_start >= limit_ind)):
                        if side == "right":
                            limit_ind = trg_end
                        elif side == "left":
                            limit_ind = trg_start
            f.write("Phase {0}|Template Alignment\n".format(phase_id))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Phase ", phase_id, ":", 
                    qry_start, qry_end, qry_len))
            f.write("{0:26}\t{1:5}-{2:5} of {3:5}\n".format("Template:",  
                    trg_start, trg_end, trg_len))
        f.write("\n")
        f.write("{0:26}\t{1}\n".format(limit_label, limit_ind))
        f.write("(End of positions considered)\n\n")
        #Write out alignment indices for phases vs phases
        phase_pairs = sorted(combinations(phases, 2))
        for phase_one, phase_two in phase_pairs:
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            qry_seq = ""
            trg_seq = ""
            if (os.path.isfile(cons_vs_cons_path.format(it, side, phase_one, 
                                                       it, side, phase_two)) and
                os.path.isfile(consensuses[(it, side, phase_two)])):
                cons_vs_cons = _read_alignment(cons_vs_cons_path.format(
                                                    it, side, phase_one, 
                                                    it, side, phase_two), 
                                               consensuses[(it, side, 
                                                            phase_two)], 
                                               CONS_ALN_RATE)
                if cons_vs_cons and cons_vs_cons[0]:
                    qry_start = cons_vs_cons[0][0].qry_start
                    qry_end = cons_vs_cons[0][0].qry_end
                    qry_len = cons_vs_cons[0][0].qry_len
                    trg_start = cons_vs_cons[0][0].trg_start
                    trg_end = cons_vs_cons[0][0].trg_end
                    trg_len = cons_vs_cons[0][0].trg_len
                    qry_seq = cons_vs_cons[0][0].qry_seq
                    trg_seq = cons_vs_cons[0][0].trg_seq
            f.write("Phase {0}|Phase {1} Alignment\n".format(phase_one, phase_two))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Phase ", phase_one, ":", 
                    qry_start, qry_end, qry_len))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Phase ", phase_two, ":", 
                    trg_start, trg_end, trg_len))
            div_rate = _calculate_divergence(qry_seq, trg_seq)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
        #Write overall position stats
        types = ["total", "sub", "del", "ins"]
        confirmed = {t:[] for t in types}
        rejected = {t:[] for t in types}
        pos = {t:[] for t in types}
        if it > 0:
            confirmed_pos_output = _read_confirmed_positions(confirmed_pos_path)
            confirmed, rejected, pos = confirmed_pos_output
        if side == "right":
            largest_pos = -1
            if confirmed["total"]:
                largest_pos = max(confirmed["total"])
            f.write("{0:26}\t{1}\n".format("Largest Confirmed Position:", 
                                           largest_pos))
        elif side == "left":
            smallest_pos = -1
            if confirmed["total"]:
                smallest_pos = min(confirmed["total"])
            f.write("{0:26}\t{1}\n".format("Smallest Confirmed Position:", 
                                           smallest_pos))
        remainings = {}
        for typ in types:
            remainings[typ] = len(pos[typ]) - (len(confirmed[typ]) + 
                                               len(rejected[typ]))
        type_strs = ["Total", "Sub", "Del", "Ins"]
        for typ, typ_str in zip(types, type_strs):
            confirmed_frac = 0.0
            rejected_frac = 0.0
            remaining_frac = 0.0
            if len(pos[typ]) != 0:
                confirmed_frac = len(confirmed[typ]) / float(len(pos[typ]))
                rejected_frac = len(rejected[typ]) / float(len(pos[typ]))
                remaining_frac = remainings[typ] / float(len(pos[typ]))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Confirmed {0} Positions:".format(typ_str), 
                            len(confirmed[typ]), 
                            len(pos[typ]), 
                            confirmed_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Rejected {0} Positions:".format(typ_str), 
                            len(rejected[typ]), 
                            len(pos[typ]), 
                            rejected_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Remaining {0} Positions:".format(typ_str), 
                            remainings[typ], 
                            len(pos[typ]), 
                            remaining_frac))
            f.write("\n")
        f.write("\n")
        #Write overall partitioning stats
        part_list = _read_partitioning_file(partitioning)
        phased_reads = {phase:0 for phase in phases}
        tied_reads = 0
        unassigned_reads = 0
        total_reads = len(part_list)
        for part in part_list:
            read_id, status, phase, top_sc, total_sc, header = part
            if status == "Partitioned" and phase != "NA":
                phased_reads[int(phase)] += 1
            elif status == "Tied":
                tied_reads += 1
            elif status == "None":
                unassigned_reads += 1
            else:
                exception_str = "Unknown status {0} in partitioning file {1}"
                raise Exception(exception_str.format(status, partitioning))
        for phase_id in sorted(phases):
            phase_frac = 0.0
            if total_reads != 0:
                phase_frac = phased_reads[phase_id]/float(total_reads)
            f.write("{0}{1}{2:13}\t{3}/{4} = {5:.4f}\n".format(
                                    "Total Phase ", phase_id, " Reads:", 
                                    phased_reads[phase_id], total_reads, 
                                    phase_frac))
        tied_frac = 0.0
        unassigned_frac = 0.0
        if total_reads != 0:
            tied_frac = tied_reads/float(total_reads)
            unassigned_frac = unassigned_reads/float(total_reads)
        f.write("{0:26}\t{1}/{2} = {3:.4f}\n".format("Total Tied Reads:", 
                                              tied_reads, total_reads, 
                                              tied_frac))
        f.write("{0:26}\t{1}/{2} = {3:.4f}\n".format("Total Unassigned Reads:", 
                                      unassigned_reads, total_reads, 
                                      unassigned_frac))
        f.write("\n")


def init_int_stats(edge, side_labels, phase_labels, zero_it, position_path, 
                   partitioning, all_reads_file, template_len, cov, 
                   int_stats_file):
    #Count phased reads
    side_reads = {}
    side_remaining_reads = {}
    for side in sorted(side_labels):
        part_list = _read_partitioning_file(partitioning.format(zero_it, side))
        total_side_reads = len(part_list)
        partitioning_outputs = _get_partitioning_info(part_list, 
                                                      phase_labels)
        side_reads[side], tied_reads, unassigned_reads = partitioning_outputs
        all_side_reads = sum(side_reads[side].values())
        side_remaining_reads[side] = total_side_reads - all_side_reads
    all_reads_n50, num_reads = _n50(all_reads_file)
    #Prepare header for iterative integrated stats
    #in/out Iter,Mean in/out/gap Len,Confirmed/Rejected Pos,Bridging Reads
    header_labels = []
    for side in sorted(side_labels):
        header_labels.extend(["{0} Iter".format(side)])
    header_labels.extend(["left Len", "right Len"])
    header_labels.extend(["Confirmed", "Rejected"])
    spaced_header = map("{:8}".format, header_labels)
    #Write to file
    with open(int_stats_file, 'w') as f:
        f.write("{0:16}\t{1}\n".format("Edge:", edge))
        f.write("{0:16}\t{1}\n".format("Template Length:", template_len))
        f.write("{0:16}\t{1:.2f}\n".format("Avg Coverage:", cov))
        f.write("{0:16}\t{1}\n".format("# All Reads:", num_reads))
        f.write("{0:16}\t{1}\n\n".format("All Reads N50:", all_reads_n50))
        phase_headers = ["Side", "Phase", "# Reads"]
        spaced_phase_header = map("{:5}".format, phase_headers)
        f.write("\t".join(spaced_phase_header))
        f.write("\n")
        for side in sorted(side_labels):
            for phase_id in phase_labels:
                phase_values = [side, phase_id, side_reads[side][phase_id]]
                spaced_values = map("{:6}".format, phase_values)
                f.write("\t".join(spaced_values))
                f.write("\n")
            rem_values = [side, "remain", side_remaining_reads[side]]
            spaced_rem = map("{:6}".format, rem_values)
            f.write("\t".join(spaced_rem))
            f.write("\n")
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")
        

def update_int_stats(edge, side_labels, phase_labels, side_it, cons_align_path, 
                     template, template_len, confirmed_pos_path,
                     int_confirmed_path, partitioning, int_stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    
    stats_out = []
    #Add side iters
    for side in sorted(side_labels):
        stats_out.extend([str(side_it[side])])
    #Find median in, out, and gap lengths
    medians = {s:0 for s in side_labels}
    for side in sorted(side_labels):
        trg_limits = []
        for phase_id in phase_labels:
            curr_cons_path = cons_align_path.format(side_it[side], 
                                                    side, phase_id)
            if os.path.isfile(curr_cons_path):
                cons_align = _read_alignment(curr_cons_path, 
                                             template, 
                                             CONS_ALN_RATE)
                if cons_align and cons_align[0]:
                    if side == "right":
                        trg_limits.append(cons_align[0][0].qry_len)
                    elif side == "left":
                        trg_limits.append(cons_align[0][0].qry_len)
        if trg_limits:
            medians[side] = _get_median(trg_limits)
    stats_out.extend([str(medians["left"]), str(medians["right"])])
    #Add confirmed and rejected reads
    left_confirmed_path = confirmed_pos_path.format(side_it["left"], "left")
    right_confirmed_path = confirmed_pos_path.format(side_it["right"], "right")
    types = ["total", "sub", "del", "ins"]
    int_confirmed = {t:[] for t in types}
    int_rejected = {t:[] for t in types}
    pos = {t:[] for t in types}
    if side_it["left"] > 0 and side_it["right"] > 0:
        all_left_pos = _read_confirmed_positions(left_confirmed_path)
        all_right_pos = _read_confirmed_positions(right_confirmed_path)
        confirmed_pos_outputs = _integrate_confirmed_pos(all_left_pos,
                                                         all_right_pos)
        int_confirmed, int_rejected, pos = confirmed_pos_outputs
    elif side_it["left"] > 0:
        all_left_pos = _read_confirmed_positions(left_confirmed_path)
        int_confirmed, int_rejected, pos = all_left_pos
    elif side_it["right"] > 0:
        all_right_pos = _read_confirmed_positions(right_confirmed_path)
        int_confirmed, int_rejected, pos = all_right_pos
    _write_confirmed_positions(int_confirmed, int_rejected, pos, 
                               int_confirmed_path.format(side_it["left"], 
                                                         side_it["right"]))
    stats_out.extend([str(len(int_confirmed["total"])), 
                      str(len(int_rejected["total"]))])
    spaced_header = map("{:8}".format, stats_out)
    #Write to file
    with open(int_stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")
        

def finalize_int_stats(edge, side_labels, phase_labels, side_it, 
                       cons_align_path, template, template_len, 
                       cons_vs_cons_path, consensuses, int_confirmed_path, 
                       partitioning, int_stats_file, phased_seq_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    PHASED_THRESH = trestle_config.vals["phased_thresh"]
    
    #Phased edge seqs to be returned, NOT written
    phased_edges = {}
    summ_vals = []
    with open(int_stats_file, "a") as f:
        f.write("\n\n")
        for side in sorted(side_labels):
            f.write("{0}'{1}'{2:8}\t{3}\n".format("Final ", side, " Iter:", 
                                              side_it[side]))
        f.write("\n\n")
        #Overall confirmed and rejected positions
        types = ["total", "sub", "del", "ins"]
        int_confirmed = {t:[] for t in types}
        int_rejected = {t:[] for t in types}
        pos = {t:[] for t in types}
        if side_it["left"] > 0 or side_it["right"] > 0:
            int_confirmed, int_rejected, pos = _read_confirmed_positions(
                int_confirmed_path.format(side_it["left"], side_it["right"]))
        remainings = {}
        for typ in types:
            remainings[typ] = len(pos[typ]) - (len(int_confirmed[typ]) + 
                                               len(int_rejected[typ]))
        type_strs = ["Total", "Sub", "Del", "Ins"]
        for typ, typ_str in zip(types, type_strs):
            confirmed_frac = 0.0
            rejected_frac = 0.0
            remaining_frac = 0.0
            if len(pos[typ]) != 0:
                confirmed_frac = len(int_confirmed[typ]) / float(len(pos[typ]))
                rejected_frac = len(int_rejected[typ]) / float(len(pos[typ]))
                remaining_frac = remainings[typ] / float(len(pos[typ]))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Confirmed {0} Positions:".format(typ_str), 
                            len(int_confirmed[typ]), 
                            len(pos[typ]), 
                            confirmed_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Rejected {0} Positions:".format(typ_str), 
                            len(int_rejected[typ]), 
                            len(pos[typ]), 
                            rejected_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Remaining {0} Positions:".format(typ_str), 
                            remainings[typ], 
                            len(pos[typ]), 
                            remaining_frac))
            f.write("\n")
        f.write("\n")
        #Basic stats for confirmed positions
        av_div = 0.0
        if template_len != 0:
            av_div = len(int_confirmed["total"]) / float(template_len)
        position_gaps = [0 for x in range(len(int_confirmed["total"]) + 1)]
        curr_pos = 0
        for i, p in enumerate(int_confirmed["total"]):
            position_gaps[i] = p - curr_pos
            curr_pos = p
        position_gaps[-1] = template_len - curr_pos
        mean_position_gap = _mean(position_gaps)
        max_position_gap = max(position_gaps)
        f.write("{0:26}\t{1}\n".format("Template Length:", template_len))
        f.write("{0:26}\t{1}\n".format("# Confirmed Positions:", 
                                       len(int_confirmed["total"])))
        f.write("{0:26}\t{1:.4f}\n".format("Confirmed Pos Avg Divergence:", 
                                              av_div))
        f.write("{0:26}\t{1:.2f}\n".format("Mean Confirmed Pos Gap:", 
                                              mean_position_gap))
        f.write("{0:26}\t{1}\n".format("Max Confirmed Pos Gap:", 
                                              max_position_gap))
        f.write("\n\n")
        summ_vals.extend([len(int_confirmed["total"]), max_position_gap])
        
        #Write mean in and out divergence rates
        div_rates = {s:[] for s in side_labels}
        for side in sorted(side_labels):
            phase_pairs = sorted(combinations(phase_labels, 2))
            for phase_one, phase_two in phase_pairs:
                cons_cons_file = cons_vs_cons_path.format(
                                    side_it[side], side, phase_one, 
                                    side_it[side], side, phase_two)
                if (os.path.isfile(cons_cons_file) and 
                    os.path.isfile(consensuses[(side_it[side], 
                                                side, phase_two)])):
                    cons_vs_cons = _read_alignment(
                            cons_cons_file, 
                            consensuses[(side_it[side], side, phase_two)], 
                            CONS_ALN_RATE)
                    if cons_vs_cons and cons_vs_cons[0]:
                        phase_rate = _calculate_divergence(
                                        cons_vs_cons[0][0].qry_seq, 
                                        cons_vs_cons[0][0].trg_seq)
                        div_rates[side].append(phase_rate)
        mean_left_div = 0.0
        if div_rates["left"]:
            mean_left_div = _mean(div_rates["left"])
        mean_right_div = 0.0
        if div_rates["right"]:
            mean_right_div = _mean(div_rates["right"])
        f.write("{0:30}\t{1}\n".format("Mean left Divergence Rate:", 
                                        mean_left_div))
        f.write("{0:30}\t{1}\n\n".format("Mean right Divergence Rate:", 
                                        mean_right_div))
        
        
        phased = True
        #Get total length of phased sequence and write lengths
        phase_side_lens = {s:[] for s in side_labels}
        phase_rem_gap = {p:template_len for p in phase_labels}
        phased_headers = []
        for phase_id in phase_labels:
            left_align = None
            right_align = None
            for side in sorted(side_labels):
                curr_cons_path = cons_align_path.format(side_it[side],
                                                         side, phase_id)
                if os.path.isfile(curr_cons_path):
                    cons_align = _read_alignment(
                            curr_cons_path, 
                            template, 
                            CONS_ALN_RATE)
                    if cons_align and cons_align[0]:
                        #collapse multiple consensus alignments
                        coll_cons_align = _collapse_cons_aln(cons_align)
                        if side == "left":
                            left_align = coll_cons_align
                        elif side == "right":
                            right_align = coll_cons_align
            f.write("Phase {0}\n".format(phase_id))
            if not left_align:
                left_start = 0
                left_end = 0
                left_trg_end = 0
                f.write("\tLeft:\n\t\tNo alignment\n")
            else:
                left_start = left_align.qry_start
                left_end = left_align.qry_end
                left_trg_end = left_align.trg_end
                phase_side_lens["left"].append(left_align.qry_len)
                phase_rem_gap[phase_id] += left_align.trg_start
                f.write("\tLeft:\n\t\tqry {0} - {1} of {2}\n".format(
                    left_align.qry_start, left_align.qry_end, left_align.qry_len))
                f.write("\t\ttrg {0} - {1} of {2}\n".format(
                    left_align.trg_start, left_align.trg_end, left_align.trg_len))
                
            if not right_align:
                right_start = 0
                right_end = 0
                right_trg_start = 0
                right_qry_seq = ""
                right_trg_seq = ""
                f.write("\tRight:\n\t\tNo alignment\n")
            else:
                right_start = right_align.qry_start
                right_end = right_align.qry_end
                right_qry_seq = right_align.qry_seq
                right_trg_seq = right_align.trg_seq
                right_trg_start = right_align.trg_start
                phase_side_lens["right"].append(right_align.qry_len)
                phase_rem_gap[phase_id] -= right_align.trg_end    
                f.write("\tRight:\n\t\tqry {0} - {1} of {2}\n".format(
                    right_align.qry_start, right_align.qry_end, right_align.qry_len))
                f.write("\t\ttrg {0} - {1} of {2}\n".format(
                    right_align.trg_start, right_align.trg_end, right_align.trg_len))
            f.write("\tRemaining unspanned:\n\t\t{0}\n\n".format(
                phase_rem_gap[phase_id]))
            if phase_rem_gap[phase_id] / float(template_len) > PHASED_THRESH:
                phased = False
            
            #in sequence used to represent overlapping segment
            adj_right_trg_start = left_trg_end
            adj_right_qry_start = None
            trg_end_ind = None
            right_aln_ind = None
            right_qry_aln, right_aln_qry = _index_mapping(right_qry_seq)
            right_trg_aln, right_aln_trg = _index_mapping(right_trg_seq)
            if adj_right_trg_start < right_trg_start:
                adj_right_qry_start = right_start
            else:
                trg_end_ind = adj_right_trg_start - right_trg_start
                if 0 <= trg_end_ind < len(right_trg_aln):
                    right_aln_ind = right_trg_aln[trg_end_ind]
                    if 0 <= right_aln_ind < len(right_aln_qry):
                        adj_right_qry_start = (right_start + 
                                                right_aln_qry[right_aln_ind])
                elif trg_end_ind == len(right_trg_aln):
                    adj_right_qry_start = right_end
            if adj_right_qry_start:
                right_start = adj_right_qry_start            
            f.write("\tAdj Right trg start\t{0}\n".format(adj_right_trg_start))
            f.write("\tAdj trg ind        \t{0}\n".format(trg_end_ind))
            f.write("\tAdj aln ind        \t{0}\n".format(right_aln_ind))
            f.write("\tNew Right qry start\t{0}\n\n".format(right_start))
            
            print "Here C", side_it["left"], side_it["right"]
            header = "edge_{0}_haplotype_{1}".format(edge, phase_id)
            copy_seq = ""
            if side_it["left"] > 0 and side_it["right"] > 0:
                copy_seq = _construct_repeat_copy(
                        consensuses[(side_it["left"], "left", phase_id)], 
                        consensuses[(side_it["right"], "right", phase_id)], 
                        left_start, left_end, 
                        right_start, right_end)
            print "Here D"
            phased_edges[header] = copy_seq
            if copy_seq:
                seq_dict = {header:copy_seq}
                fp.write_fasta_dict(seq_dict, 
                                    phased_seq_file.format(edge, phase_id))
            phased_headers.append(header)
        summ_vals.extend([phased, ":".join(phased_headers)])
        
        #Calculate and write weighted mean div rate
        mean_left_len = 0.0
        if phase_side_lens["left"]:
            mean_left_len = _mean(phase_side_lens["left"])
        mean_right_len = 0.0
        if phase_side_lens["right"]:
            mean_right_len = _mean(phase_side_lens["right"])
        weighted_mean_div = 0.0
        if mean_left_len + mean_right_len != 0:
            weighted_mean_div = (mean_left_div*mean_left_len + 
                                 mean_right_div*mean_right_len) / float(
                                 mean_left_len + mean_right_len)
        f.write("{0:30}\t{1}\n\n".format("Weighted Mean Divergence Rate:", 
                                      weighted_mean_div))
    return phased, phased_edges, summ_vals


def int_stats_postscript(edge, side_labels, phase_labels, integrated_stats, 
                         phased_edge_path, res_vs_res):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    
    divs = []
    with open(integrated_stats, "a") as f:
        res_inds = phase_labels
        f.write("Phased Edge Sequence Alignments\n")
        for res_one, res_two in sorted(combinations(res_inds, 2)):
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            qry_seq = ""
            trg_seq = ""
            if os.path.isfile(res_vs_res.format(edge, res_one, res_two) and
                phased_edge_path.format(edge, res_two)):
                res_align = _read_alignment(res_vs_res.format(edge, res_one, 
                                                              res_two), 
                                            phased_edge_path.format(edge, 
                                                                     res_two), 
                                            CONS_ALN_RATE)
                if res_align and res_align[0]:
                    qry_start = res_align[0][0].qry_start
                    qry_end = res_align[0][0].qry_end
                    qry_len = res_align[0][0].qry_len
                    trg_start = res_align[0][0].trg_start
                    trg_end = res_align[0][0].trg_end
                    trg_len = res_align[0][0].trg_len
                    qry_seq = res_align[0][0].qry_seq
                    trg_seq = res_align[0][0].trg_seq
            f.write("Phase {0}|Phase {1}\n".format(res_one, res_two))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Phase ", res_one, ":", 
                    qry_start, qry_end, qry_len))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Phase ", res_two, ":",   
                    trg_start, trg_end, trg_len))
            div_rate = _calculate_divergence(qry_seq, trg_seq)
            divs.append(div_rate)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
    return _mean(divs)
        

def _get_partitioning_info(part_list, phases):
    phased_reads = {phase:0 for phase in phases}
    tied_reads = 0
    unassigned_reads = 0
    for part in part_list:
        read_id, status, phase, top_sc, total_sc, header = part
        if status == "Partitioned" and phase != "NA":
            phased_reads[int(phase)] += 1
        elif status == "Tied":
            tied_reads += 1
        elif status == "None":
            unassigned_reads += 1
        else:
            exception_str = "Unknown status {0} in partitioning file"
            raise Exception(exception_str.format(status))
    return phased_reads, tied_reads, unassigned_reads


def _calculate_divergence(qry_seq, trg_seq):
    if not qry_seq or not trg_seq:
        return 0.0
        
    curr_del = 0
    curr_ins = 0
        
    match_count = 0
    mis_count = 0
    del_count = 0
    ins_count = 0
    
    for q, t in zip(qry_seq, trg_seq):
        if q == t:
            match_count += 1
            if curr_del != 0:
                del_count += 1
                curr_del = 0
            if curr_ins != 0:
                ins_count += 1
                curr_ins = 0
        elif q == "-" and t != "-":
            curr_del += 1
            if curr_ins != 0:
                ins_count += 1
                curr_ins = 0
        elif q != "-" and t == "-":
            curr_ins += 1
            if curr_del != 0:
                del_count += 1
                curr_del = 0
        elif q != t:
            mis_count += 1
            if curr_del != 0:
                del_count += 1
                curr_del = 0
            if curr_ins != 0:
                ins_count += 1
                curr_ins = 0
        else:
            raise Exception("No alignment conditions fit, {0} {1}".format(q, 
                                                                          t))
    if curr_del != 0:
        del_count += 1
        curr_del = 0
    if curr_ins != 0:
        ins_count += 1
        curr_ins = 0
    
    indel_sim_rate = 0.0
    total = match_count + mis_count + del_count + ins_count
    if total != 0:
        indel_sim_rate = match_count / float(total)
    return 1 - indel_sim_rate


def _n50(reads_file):
    reads_dict = fp.read_sequence_dict(reads_file)
    read_lengths = sorted([len(x) for x in reads_dict.values()], reverse=True)
    num_reads = len(reads_dict)
    summed_len = 0
    n50 = 0
    for l in read_lengths:
        summed_len += l
        if summed_len >= sum(read_lengths)/2:
            n50 = l
            break
    return n50, num_reads


def _get_median(lst):
    if not lst:
        raise ValueError("_get_median() arg is an empty sequence")
    sorted_list = sorted(lst)
    if len(lst) % 2 == 1:
        return sorted_list[len(lst)/2]
    else:
        mid1 = sorted_list[(len(lst)/2) - 1]
        mid2 = sorted_list[(len(lst)/2)]
        return float(mid1 + mid2) / 2


def _integrate_confirmed_pos(all_left_pos, all_right_pos):
    left_conf, left_rej, left_pos = all_left_pos
    right_conf, right_rej, right_pos = all_right_pos
    
    integrated_confirmed = {"total":[], "sub":[], "ins":[], "del":[]}
    integrated_rejected = {"total":[], "sub":[], "ins":[], "del":[]}
    
    for pos in sorted(left_pos["total"]):
        for pos_type in left_conf:
            if pos in left_conf[pos_type] or pos in right_conf[pos_type]:
                integrated_confirmed[pos_type].append(pos)
            elif pos in left_rej[pos_type] or pos in right_rej[pos_type]:
                integrated_rejected[pos_type].append(pos)
    return integrated_confirmed, integrated_rejected, left_pos


def _get_combos(left_list, right_list):
    all_combos = []
    for combo in _combo_helper(left_list, right_list):
        all_combos.append(combo)
    return all_combos


def _combo_helper(left_list, right_list):
    if not left_list or not right_list:
        yield []
        return
    else:
        left1 = left_list[0]
        for j in range(len(right_list)):
            combo = (left1, right_list[j])
            for rest in _combo_helper(left_list[1:], 
                                      right_list[:j] + right_list[j + 1:]):
                 yield [combo] + rest


def _get_aln_end(aln_start, aln_seq):
    return aln_start+len(aln_seq.replace("-",""))


def _check_overlap(in_file, temp_file, out_file, overlap, in_start, in_end, 
                   temp_start, temp_end, out_start, out_end, new_out_start, 
                   in_qry, in_trg, out_qry, out_trg, out_trg_aln, out_aln_trg, 
                   out_qry_aln, out_aln_qry, out_trg_end, out_qry_end, 
                   in_align, out_align):
    in_dict = fp.read_sequence_dict(in_file)
    in_seq = in_dict.values()[0]
    temp_dict = fp.read_sequence_dict(temp_file)
    temp_seq = temp_dict.values()[0]
    out_dict = fp.read_sequence_dict(out_file)
    out_seq = out_dict.values()[0]  
    for i in range(len(out_qry)/50-1, len(out_qry)/50+1):
        aln_ind_st = i*50
        aln_ind_end = (i+1)*50
        if aln_ind_end > len(out_qry):
            aln_ind_end = len(out_qry)
        print 'ALN inds', aln_ind_st, aln_ind_end
        qry_ind_st = out_aln_qry[aln_ind_st]
        if aln_ind_end < len(out_aln_qry):
            qry_ind_end = out_aln_qry[aln_ind_end]
        else:
            qry_ind_end = out_aln_qry[-1]
        print 'QRY inds', qry_ind_st, qry_ind_end
        trg_ind_st = out_aln_trg[aln_ind_st]
        if aln_ind_end < len(out_aln_trg):
            trg_ind_end = out_aln_trg[aln_ind_end]
        else:
            trg_ind_end = out_aln_trg[-1]
        print 'TRG inds', trg_ind_st, trg_ind_end  
        
        print "TRG ALN", out_trg_aln[trg_ind_st:trg_ind_end]
        print "ALN TRG", out_aln_trg[aln_ind_st:aln_ind_end]
        print "QRY ALN", out_qry_aln[qry_ind_st:qry_ind_end]
        print "ALN QRY", out_aln_qry[aln_ind_st:aln_ind_end]
        print "QRY SEQ", out_qry[aln_ind_st:aln_ind_end]
        print "TRG SEQ", out_trg[aln_ind_st:aln_ind_end]
        print
    print 'In end, in template end',in_end,temp_start
    print 'AR In qry end',in_qry[-10:]
    print 'AR In trg end',in_trg[-10:]
    print 'Out old start, old end, new start, out template start', out_start, 
    print out_end, new_out_start, temp_end
    print "Out_trg_end", out_trg_end
    print "Out_qry_end", out_qry_end
    print "In align qry inds", in_align.qry_start, in_align.qry_end, 
    print in_align.qry_len
    print "In align trg inds", in_align.trg_start, in_align.trg_end, 
    print in_align.trg_len
    print "Out align qry inds", out_align.qry_start, out_align.qry_end, 
    print out_align.qry_len
    print "Out align trg inds", out_align.trg_start, out_align.trg_end, 
    print out_align.trg_len
    print
    print "Overlap:\t{0}".format(overlap)
    print "In seq(-30 to end):\t{0}".format(in_seq[in_end-30:in_end])
    temp_end_seq = temp_seq[temp_start-30:temp_start]
    print "Template seq(-30 to end):\t{0}".format(temp_end_seq)
    #print "Out seq:\t{0}".format(out_seq[out_start:out_end])
    #print "AR In seq:\t{0}".format(in_seq[in_start-10:in_end+10])
    #print "AR Template seq:\t{0}".format(temp_seq[temp_end:temp_start+10])
    #print "AR Out seq:\t{0}".format(out_seq[out_start:out_end+10])
    pre_new_out = out_seq[new_out_start-30:new_out_start]
    post_new_out = out_seq[new_out_start:new_out_start+30]
    print "New out seq(-30 to new start):\t{0}".format(pre_new_out)
    print "New out seq(new_start to +30):\t{0}".format(post_new_out)
    print

def _construct_repeat_copy(left_file, right_file, left_start, left_end, 
                           right_start, right_end):
    if (not os.path.isfile(left_file) or 
        not os.path.isfile(right_file)):
        return ""
    left_dict = fp.read_sequence_dict(left_file)
    left_seq = left_dict.values()[0]
    right_dict = fp.read_sequence_dict(right_file)
    right_seq = right_dict.values()[0]
    seq = ""
    if (0 <= left_start <= len(left_seq) and
            0 <= left_end <= len(left_seq) and
            0 <= right_start <= len(right_seq) and
            0 <= right_end <= len(right_seq)):
        seq = ''.join([left_seq[left_start:left_end], 
                       right_seq[right_start:right_end]])
    return seq

def init_summary(summary_file):
    with open(summary_file, "w") as f:
        summ_header_labels = ["Edge_ID", "Template", "Cov", 
                              "#Conf_Pos", "Max_Pos_Gap", "Phased?", "Avg_Div", 
                              "Seq_Headers"]
        #spaced_header = map("{:13}".format, summ_header_labels)
        f.write(" ".join(map(lambda x: "{:<13}".format(str(x)),
                             summ_header_labels)))
        f.write("\n")


def update_summary(summ_items, summary_file):
    (edge_id, template_len, avg_cov, summ_vals, 
     avg_div, both_phased_present) = summ_items
    
    (confirmed_pos, max_pos_gap, phased, seq_headers) = tuple(summ_vals)

    avg_cov = "{:.1f}".format(avg_cov)
    avg_div = "{:.4f}".format(avg_div)
    phased_present = phased and both_phased_present

    summ_out = [edge_id, template_len, avg_cov, confirmed_pos, 
                max_pos_gap, phased_present, avg_div, seq_headers]
    with open(summary_file, "a") as f:
        f.write(" ".join(map(lambda x: "{:<13}".format(str(x)),
                             summ_out)))
        f.write("\n")

def remove_unneeded_files(repeat_edges, rep, side_labels, side_it, orient_dir, 
                          template, extended, pol_temp_dir, pol_ext_dir, 
                          pre_edge_reads, pre_partitioning, pre_read_align, 
                          partitioning, cons_align, cut_cons_align, 
                          read_align, confirmed_pos_path, edge_reads, 
                          cut_cons, polishing_dir, cons_vs_cons, 
                          int_confirmed_path, repeat_reads, frequency_path, 
                          alignment_file, num_pol_iters, iter_pairs):
    add_dir_name = "additional_output"
    add_dir = os.path.join(orient_dir, add_dir_name)
    if not os.path.isdir(add_dir):
        os.mkdir(add_dir)
    pol_name = "polished_{0}.fasta".format(num_pol_iters)
    pol_template = "polished_template.fasta"
    pol_ext = "polished_extended.{0}.{1}.fasta"
    pol_temp_file = os.path.join(pol_temp_dir, pol_name)
    if os.path.exists(pol_temp_file):
        os.rename(pol_temp_file, os.path.join(add_dir, pol_template))
    for side in side_labels:
        for edge_id in repeat_edges[rep][side]:
            pol_ext_file = os.path.join(pol_ext_dir.format(side, edge_id), 
                                        pol_name)
            if os.path.exists(pol_ext_file):
                os.rename(pol_ext_file, 
                          os.path.join(add_dir,
                                       pol_ext.format(side, edge_id)))
    
    files_to_remove = [template]
    dirs_to_remove = [pol_temp_dir]
    files_to_move = [repeat_reads, frequency_path, alignment_file]
    if os.path.exists(pol_temp_dir):
        for fil in os.listdir(pol_temp_dir):
            files_to_remove.append(os.path.join(pol_temp_dir, fil))
    
    for side in side_labels:
        for edge_id in repeat_edges[rep][side]:
            files_to_remove.append(extended.format(side, edge_id))
            curr_pol_ext_dir = pol_ext_dir.format(side, edge_id)
            dirs_to_remove.append(curr_pol_ext_dir)
            if os.path.exists(curr_pol_ext_dir):
                for fil in os.listdir(curr_pol_ext_dir):
                    files_to_remove.append(os.path.join(curr_pol_ext_dir, fil))
            files_to_remove.append(pre_edge_reads.format(side, edge_id))
            files_to_remove.append(pre_read_align.format(side, edge_id))
            for it in range(1, side_it[side] + 1):
                files_to_remove.append(cons_align.format(it, side, edge_id))
                files_to_remove.append(read_align.format(it, side, edge_id))
                files_to_remove.append(edge_reads.format(it, side, edge_id))
                pol_cons = polishing_dir.format(it, side, edge_id)
                dirs_to_remove.append(pol_cons)
                if os.path.exists(pol_cons):
                    for fil in os.listdir(pol_cons):
                        files_to_remove.append(os.path.join(pol_cons, fil))
            for it in range(1, side_it[side]):
                files_to_remove.append(cut_cons_align.format(it, side, edge_id))
                files_to_remove.append(cut_cons.format(it, side, edge_id))
            it = side_it[side]
            files_to_move.append(cut_cons_align.format(it, side, edge_id))
            files_to_move.append(cut_cons.format(it, side, edge_id))
        
        edge_pairs = sorted(combinations(repeat_edges[rep][side], 2))
        for edge_one, edge_two in edge_pairs:
            for it in range(1, side_it[side]):
                cons_cons_file = cons_vs_cons.format(it, side, edge_one, 
                                                     it, side, edge_two)
                files_to_remove.append(cons_cons_file)
            it = side_it[side]
            cons_cons_file = cons_vs_cons.format(it, side, edge_one, 
                                                 it, side, edge_two)
            files_to_move.append(cons_cons_file)
        files_to_remove.append(pre_partitioning.format(side))
        for it in range(1, side_it[side]):
            files_to_remove.append(partitioning.format(it, side))
            files_to_remove.append(confirmed_pos_path.format(it, side))
        for it in [0, side_it[side]]:
            files_to_move.append(partitioning.format(it, side))
        it = side_it[side]
        files_to_move.append(confirmed_pos_path.format(it, side))                     
    
    last_conf_pos = int_confirmed_path.format(side_it[side_labels[0]], 
                                             side_it[side_labels[1]])
    for it1, it2 in iter_pairs:
        curr_conf_pos = int_confirmed_path.format(it1, it2)
        if curr_conf_pos != last_conf_pos:
            files_to_remove.append(curr_conf_pos)
        else:
            files_to_move.append(curr_conf_pos)
    
    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)
    for d in dirs_to_remove:
        if os.path.exists(d):
            os.rmdir(d)
    for f in files_to_move:
        if os.path.exists(f):
            split_path = os.path.split(f)
            new_file = os.path.join(split_path[0], add_dir_name, split_path[1])
            os.rename(f, new_file)


def _mean(lst):
    if not lst:
        return 0
    return float(sum(lst)) / len(lst)
