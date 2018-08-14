# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 03:50:31 2017

@author: jeffrey_yuan
"""

import os
import logging
import numpy as np
from itertools import combinations, product

import flye.alignment as flye_aln
import flye.fasta_parser as fp
import flye.config as config
import flye.bubbles as bbl
import flye.polish as pol

import flye.trestle.divergence as div
import flye.trestle.trestle_config as trestle_config

logger = logging.getLogger()

def resolve_repeats(args, trestle_dir, repeats_dump, graph_edges, summ_file):
    SUB_THRESH = trestle_config.vals["sub_thresh"]
    DEL_THRESH = trestle_config.vals["del_thresh"]
    INS_THRESH = trestle_config.vals["ins_thresh"]
    MAX_ITER = trestle_config.vals["max_iter"]
    
    #Defining directory and file names for trestle output
    repeat_label = "repeat_{0}"
    orient_labels = ["forward", "reverse"]
    side_labels = ["in", "out"]
    pol_temp_name = "Polishing.Template"
    pol_ext_name = "Polishing.Extended.{0}.{1}"
    pol_cons_dir_name = "Polishing.Consensus.{0}.{1}.{2}"
    template_name = "template.fasta"
    extended_name = "extended_templates.{0}.{1}.fasta"
    repeat_reads_name = "repeat_reads.fasta"
    partitioning_name = "partitioning.{0}.{1}.txt"
    div_freq_name = "divergence_frequencies.txt"
    div_pos_name = "divergent_positions.txt"
    div_summ_name = "divergence_summary.txt"
    pol_cons_name = "polished_consensus.{0}.{1}.fasta"
    overext_aln_name = "overext.{0}.{1}.{2}.vs.template.minimap.sam"
    reads_template_aln_name = "reads.vs.template.minimap.sam"
    cons_template_aln_name = "consensus.{0}.{1}.{2}.vs.template.minimap.sam"
    reads_cons_aln_name = "reads.vs.consensus.{0}.{1}.{2}.minimap.sam"
    confirmed_pos_name = "confirmed_positions.{0}.{1}.txt"
    edge_reads_name = "edge_reads.{0}.{1}.{2}.fasta"
    cons_vs_cons_name = "".join(["consensus.{0}.{1}.{2}.vs.",
                                 "consensus.{3}.{4}.{5}.minimap.sam"])
    side_stats_name = "{0}_stats.txt"
    int_stats_name = "integrated_stats.txt"
    int_confirmed_pos_name = "integrated_confirmed_positions.{0}.{1}.txt"
    resolved_rep_name = "resolved_repeat_{0}.copy.{1}.fasta"
    res_vs_res_name = "resolved_repeat_{0}.copy.{1}.vs.{2}.minimap.sam"
    """"""
    test_pos_name = "test_pos.{0}.{1}.txt"
    
    zero_it = 0
    all_resolved_reps_dict = {}
    init_summary(summ_file)
    
    #1. Process repeats from graph - generates a folder for each repeat
    logger.debug("Finding unbridged repeats")
    process_outputs = process_repeats(args.reads, repeats_dump, graph_edges, 
                                      trestle_dir, repeat_label, orient_labels, 
                                      template_name, extended_name, 
                                      repeat_reads_name, partitioning_name, 
                                      side_labels, zero_it)
    repeat_list, repeat_edges, all_edge_headers = process_outputs
    logger.info("Repeats to be resolved: {0}".format(len(repeat_list)))
    for rep_id in sorted(repeat_list):
        logger.info("Resolving repeat '{0}'".format(rep_id))
        repeat_dir = os.path.join(trestle_dir, 
                                  repeat_label.format(rep_id))
        orient_reps = [rep_id, -rep_id]
        
        for orientation, rep in zip(orient_labels, orient_reps):
            orient_dir = os.path.join(repeat_dir, orientation)
            template = os.path.join(orient_dir, template_name)
            extended = os.path.join(orient_dir, extended_name)
            repeat_reads = os.path.join(orient_dir, repeat_reads_name)
            term_bool = {s:False for s in side_labels}
            
            #2. Polish template and extended templates
            logger.debug("Polishing templates")
            pol_temp_dir = os.path.join(orient_dir, pol_temp_name)
            polished_template = _run_polishing(args, [repeat_reads], template, 
                                               pol_temp_dir)

            if not os.path.isfile(polished_template):
                for side in side_labels:
                    term_bool[side] = True
            
            polished_extended = {}
            for side in side_labels:
                for edge_id in repeat_edges[rep][side]:
                    pol_ext_dir = os.path.join(orient_dir, pol_ext_name)
                    pol_output = _run_polishing(args, [repeat_reads], 
                                                extended.format(side, edge_id), 
                                                pol_ext_dir.format(side, edge_id))                    
                    polished_extended[(side, edge_id)] = pol_output
                    if not os.path.isfile(pol_output):
                        term_bool[side] = True
            
            #3. Find divergent positions
            logger.debug("Estimating divergence")
            frequency_path = os.path.join(orient_dir, div_freq_name)
            position_path = os.path.join(orient_dir, div_pos_name)
            summary_path = os.path.join(orient_dir, div_summ_name)
            
            #logger.info("running Minimap2")
            alignment_file = os.path.join(orient_dir, reads_template_aln_name)
            template_len = 0.0
            if os.path.isfile(polished_template):
                flye_aln.make_alignment(polished_template, [repeat_reads], 
                           args.threads, orient_dir, args.platform, 
                           alignment_file)
                template_info = flye_aln.get_contigs_info(polished_template)
                template_len = template_info[str(rep)].length
            
            logger.debug("Finding tentative divergence positions")
            div.find_divergence(alignment_file, polished_template, 
                                template_info, frequency_path, position_path, 
                                summary_path, config.vals["min_aln_rate"], 
                                args.platform, args.threads,
                                SUB_THRESH, DEL_THRESH, INS_THRESH) 
            avg_cov = find_coverage(frequency_path)
            
            #4. Initialize paths, variables, and stats
            partitioning = os.path.join(orient_dir, partitioning_name)
            cons_align = os.path.join(orient_dir, cons_template_aln_name)
            read_align = os.path.join(orient_dir, reads_cons_aln_name)
            confirmed_pos_path = os.path.join(orient_dir, confirmed_pos_name)
            edge_reads = os.path.join(orient_dir, edge_reads_name)
            polishing_dir = os.path.join(orient_dir, pol_cons_dir_name)
            polished_cons = os.path.join(orient_dir, pol_cons_name)
            overext_aln_path = os.path.join(orient_dir, overext_aln_name)
            cons_vs_cons = os.path.join(orient_dir, cons_vs_cons_name)
            side_stats = os.path.join(orient_dir, side_stats_name)
            integrated_stats = os.path.join(orient_dir, int_stats_name)
            int_confirmed_path = os.path.join(orient_dir, 
                                              int_confirmed_pos_name)
            resolved_rep_path = os.path.join(orient_dir, resolved_rep_name)
            res_vs_res = os.path.join(orient_dir, res_vs_res_name)
            """""" #define test variables
            test_pos = os.path.join(orient_dir, test_pos_name)
            num_test = 10
            
            overext_consensus = {}
            side_it = {s:0 for s in side_labels}
            edge_below_cov = {s:False for s in side_labels}
            dup_part = {s:False for s in side_labels}
            prev_partitionings = {s:set() for s in side_labels}
            #Initialize stats
            for side in side_labels:
                edge_below_cov[side] = init_side_stats(
                                    rep, side, repeat_edges, args.min_overlap, 
                                    position_path, 
                                    partitioning.format(zero_it, side), 
                                    prev_partitionings[side], 
                                    template_len, 
                                    side_stats.format(side))
            init_int_stats(rep, repeat_edges, zero_it, position_path, 
                           partitioning, repeat_reads, template_len, 
                           avg_cov, integrated_stats)
            #5. Start iterations
            logger.debug("Iterative procedure")
            for it in range(1, MAX_ITER + 1):
                both_break = True
                for side in side_labels:
                    if (edge_below_cov[side] or dup_part[side] or 
                        term_bool[side]):
                        continue
                    else:
                        logger.debug("iteration {0}, '{1}'".format(it, side))
                        both_break = False
                    #5a. Call consensus on partitioned reads
                    for edge_id in sorted(repeat_edges[rep][side]):
                        pol_con_dir = polishing_dir.format(
                                    it, side, edge_id)
                        curr_reads = edge_reads.format(it, side, edge_id)
                        write_edge_reads(
                                    it, side, edge_id,
                                    repeat_reads, 
                                    partitioning.format(it-1, side), 
                                    curr_reads)
                        curr_extended = polished_extended[(side, edge_id)]
                        pol_reads_str = "polishing '{0} {1}' reads"
                        logger.debug(pol_reads_str.format(side, edge_id))
                        pol_con_out = _run_polishing(args, [curr_reads], 
                                                curr_extended, 
                                                pol_con_dir)
                        overext_consensus[(it, side, edge_id)] = pol_con_out
                    #5a+. Add to consensus sequences
                    for edge_id in sorted(repeat_edges[rep][side]):
                        curr_overext = overext_consensus[(it, side, edge_id)]
                        if os.path.isfile(curr_overext):
                            overext_aln = overext_aln_path.format(it, side, 
                                                                  edge_id)
                            flye_aln.make_alignment(polished_template, 
                                                    [curr_overext],
                                                    args.threads, orient_dir,
                                                    args.platform, overext_aln)
                            _add_to_consensus(polished_cons.format(side, 
                                                                   edge_id), 
                                              curr_overext, polished_template, 
                                              overext_aln, it, side, edge_id)
                        else:
                            term_bool[side] = True
                    #5b. Partition reads using divergent positions
                    for edge_id in sorted(repeat_edges[rep][side]):
                        curr_cons = polished_cons.format(side, edge_id)
                        if os.path.isfile(curr_cons):
                            cons_al_file = cons_align.format(it, side, edge_id)
                            flye_aln.make_alignment(polished_template, 
                                               [curr_cons], args.threads, 
                                                orient_dir, args.platform, 
                                                cons_al_file)
                            read_al_file = read_align.format(it, side, edge_id)
                            flye_aln.make_alignment(curr_cons, [repeat_reads], 
                                               args.threads, orient_dir, 
                                               args.platform, read_al_file)
                        else:
                            term_bool[side] = True
                    logger.debug("partitioning '{0}' reads".format(side))
                    partition_reads(repeat_edges[rep][side], it, side, 
                                       position_path, cons_align, 
                                       polished_template, read_align, 
                                       polished_cons,
                                       confirmed_pos_path, 
                                       partitioning, 
                                       all_edge_headers[rep], 
                                       test_pos.format(it, side), num_test)
                    #5c. Write stats file for current iteration
                    edge_pairs = sorted(combinations(repeat_edges[rep][side], 
                                                     2))
                    for edge_one, edge_two in edge_pairs:
                        cons_one = polished_cons.format(side, edge_one)
                        cons_two = polished_cons.format(side, edge_two)
                        if (not os.path.isfile(cons_one) or 
                            not os.path.isfile(cons_two)):
                            continue
                        cons_cons_file = cons_vs_cons.format(
                                                it, side, edge_one, 
                                                it, side, edge_two)
                        flye_aln.make_alignment(cons_two, 
                                           [cons_one], 
                                            args.threads, 
                                            orient_dir, 
                                            args.platform, 
                                            cons_cons_file)
                    side_stat_outputs = update_side_stats(
                                        repeat_edges[rep][side], it, side, 
                                        cons_align, polished_template, 
                                        confirmed_pos_path.format(it, side), 
                                        partitioning.format(it, side), 
                                        prev_partitionings[side], 
                                        side_stats.format(side))
                    edge_below_cov[side], dup_part[side] = side_stat_outputs
                    side_it[side] = it
                update_int_stats(rep, repeat_edges, side_it, cons_align, 
                                    polished_template, 
                                    template_len,
                                    confirmed_pos_path, int_confirmed_path, 
                                    partitioning, integrated_stats)
                if both_break:
                    break
            #6. Finalize stats files
            logger.debug("Writing stats files")
            for side in side_labels:
                finalize_side_stats(repeat_edges[rep][side], side_it[side], 
                                    side, cons_align, polished_template, 
                                    cons_vs_cons, polished_cons, 
                                    confirmed_pos_path.format(side_it[side], 
                                                              side), 
                                    partitioning.format(side_it[side], side), 
                                    edge_below_cov[side],
                                    dup_part[side], term_bool[side], 
                                    side_stats.format(side))
            final_int_outputs = finalize_int_stats(rep, repeat_edges, side_it, 
                                                   cons_align, 
                                                   polished_template, 
                                                   template_len, cons_vs_cons, 
                                                   polished_cons, 
                                                   int_confirmed_path, 
                                                   partitioning, 
                                                   integrated_stats, 
                                                   resolved_rep_path)
            bridged, repeat_seqs, summ_vals = final_int_outputs
            #7. Generate summary and resolved repeat file
            logger.debug("Generating summary and resolved repeat file")
            avg_div = 0.0
            if bridged:
                res_inds = range(len(repeat_edges[rep]["in"]))
                for res_one, res_two in sorted(combinations(res_inds, 2)):
                    res_one_path = resolved_rep_path.format(rep, res_one)
                    res_two_path = resolved_rep_path.format(rep, res_two)
                    if (os.path.isfile(res_one_path) and
                        os.path.isfile(res_two_path)):
                        flye_aln.make_alignment(res_two_path, 
                                           [res_one_path], 
                                           args.threads, 
                                           orient_dir, 
                                           args.platform, 
                                           res_vs_res.format(rep, res_one, 
                                                             res_two))
                avg_div = int_stats_postscript(rep, repeat_edges, 
                                               integrated_stats, 
                                               resolved_rep_path, res_vs_res)
            all_resolved_reps_dict.update(repeat_seqs)
            update_summary(rep, template_len, avg_cov, summ_vals, avg_div, 
                           summ_file)
    return all_resolved_reps_dict

#Process Repeats functions
class ProcessingException(Exception):
    pass
       
def process_repeats(reads, repeats_dump, graph_edges, work_dir, repeat_label, 
                    orient_labels, template_name, extended_name, 
                    repeat_reads_name, partitioning_name, side_labels, 
                    zero_it):
    """Generates repeat dirs and files given reads, repeats_dump and
    graph_edges files. Only returns repeats between min_mult and max_mult"""
    MIN_MULT = trestle_config.vals["min_mult"]
    MAX_MULT = trestle_config.vals["max_mult"]
    FLANKING_LEN = trestle_config.vals["flanking_len"]
    
    #Reads input files
    repeats_dict = _read_repeats_dump(repeats_dump)
    if not repeats_dict:
        logger.debug("Empty repeats_dump file: {0}".format(repeats_dump))
        return [], {}, {}
        
    reads_dict = {}
    for read_file in reads:
        reads_dict.update(fp.read_fast_file(read_file))
    orig_graph = fp.read_fasta_dict(graph_edges)
    graph_dict = {int(h.split('_')[1]):orig_graph[h] for h in orig_graph}
    
    if not reads_dict:
        raise ProcessingException("No reads found from {0}".format(reads))
    if not graph_dict:
        raise ProcessingException("No edges found from {0}".format(
            graph_edges))
        
    repeat_list = []
    repeat_edges = {}
    all_edge_headers = {}
    for rep in sorted(repeats_dict, reverse=True):
        #Checks multiplicity of repeat and presence of reverse strand
        #One run processes both forward and reverse strand of repeat
                
        if rep <= 0:
            continue
        
        valid_repeat = True
        if -rep not in repeats_dict:
            logger.debug("Repeat {0} missing reverse strand".format(rep))
            valid_repeat = False
        if (repeats_dict[rep][0] < MIN_MULT or
                 repeats_dict[rep][0] > MAX_MULT or
                 repeats_dict[-rep][0] < MIN_MULT or
                 repeats_dict[-rep][0] > MAX_MULT):
            logger.debug("Repeat {0} multiplicity not in range: {1}".format(
                                                 rep, repeats_dict[rep][0]))
            valid_repeat = False
        if rep not in graph_dict:
            logger.debug("Repeat {0} missing from graph file".format(rep))
            valid_repeat = False
        if not valid_repeat:
            continue
        
        #Makes repeat dirs
        repeat_dir = os.path.join(work_dir, repeat_label.format(rep))
        if not os.path.isdir(repeat_dir):
            os.mkdir(repeat_dir)
        repeat_list.append(rep)
                
        orient_reps = [rep, -rep]
        for curr_label, curr_rep in zip(orient_labels, orient_reps):
            orient_path = os.path.join(repeat_dir, curr_label)
            if not os.path.isdir(orient_path):
                os.mkdir(orient_path)
            template_path = os.path.join(orient_path, template_name)
            extended_path = os.path.join(orient_path, extended_name)
            repeat_reads_path = os.path.join(orient_path, repeat_reads_name)
            partitioning_path = os.path.join(orient_path, partitioning_name)
            
            in_label = side_labels[0]
            out_label = side_labels[1]
            repeat_edges[curr_rep] = {in_label:[], out_label:[]}
            
            repeat_parts = repeats_dict[curr_rep]
            mult, all_reads_list, inputs_dict, outputs_dict = repeat_parts
            
            template_dict = {}
            extended_dicts = {}
            repeat_reads_dict = {}
            #Partitioning parts: id_num, Partitioned/Tied/None, 
            #edge_id, top_score, total_score, Header
            partitioning = {in_label:[], out_label:[]}
            read_id = 0
            
            template_seq = graph_dict[rep]
            if curr_label == "reverse":
                template_seq = fp.reverse_complement(graph_dict[rep])
            template_dict[curr_rep] = template_seq
            
            all_edge_headers[curr_rep] = {}
            out_headers = set()
            #Headers will be in the form -h or +h,
            #edge_dict is in the form >[Input,Output]_edge##_h,
            #rev_comp of read will be written if the header is -h
            for edge_id in inputs_dict:         
                repeat_edges[curr_rep][in_label].append(edge_id)
                extended_dicts[(in_label, edge_id)] = {}
                
                headers = inputs_dict[edge_id]
                for header in headers:
                    if (not header) or (header[0] != '+' and header[0] != '-'):
                        raise ProcessingException(
                            "Input read format not recognized: {0}".format(
                                header))
                    if header[1:] not in reads_dict:
                        raise ProcessingException(
                            "Read header {0} not in any of {1}".format(
                                header[1:], reads))
                    
                    if header[1:] not in all_edge_headers[curr_rep]:
                        status_label = "Partitioned"
                        edge_label = str(edge_id)
                        score = 1
                        total_score = 0
                        partitioning[in_label].append((read_id, status_label, 
                                                       edge_label, score, 
                                                       total_score, 
                                                       header[1:]))
                        all_edge_headers[curr_rep][header[1:]] = read_id
                        read_id += 1
                
                extend_in_header = "Extended_Template_Input_{0}".format(
                    edge_id)
                if edge_id > 0:
                    edge_seq = graph_dict[edge_id]
                elif edge_id < 0:
                    edge_seq = fp.reverse_complement(graph_dict[-edge_id])
                extended_seq = edge_seq[-FLANKING_LEN:]
                extended_dicts[(in_label, edge_id)][extend_in_header] = (
                                        extended_seq + template_seq)
                
            
            for edge_id in outputs_dict:         
                repeat_edges[curr_rep][out_label].append(edge_id)
                extended_dicts[(out_label, edge_id)] = {}
                
                headers = outputs_dict[edge_id]
                for header in headers:
                    if (not header) or (header[0] != '+' and header[0] != '-'):
                        raise ProcessingException(
                        "Output read format not recognized: {0}".format(
                            header))
                    if header[1:] not in reads_dict:
                        raise ProcessingException(
                        "Read header {0} not in any of {1}".format(
                            header[1:], reads))
                    
                    curr_read_id = read_id
                    if header[1:] not in all_edge_headers[curr_rep]:
                        status_label = "None"
                        edge_label = "NA"
                        score = 0
                        total_score = 0
                        partitioning[in_label].append((read_id, status_label, 
                                                       edge_label, score, 
                                                       total_score,
                                                       header[1:]))
                        
                        all_edge_headers[curr_rep][header[1:]] = read_id
                        read_id += 1
                    else:
                        curr_read_id = all_edge_headers[curr_rep][header[1:]]
                    
                    if header[1:] not in out_headers:
                        status_label = "Partitioned"
                        edge_label = str(edge_id)
                        score = 1
                        total_score = 0
                        partitioning[out_label].append((curr_read_id, 
                                                        status_label, 
                                                        edge_label, score, 
                                                        total_score, 
                                                        header[1:]))
                        out_headers.add(header[1:])   
                
                extend_out_header = "Extended_Template_Output_{0}".format(
                    edge_id)
                if edge_id > 0:
                    edge_seq = graph_dict[edge_id]
                elif edge_id < 0:
                    edge_seq = fp.reverse_complement(graph_dict[-edge_id])
                extended_seq = edge_seq[:FLANKING_LEN]
                extended_dicts[(out_label, edge_id)][extend_out_header] = (
                                        template_seq + extended_seq)
            
            #Need to reiterate over in_headers to add in_headers to 
            #out-partitioning while avoiding double-adding ones in both
            for edge_id in inputs_dict:
                headers = inputs_dict[edge_id]
                for header in headers:
                    if header[1:] not in out_headers:
                        curr_read_id = all_edge_headers[curr_rep][header[1:]]
                        status_label = "None"
                        edge_label = "NA"
                        score = 0
                        total_score = 0
                        partitioning[out_label].append((curr_read_id, 
                                                        status_label, 
                                                        edge_label, score, 
                                                        total_score, 
                                                        header[1:]))
                
            
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
                repeat_reads_dict[header[1:]] = seq
                
                curr_read_id = read_id
                if header[1:] not in all_edge_headers[curr_rep]:
                    all_edge_headers[curr_rep][header[1:]] = read_id
                    read_id += 1
                    
                    status_label = "None"
                    edge_label = "NA"
                    score = 0
                    total_score = 0
                    partitioning[in_label].append((curr_read_id, status_label, 
                                                   edge_label, score, 
                                                   total_score, header[1:]))
                                                   
                    status_label = "None"
                    edge_label = "NA"
                    score = 0
                    total_score = 0
                    partitioning[out_label].append((curr_read_id, status_label, 
                                                   edge_label, score, 
                                                   total_score, header[1:]))
            
            if template_dict and template_dict.values()[0]:
                fp.write_fasta_dict(template_dict, template_path)
            for edge in extended_dicts:
                if extended_dicts[edge] and extended_dicts[edge].values()[0]:
                    extended_edge_path = extended_path.format(edge[0], 
                                                              edge[1])
                    fp.write_fasta_dict(extended_dicts[edge], 
                                        extended_edge_path)
            if repeat_reads_dict and repeat_reads_dict.values()[0]:
                fp.write_fasta_dict(repeat_reads_dict, repeat_reads_path)
            for side in side_labels:
                _write_partitioning_file(partitioning[side], 
                                         partitioning_path.format(zero_it, 
                                                                  side))
            
            if not template_dict:
                raise ProcessingException("No template {0} found".format(
                                                curr_rep))
            for edge in extended_dicts:
                if not template_dict:
                    raise ProcessingException(
                        "No extended template {0} {1} {2} found".format(
                            curr_rep, edge[0], edge[1]))
            if not repeat_reads_dict:
                raise ProcessingException("No repeat reads {0} found".format(
                                                curr_rep))
            for side in side_labels:
                if not partitioning[side]:
                    raise ProcessingException(
                        "Empty partitioning file {0}".format(
                            partitioning_path.format(zero_it, side)))
    return repeat_list, repeat_edges, all_edge_headers

def _read_repeats_dump(repeats_dump):
    """Read repeats_dump.txt file following format guidelines"""
    repeats_dict = {}
    
    curr_repeat = 0
    mult = 0
    all_reads_list = []
    all_bool = False
    inputs_dict = {}
    curr_input = 0
    curr_input_list = []
    input_bool = False
    outputs_dict = {}
    curr_output = 0
    curr_output_list = []
    output_bool = False
    with open(repeats_dump,'r') as rf:
        for i,line in enumerate(rf):
            line = line.strip()
            if line:
                if line[0] == '#':
                    all_bool = False
                    if input_bool:
                        inputs_dict[curr_input] = curr_input_list[:]
                        curr_input_list = []
                        curr_input = 0
                        input_bool = False
                    if output_bool:
                        outputs_dict[curr_output] = curr_output_list[:]
                        curr_output_list = []
                        curr_output = 0
                        output_bool = False
                    
                    parts = line.split('\t')
                    header = parts[0]
                    head_parts = header.split(' ')
                    if head_parts[0] == '#Repeat':
                        if curr_repeat != 0:
                            repeats_dict[curr_repeat] = (
                                mult, all_reads_list,
                                inputs_dict, outputs_dict)
                            all_reads_list = []
                            inputs_dict = {}
                            outputs_dict = {}
                        
                        curr_repeat = int(head_parts[1])
                        mult = int(parts[1])
                    elif head_parts[0] == '#All':
                        all_bool = True
                    elif head_parts[0] == '#Input':                        
                        curr_input = int(head_parts[1])
                        input_bool = True
                    elif head_parts[0] == '#Output':
                        curr_output = int(head_parts[1])
                        output_bool = True
                    else:
                        raise ProcessingException(
                        "Unexpected repeats dump file format")
                else:
                    if all_bool:
                        all_reads_list.append(line)
                    elif input_bool:
                        curr_input_list.append(line)
                    elif output_bool:
                        curr_output_list.append(line)
                    else:
                        raise ProcessingException(
                        "Unexpected repeats dump file format")
        if output_bool:
            outputs_dict[curr_output] = curr_output_list[:]
        if input_bool:
            inputs_dict[curr_input] = curr_input_list[:]
        if curr_repeat != 0:
            repeats_dict[curr_repeat] = (
                mult, all_reads_list, inputs_dict, outputs_dict)
    return repeats_dict        
        
def _write_partitioning_file(part_list, part_path):
    with open(part_path, "w") as f:
        header_labels = ["Read_ID", "Status", "Edge", "Top Score", 
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
        coverage = np.mean(all_covs)
        #print min(all_covs), np.mean(all_covs), max(all_covs)
    return coverage

def write_edge_reads(it, side, edge_id, all_reads, partitioning, out_file):
    all_reads_dict = fp.read_fasta_dict(all_reads)
    part_list = _read_partitioning_file(partitioning)
    edge_header_name = "Read_{0}_Iter_{1}_{2}_{3}_edge_{4}"
    edge_reads = {}
    for part in part_list:
        read_id, status, edge, top_sc, total_sc, header = part
        if status == "Partitioned" and edge != "NA" and int(edge) == edge_id:
            edge_seq = all_reads_dict[header]
            edge_header = edge_header_name.format(read_id, it, 
                                                  side, edge_id, header)
            edge_reads[edge_header] = edge_seq
    if edge_reads and edge_reads.values()[0]:
        fp.write_fasta_dict(edge_reads, out_file)

#Polishing Functions

def _run_polishing(args, reads, seqs, polish_dir):
    MIN_ALN_RATE = config.vals["min_aln_rate"]
    if not os.path.isdir(polish_dir):
        os.mkdir(polish_dir)
    polished_seqs = os.path.join(polish_dir, 
                                 "polished_{0}.fasta".format(args.num_iters))
    polished_stats = os.path.join(polish_dir, "contigs_stats.txt")
    
    prev_template = seqs
    contig_lengths = None
    for i in xrange(args.num_iters):
        alignment_file = os.path.join(polish_dir,
                                      "minimap_{0}.sam".format(i + 1))
        flye_aln.make_alignment(prev_template, reads, args.threads,
                               polish_dir, args.platform, alignment_file)
        contigs_info = flye_aln.get_contigs_info(prev_template)
        bubbles_file = os.path.join(polish_dir,
                                    "bubbles_{0}.fasta".format(i + 1))
        coverage_stats, err_rate = \
            bbl.make_bubbles(alignment_file, contigs_info, prev_template,
                             args.platform, args.threads, MIN_ALN_RATE, 
                             bubbles_file)
        polished_file = os.path.join(polish_dir,
                                     "polished_{0}.fasta".format(i + 1))
        contig_lengths = pol.polish(bubbles_file, args.threads, args.platform, 
                                    polish_dir, i + 1, polished_file,
                                    output_progress=False)
        prev_template = polished_file
        
    with open(polished_stats, "w") as f:
        f.write("seq_name\tlength\tcoverage\n")
        for ctg_id in contig_lengths:
            f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                    contig_lengths[ctg_id], coverage_stats[ctg_id]))
    
    return polished_seqs

def _add_to_consensus(polished_cons, overext_cons, template, 
                      overext_aln_file, it, side, edge_id):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    EXTEND_LEN = trestle_config.vals["extend_len"]
    polished_dict = {}
    polished_seq = []
    polished_header = "edge_{0}_polished_consensus".format(edge_id)
    if it > 1 and os.path.isfile(polished_cons):
        prev_dict = fp.read_fasta_dict(polished_cons)
        prev_header, prev_seq = prev_dict.items()[0]
        polished_header = prev_header
        polished_seq.append(prev_seq)
    
    skip_bool = False
    overext_aln = []
    if not os.path.isfile(overext_aln_file):
        skip_bool = True
    else:
        overext_aln = _read_alignment(overext_aln_file, template, 
                                      CONS_ALN_RATE)
    if not overext_aln or not overext_aln[0]:
        logger.debut("No overext alignment")
        skip_bool = True
    if not skip_bool:
        if side == "in":
            extend_start = EXTEND_LEN * (it - 1)
            extend_end = EXTEND_LEN * it
        elif side == "out":
            extend_start = overext_aln[0][0].trg_len - (EXTEND_LEN * it)
            extend_end = overext_aln[0][0].trg_len - (EXTEND_LEN * (it - 1))
        trg_aln, aln_trg = _index_mapping(overext_aln[0][0].trg_seq)
        qry_aln, aln_qry = _index_mapping(overext_aln[0][0].qry_seq)
        if side == "in" and extend_end > overext_aln[0][0].trg_end:
            extend_end = overext_aln[0][0].trg_end
        elif side == "out" and extend_start <= overext_aln[0][0].trg_start:
            extend_start = overext_aln[0][0].trg_start
        
        qry_st = 0
        qry_end = overext_aln[0][0].qry_len
        if it == 1 and side == "in":
            ext_end_minus_st = extend_end - overext_aln[0][0].trg_start
            if ext_end_minus_st == len(trg_aln):
                aln_end_ind = trg_aln[-1] + 1
                qry_end = aln_qry[-1] + 1 + overext_aln[0][0].qry_start
            else:
                aln_end_ind = trg_aln[ext_end_minus_st]
                qry_end = aln_qry[aln_end_ind] + overext_aln[0][0].qry_start
        elif it == 1 and side =="out":
            ext_st_minus_st = extend_start - overext_aln[0][0].trg_start
            aln_st_ind = trg_aln[ext_st_minus_st]
            qry_st = aln_qry[aln_st_ind] + overext_aln[0][0].qry_start
        elif ((side == "in" and extend_start >= overext_aln[0][0].trg_start) or 
              (side == "out" and extend_end <= overext_aln[0][0].trg_end)):
            ext_st_minus_st = extend_start - overext_aln[0][0].trg_start
            ext_end_minus_st = extend_end - overext_aln[0][0].trg_start
            aln_st_ind = trg_aln[ext_st_minus_st]
            qry_st = aln_qry[aln_st_ind] + overext_aln[0][0].qry_start
            if ext_end_minus_st == len(trg_aln):
                aln_end_ind = trg_aln[-1] + 1
                qry_end = aln_qry[-1] + 1 + overext_aln[0][0].qry_start
            else:
                aln_end_ind = trg_aln[ext_end_minus_st]
                qry_end = aln_qry[aln_end_ind] + overext_aln[0][0].qry_start
        else:
            skip_bool = True
    if not skip_bool:
        overext_dict = fp.read_fasta_dict(overext_cons)
        overext_seq = overext_dict.values()[0]
        extend_seq = overext_seq[qry_st:qry_end]
        if side == "in":
            polished_seq.append(extend_seq)
        elif side == "out":
            polished_seq.insert(0, extend_seq)
        polished_header += "|it_{0}_template_{1}_{2}".format(
                                            it, extend_start, extend_end)
    polished_dict[polished_header] = "".join(polished_seq)
    fp.write_fasta_dict(polished_dict, polished_cons)

#Partition Reads Functions

def partition_reads(edges, it, side, position_path, cons_align_path, 
                    template, read_align_path, consensuses, 
                    confirmed_pos_path, part_file, 
                    headers_to_id, pos_test, num_test):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    BUFFER_COUNT = trestle_config.vals["buffer_count"]
    
    skip_bool = False
    pos_headers, pos = div.read_positions(position_path)
    """"""
    test_pos = pos["total"][:num_test]
    with open(pos_test, 'w') as f:
        f.write("Test positions ({0})\n".format(num_test))
        f.write("{0}\n\n".format(test_pos))
        
    cons_aligns = {}
    for edge_id in edges:
        if not os.path.isfile(cons_align_path.format(it, side, edge_id)):
            skip_bool = True
        else:
            cons_aligns[edge_id] = _read_alignment(cons_align_path.format(it, 
                                                        side, 
                                                        edge_id), 
                                                   template, 
                                                   CONS_ALN_RATE)
        if (not cons_aligns or 
                not cons_aligns[edge_id] or 
                not cons_aligns[edge_id][0]):
            logger.debug("No cons alignment found for edge {0}".format(
                            edge_id))
            skip_bool = True
    if skip_bool:
        if it == 1:
            confirmed_pos = {"total":[], "sub":[], "ins":[], "del":[]}
            rejected_pos = {"total":[], "sub":[], "ins":[], "del":[]}
            consensus_pos = pos
        else:
            previous_pos = _read_confirmed_positions(
                                confirmed_pos_path.format(it - 1, side))
            confirmed_pos, rejected_pos, consensus_pos = previous_pos
    else:
        curr_pos = _evaluate_positions(pos, cons_aligns, pos_test, 
                                       test_pos, side)
        confirmed_pos, rejected_pos, consensus_pos = curr_pos
    _write_confirmed_positions(confirmed_pos, rejected_pos, pos, 
                               confirmed_pos_path.format(it, side))
    read_aligns = {}
    for edge_id in edges:
        if not os.path.isfile(read_align_path.format(it, side, edge_id)):
            skip_bool = True
        elif not skip_bool:
            read_aligns[edge_id] = _read_alignment(
                                        read_align_path.format(it, side, 
                                                               edge_id), 
                                        consensuses.format(side, edge_id), 
                                        CONS_ALN_RATE)
        if (not read_aligns or 
                not read_aligns[edge_id] or 
                not read_aligns[edge_id][0]):
            read_aln_str = "No read alignment found for edge {0}"
            logger.debug(read_aln_str.format(edge_id))
            skip_bool = True
    if skip_bool:
        partitioning = _read_partitioning_file(part_file.format(it - 1, side))
    else:
        partitioning = _classify_reads(read_aligns, consensus_pos, 
                        headers_to_id, BUFFER_COUNT, pos_test, num_test)
    _write_partitioning_file(partitioning, part_file.format(it, side))
    
def _read_alignment(alignment, target_path, min_aln_rate):
    alignments = []
    aln_reader = flye_aln.SynchronizedSamReader(alignment,
                                       fp.read_fasta_dict(target_path),
                                       min_aln_rate)
    aln_reader.init_reading()
    while not aln_reader.is_eof():
        ctg_id, ctg_aln = aln_reader.get_chunk()
        if ctg_id is None:
            break
        alignments.append(ctg_aln)
    
    return alignments
    
    
class EdgeAlignment:
    __slots__ = ("edge_id", "qry_seq", "trg_seq", "qry_start", "trg_start", 
                 "trg_end", "in_alignment", "curr_aln_ind", "curr_qry_nuc", 
                 "curr_trg_nuc", "curr_ins_nuc")

    def __init__(self, edge_id, qry_seq, trg_seq, qry_start, trg_start, 
                 trg_end):
        self.edge_id = edge_id
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

def _evaluate_positions(pos, cons_aligns, pos_test, test_pos, side):
    #Includes insertions!
    confirmed_pos = {"total":[], "sub":[], "ins":[], "del":[]}
    rejected_pos = {"total":[], "sub":[], "ins":[], "del":[]}
    consensus_pos = {e:[] for e in cons_aligns}    
    
    alns = {}
    for edge_id in cons_aligns:
        orig_aln = cons_aligns[edge_id][0][0]
        alns[edge_id] = EdgeAlignment(edge_id, orig_aln.qry_seq, 
                                      orig_aln.trg_seq, orig_aln.qry_start, 
                                      orig_aln.trg_start, orig_aln.trg_end)
    
    min_start_edge = min([alns[e].trg_start for e in alns])
    max_end_edge = max([alns[e].trg_end for e in alns])
    #end indices for conservatively defining confirmed positions
    min_end_edge = min([alns[e].trg_end for e in alns])
    max_start_edge = max([alns[e].trg_start for e in alns])
    
    for trg_ind in range(min_start_edge, max_end_edge):
        for edge_id in alns:
            aln = alns[edge_id]
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
            if ((side == "in" and trg_ind < min_end_edge) or
                (side == "out" and trg_ind >= max_start_edge)):
                """"""
                if trg_ind in test_pos:
                    with open(pos_test, "a") as f:
                        f.write("Pos {0}:\n".format(trg_ind))
                ins_confirmed = False
                del_confirmed = False
                sub_confirmed = False
                qry_nuc = ""
                trg_nuc = ""
                for edge_id in alns:
                    aln = alns[edge_id]
                    """"""
                    if trg_ind in test_pos:
                        with open(pos_test, "a") as f:
                            f.write("Edge {0} in_alignment:\t{1}\n".format(
                                            edge_id, aln.in_alignment))
                    if aln.in_alignment:
                        #Directly add positions only to consensuses 
                        # where insertions occur
                        #Add the position prior to curr_qry_ind to 
                        # account for insertion
                        if aln.curr_ins_nuc:
                            ins_confirmed = True
                            consensus_pos[edge_id].append(aln.curr_qry_ind - 1)
                            """"""
                            if trg_ind in test_pos:
                                with open(pos_test, "a") as f:
                                    ins_str = "ins_confirmed, cons_pos:\t{0}\n"
                                    f.write(ins_str.format(aln.curr_qry_ind-1))
                            
                        if qry_nuc and qry_nuc != aln.curr_qry_nuc:
                            if qry_nuc == "-":
                                del_confirmed = True
                            else:
                                sub_confirmed = True
                        else:
                            qry_nuc = aln.curr_qry_nuc
                        if trg_nuc and trg_nuc != aln.curr_trg_nuc:
                            incon_str = "Inconsistent trg_nuc, {0} {1} {2} {3}"
                            logger.debug(incon_str.format(edge_id, 
                                                          trg_ind, trg_nuc, 
                                                          aln.curr_trg_nuc))
                        trg_nuc = aln.curr_trg_nuc
                        """"""
                        if trg_ind in test_pos:
                            with open(pos_test, "a") as f:
                                f.write("aln_ind:\t{0}\n".format(
                                                        aln.curr_aln_ind))
                                f.write("qry_ind:\t{0}\n".format(
                                                        aln.curr_qry_ind))
                                f.write("qry_nuc:\t{0}\n".format(
                                                        aln.curr_qry_nuc))
                                f.write("trg_nuc:\t{0}\n".format(
                                                        aln.curr_trg_nuc))
                                f.write("ins_nuc:\t{0}\n".format(
                                                        aln.curr_ins_nuc))
                                ar_st = aln.curr_aln_ind - 5
                                ar_end = aln.curr_aln_ind + 5
                                ar_qry = aln.qry_seq[ar_st:ar_end]
                                ar_trg = aln.trg_seq[ar_st:ar_end]
                                f.write("around qry_seq:\t{0}\n".format(
                                                        ar_qry))
                                f.write("around trg_seq:\t{0}\n".format(
                                                        ar_trg))
                                f.write("current trg_nuc:\t{0}\n".format(
                                                        trg_nuc))
                                f.write("current ins_confirmed:\t{0}\n".format(
                                                        ins_confirmed))
                                f.write("current del_confirmed:\t{0}\n".format(
                                                        del_confirmed))
                                f.write("current sub_confirmed:\t{0}\n".format(
                                                        sub_confirmed))
                                f.write("\n")
                if ins_confirmed or del_confirmed or sub_confirmed:
                    confirmed_pos["total"].append(trg_ind)
                    #Add positions to consensuses for only subs/deletions
                    if del_confirmed or sub_confirmed:
                        """"""
                        if trg_ind in test_pos:
                            with open(pos_test, "a") as f:
                                f.write("qry_confirmed\n")
                        for edge_id in alns:
                            aln = alns[edge_id]
                            if aln.in_alignment:
                                consensus_pos[edge_id].append(aln.curr_qry_ind)
                                """"""
                                if trg_ind in test_pos:
                                    with open(pos_test, "a") as f:
                                        edge_str = "Edge {0}, cons_pos:\t{1}\n"
                                        f.write(edge_str.format(edge_id, 
                                                            aln.curr_qry_ind))
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
        
        for edge_id in alns:
            aln = alns[edge_id]
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
                    headers_to_id, buffer_count, pos_test, num_test):
    #Includes insertion positions where an insertion occurs right before the 
    #position for the read.
    #partitioning format same as above:
    #list of (read_id, status, edge_id, top_score, total_score, header)
    partitioning = []
    
    read_scores = {}    
    for edge_id in read_aligns:
        for aln in read_aligns[edge_id][0]:
            read_header = aln.qry_id
            cons_header = aln.trg_id
            #Unmapped segments will not be scored
            if cons_header == "*":
                continue
            if read_header not in read_scores:
                read_scores[read_header] = {}
            read_scores[read_header][edge_id] = 0
            
            positions = consensus_pos[edge_id]
            trg_aln, aln_trg = _index_mapping(aln.trg_seq)
            for pos in positions:
                if pos >= aln.trg_start and pos < aln.trg_end:
                    pos_minus_start = pos - aln.trg_start
                    aln_ind = trg_aln[pos_minus_start]
                    if aln.qry_seq[aln_ind] == aln.trg_seq[aln_ind]:
                        read_scores[read_header][edge_id] += 1
                    """"""
                    if pos in positions[:num_test]:
                        with open(pos_test, "a") as f:
                            nedge_str = "\nEdge {0}, cons_pos {1}, read {2}\n"
                            f.write(nedge_str.format(edge_id, 
                                            pos, headers_to_id[read_header]))
                            f.write("pos_minus_start:\t{0}\n".format(
                                            pos_minus_start))
                            f.write("qry_seq:\t{0}\n".format(
                                            aln.qry_seq[aln_ind]))
                            f.write("trg_seq:\t{0}\n".format(
                                            aln.trg_seq[aln_ind]))
                            f.write("around qry_seq:\t{0}\n".format(
                                            aln.qry_seq[aln_ind-5:aln_ind+5]))
                            f.write("around trg_seq:\t{0}\n".format(
                                            aln.trg_seq[aln_ind-5:aln_ind+5]))
    #Iterate through all read_headers so partitioning will be a complete set
    for read_header in headers_to_id:
        read_id = headers_to_id[read_header]
        if read_header in read_scores:
            tie_bool = False
            top_edge = 0
            top_score = 0
            total_score = 0
            for edge_id in read_scores[read_header]:
                edge_score = read_scores[read_header][edge_id]
                #print edge_id, edge_score, top_score
                if edge_score - buffer_count > top_score:
                    top_edge = edge_id
                    top_score = edge_score
                    tie_bool = False
                elif (edge_score - buffer_count <= top_score and 
                      edge_score >= top_score):
                    top_score = edge_score
                    tie_bool = True
                total_score += edge_score
            
            if total_score == 0:
                status_label = "None"
                edge_label = "NA"
            elif tie_bool:
                status_label = "Tied"
                edge_label = "NA"
            else:
                status_label = "Partitioned"
                edge_label = str(top_edge)
            partitioning.append((read_id, status_label, edge_label, 
                                 top_score, total_score, read_header))
        else:
            status_label = "None"
            edge_label = "NA"
            top_score = 0
            total_score = 0
            partitioning.append((read_id, status_label, edge_label, 
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

def init_side_stats(rep, side, repeat_edges, min_overlap, position_path, 
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
    edge_below_cov = False
    part_list = _read_partitioning_file(partitioning)
    edge_reads, tied_reads, unassigned_reads = _get_partitioning_info(
                                                    part_list, 
                                                    repeat_edges[rep][side])
    #Check break condition for iteration loop
    for edge in repeat_edges[rep][side]:
        if edge_reads[edge] < MIN_EDGE_COV:
            edge_below_cov = True
    prev_parts.add(tuple(part_list))
    #Prepare header for iteration stats
    #Iter,Rep Lens,Confirmed/Rejected Pos,Partitioned Reads
    header_labels = ["Iter"]
    for edge in sorted(repeat_edges[rep][side]):
        header_labels.extend(["Rep Len {0}".format(edge)])
    header_labels.extend(["Confirmed Pos", "Rejected Pos"])
    for edge in sorted(repeat_edges[rep][side]):
        header_labels.extend(["#Reads {0}".format(edge)])
    header_labels.extend(["#Tied", "#Unassigned"])
    spaced_header = map("{:11}".format, header_labels)
    #Write stats output
    with open(stats_file, 'w') as f:
        f.write("{0:25}\t{1}\n".format("Repeat:", rep))
        f.write("{0:25}\t'{1}'\n".format("Side:", side))
        f.write("{0:25}\t".format("Edges:"))
        f.write(", ".join(map(str,sorted(repeat_edges[rep][side]))))
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
        f.write("into the repeat from the '{0}' direction\n\n".format(side))
        f.write("{0}\n".format("Divergent Positions:"))
        f.write("{0:25}\t{1}\n".format("Total", len(pos["total"])))
        f.write("{0:25}\t{1}\n".format("Substitutions", len(pos["sub"])))
        f.write("{0:25}\t{1}\n".format("Deletions", len(pos["del"])))
        f.write("{0:25}\t{1}\n".format("Insertions", len(pos["ins"])))
        f.write("\n")
        f.write("{0:25}\t{1}\n".format("Total Starting Reads:", 
                                       sum(edge_reads.values())))
        for edge in sorted(repeat_edges[rep][side]):
            f.write("{0}{1}{2:18}\t{3}\n".format("Edge ", edge, 
                                                 " starting reads:", 
                                                 edge_reads[edge]))
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")
    
    return edge_below_cov

def update_side_stats(edges, it, side, cons_align_path, template, 
                      confirmed_pos_path, partitioning, prev_parts, 
                      stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MIN_EDGE_COV = trestle_config.vals["min_edge_cov"]
    #Write stats for each iteration
    #Iter,Rep Lens,Confirmed/Rejected Pos,Partitioned Reads
    stats_out = [str(it)]
    for edge_id in sorted(edges):
        rep_len = 0
        if os.path.isfile(cons_align_path.format(it, side, edge_id)):
            cons_align = _read_alignment(cons_align_path.format(it, side, 
                                                                edge_id), 
                                         template, 
                                         CONS_ALN_RATE)
            if cons_align and cons_align[0]:
                if side == "in":
                    rep_len = (cons_align[0][0].qry_len - 
                                cons_align[0][0].qry_start)
                elif side == "out":
                    rep_len = cons_align[0][0].qry_end
        stats_out.extend([str(rep_len)])
    confirmed, rejected, pos = _read_confirmed_positions(confirmed_pos_path)
    stats_out.extend([str(len(confirmed["total"])), 
                      str(len(rejected["total"]))])
    edge_below_cov = False
    dup_part = False
    part_list = _read_partitioning_file(partitioning)
    edge_reads, tied_reads, unassigned_reads = _get_partitioning_info(
                                                        part_list, edges)
    for edge_id in sorted(edges):
        stats_out.extend([str(edge_reads[edge_id])])
    stats_out.extend([str(tied_reads), str(unassigned_reads)])
    #Check break conditions for iteration loop
    for edge in edges:
        if edge_reads[edge] < MIN_EDGE_COV:
            edge_below_cov = True
    if tuple(part_list) in prev_parts:
        dup_part = True
    else:
        prev_parts.add(tuple(part_list))
    spaced_header = map("{:11}".format, stats_out)
    with open(stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")
        
    return edge_below_cov, dup_part

def finalize_side_stats(edges, it, side, cons_align_path, template, 
                 cons_vs_cons_path, consensuses, confirmed_pos_path, 
                 partitioning, edge_below_cov, dup_part, term_bool, 
                 stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MAX_ITER = trestle_config.vals["max_iter"]

    with open(stats_file, "a") as f:
        f.write("\n\n")
        f.write("{0:26}\t{1}\n\n".format("Final Iter:", it))
        f.write("Iteration terminated because:\n")
        if it == MAX_ITER:
            f.write("Max iter reached\n")
        if edge_below_cov:
            f.write("Edge coverage fell below min_edge_cov\n")
        if dup_part:
            f.write("Partitioning was identical to a previous iteration\n")
        if term_bool:
            f.write("Encountered empty consensus sequence or alignment\n")
        f.write("\n")
        #Write out alignment indices for edges vs template
        limit_ind = None
        limit_label = ""
        if side == "in":
            limit_label = "Min Template End"
        elif side == "out":
            limit_label = "Max Template Start"
        for edge_id in sorted(edges):
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            curr_cons_path = cons_align_path.format(it, side, edge_id)
            if os.path.isfile(curr_cons_path):
                cons_align = _read_alignment(curr_cons_path, 
                                             template, 
                                             CONS_ALN_RATE)
            if cons_align and cons_align[0]:
                qry_start = cons_align[0][0].qry_start
                qry_end = cons_align[0][0].qry_end
                qry_len = cons_align[0][0].qry_len
                trg_start = cons_align[0][0].trg_start
                trg_end = cons_align[0][0].trg_end
                trg_len = cons_align[0][0].trg_len
                if limit_ind is None or (
                        (side == "in" and trg_end < limit_ind) or
                        (side == "out" and trg_start >= limit_ind)):
                    if side == "in":
                        limit_ind = trg_end
                    elif side == "out":
                        limit_ind = trg_start
            f.write("Edge {0}|Template Alignment\n".format(edge_id))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Edge ", edge_id, ":", 
                    qry_start, qry_end, qry_len))
            f.write("{0:26}\t{1:5}-{2:5} of {3:5}\n".format("Template:",  
                    trg_start, trg_end, trg_len))
        f.write("\n")
        f.write("{0:26}\t{1}\n".format(limit_label, limit_ind))
        f.write("(End of positions considered)\n\n")
        #Write out alignment indices for edges vs edges
        edge_pairs = sorted(combinations(edges, 2))
        for edge_one, edge_two in edge_pairs:
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            qry_seq = ""
            trg_seq = ""
            if os.path.isfile(cons_vs_cons_path.format(it, side, edge_one, 
                                                       it, side, edge_two)):
                cons_vs_cons = _read_alignment(cons_vs_cons_path.format(
                                                    it, side, edge_one, 
                                                    it, side, edge_two), 
                                               consensuses.format(side, 
                                                                  edge_two), 
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
            f.write("Edge {0}|Edge {1} Alignment\n".format(edge_one, edge_two))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Edge ", edge_one, ":", 
                    qry_start, qry_end, qry_len))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Edge ", edge_two, ":", 
                    trg_start, trg_end, trg_len))
            div_rate = _calculate_divergence(qry_seq, trg_seq)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
        #Write overall position stats
        confirmed_pos_output = _read_confirmed_positions(confirmed_pos_path)
        confirmed, rejected, pos = confirmed_pos_output
        if side == "in":
            largest_pos = -1
            if confirmed["total"]:
                largest_pos = max(confirmed["total"])
            f.write("{0:26}\t{1}\n".format("Largest Confirmed Position:", 
                                           largest_pos))
        elif side == "out":
            smallest_pos = -1
            if confirmed["total"]:
                smallest_pos = min(confirmed["total"])
            f.write("{0:26}\t{1}\n".format("Smallest Confirmed Position:", 
                                           smallest_pos))
        types = ["total", "sub", "del", "ins"]
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
        edge_reads = {edge:0 for edge in edges}
        tied_reads = 0
        unassigned_reads = 0
        total_reads = len(part_list)
        for part in part_list:
            read_id, status, edge, top_sc, total_sc, header = part
            if status == "Partitioned" and edge != "NA":
                edge_reads[int(edge)] += 1
            elif status == "Tied":
                tied_reads += 1
            elif status == "None":
                unassigned_reads += 1
            else:
                exception_str = "Unknown status {0} in partitioning file {1}"
                raise Exception(exception_str.format(status, partitioning))
        for edge_id in sorted(edges):
            f.write("{0}{1}{2:13}\t{3}/{4} = {5:.4f}\n".format(
                                    "Total Edge ", edge_id, " Reads:", 
                                    edge_reads[edge_id], total_reads, 
                                    edge_reads[edge_id]/float(total_reads)))
        f.write("{0:26}\t{1}/{2} = {3:.4f}\n".format("Total Tied Reads:", 
                                              tied_reads, total_reads, 
                                              tied_reads/float(total_reads)))
        f.write("{0:26}\t{1}/{2} = {3:.4f}\n".format("Total Unassigned Reads:", 
                                      unassigned_reads, total_reads, 
                                      unassigned_reads/float(total_reads)))
        f.write("\n")

def init_int_stats(rep, repeat_edges, zero_it, position_path, partitioning, 
                   all_reads_file, template_len, cov, int_stats_file):
    #Count edge reads
    side_reads = {}
    total_reads = 0
    all_side_reads = 0
    internal_reads = 0
    for side in sorted(repeat_edges[rep]):
        part_list = _read_partitioning_file(partitioning.format(zero_it, side))
        total_reads = len(part_list)
        partitioning_outputs = _get_partitioning_info(part_list, 
                                                      repeat_edges[rep][side])
        side_reads[side], tied_reads, unassigned_reads = partitioning_outputs
        all_side_reads += sum(side_reads[side].values())
    internal_reads = total_reads - all_side_reads
    all_reads_n50 = _n50(all_reads_file)
    #Prepare header for iterative integrated stats
    #in/out Iter,Mean in/out/gap Len,Confirmed/Rejected Pos,Bridging Reads
    header_labels = []
    for side in sorted(repeat_edges[rep]):
        header_labels.extend(["{0} Iter".format(side)])
    header_labels.extend(["in Len", "Gap Len", "out Len"])
    header_labels.extend(["Confirmed", "Rejected"])
    side_edges = []
    for side in sorted(repeat_edges[rep]):
        side_edges.append([])
        for edge in sorted(repeat_edges[rep][side]):
            side_edges[-1].append("{0}{1}".format(side,edge))
    for edge_pair in sorted(product(*side_edges)):
        header_labels.extend(["{0}".format("|".join(edge_pair))])
    spaced_header = map("{:8}".format, header_labels)
    #Write to file
    with open(int_stats_file, 'w') as f:
        f.write("{0:16}\t{1}\n".format("Repeat:", rep))
        f.write("{0:16}\t{1}\n".format("Template Length:", template_len))
        f.write("{0:16}\t{1:.2f}\n".format("Avg Coverage:", cov))
        f.write("{0:16}\t{1}\n".format("# All Reads:", total_reads))
        f.write("{0:16}\t{1}\n\n".format("All Reads N50:", all_reads_n50))
        edge_headers = ["Side", "    Edge", "# Reads"]
        spaced_edge_header = map("{:5}".format, edge_headers)
        f.write("\t".join(spaced_edge_header))
        f.write("\n")
        for side in sorted(repeat_edges[rep]):
            for edge_id in sorted(repeat_edges[rep][side]):
                edge_values = [side, edge_id, side_reads[side][edge_id]]
                spaced_values = map("{:6}".format, edge_values)
                f.write("\t".join(spaced_values))
                f.write("\n")
        f.write("{0:12}\t  {1}\n".format("Internal", internal_reads))
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")
        

def update_int_stats(rep, repeat_edges, side_it, cons_align_path, template, 
                     template_len, confirmed_pos_path, int_confirmed_path, 
                     partitioning, int_stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    
    stats_out = []
    #Add side iters
    for side in sorted(repeat_edges[rep]):
        stats_out.extend([str(side_it[side])])
    #Find median in, out, and gap lengths
    medians = {s:0 for s in repeat_edges[rep]}
    for side in sorted(repeat_edges[rep]):
        trg_limits = []
        for edge_id in sorted(repeat_edges[rep][side]):
            curr_cons_path = cons_align_path.format(side_it[side], 
                                                    side, edge_id)
            if os.path.isfile(curr_cons_path):
                cons_align = _read_alignment(curr_cons_path, 
                                             template, 
                                             CONS_ALN_RATE)
            if  cons_align and cons_align[0]:
                if side == "in":
                    trg_limits.append(cons_align[0][0].trg_end)
                elif side == "out":
                    trg_limits.append(template_len - 
                                        cons_align[0][0].trg_start)
        if trg_limits:
            medians[side] = _get_median(trg_limits)
    gap_len = template_len - (medians["in"] + medians["out"])
    stats_out.extend([str(medians["in"]), str(gap_len), str(medians["out"])])
    #Add confirmed and rejected reads
    in_confirmed_path = confirmed_pos_path.format(side_it["in"], "in")
    out_confirmed_path = confirmed_pos_path.format(side_it["out"], "out")
    confirmed_pos_outputs = _integrate_confirmed_pos(in_confirmed_path,
                                                     out_confirmed_path)
    int_confirmed, int_rejected, pos = confirmed_pos_outputs
    _write_confirmed_positions(int_confirmed, int_rejected, pos, 
                               int_confirmed_path.format(side_it["in"], 
                                                         side_it["out"]))
    stats_out.extend([str(len(int_confirmed["total"])), 
                      str(len(int_rejected["total"]))])
    #Get bridging reads for each pair of in/out edges
    side_headers_dict = {}
    all_headers = set()
    for side in sorted(repeat_edges[rep]):
        side_headers_dict[side] = {}
        part_list = _read_partitioning_file(partitioning.format(side_it[side], 
                                                                side))
        for part in part_list:
            read_id, status, edge, top_sc, total_sc, header = part
            all_headers.add(header)
            if status == "Partitioned" and edge != "NA":
                side_headers_dict[side][header] = (side, int(edge))
    bridging_reads = {}
    side_edges = []
    for side in sorted(repeat_edges[rep]):
        side_edges.append([])
        for edge in sorted(repeat_edges[rep][side]):
            side_edges[-1].append((side, edge))
    for edge_pair in sorted(product(*side_edges)):
        bridging_reads[edge_pair] = 0
    for header in all_headers:
        if (header in side_headers_dict["in"] and 
            header in side_headers_dict["out"]):
            in_edge = side_headers_dict["in"][header]
            out_edge = side_headers_dict["out"][header]
            bridging_reads[(in_edge, out_edge)] += 1
    for edge_pair in sorted(bridging_reads):
        #stats_out.extend(["{0}".format(edge_pair)])
        stats_out.extend([str(bridging_reads[edge_pair])])
    spaced_header = map("{:8}".format, stats_out)
    #Write to file
    with open(int_stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")
        
def finalize_int_stats(rep, repeat_edges, side_it, cons_align_path, template, 
                       template_len, cons_vs_cons_path, consensuses, 
                       int_confirmed_path, partitioning, int_stats_file, 
                       resolved_seq_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MIN_BRIDGE_COUNT = trestle_config.vals["min_bridge_count"]
    MIN_BRIDGE_FACTOR = trestle_config.vals["min_bridge_factor"]
    
    #Resolved repeat seqs to be returned, NOT written
    resolved_repeats = {}
    summ_vals = []
    with open(int_stats_file, "a") as f:
        f.write("\n\n")
        for side in sorted(repeat_edges[rep]):
            f.write("{0}'{1}'{2:8}\t{3}\n".format("Final ", side, " Iter:", 
                                              side_it[side]))
        f.write("\n\n")
        #Overall confirmed and rejected positions
        int_confirmed, int_rejected, pos = _read_confirmed_positions(
                int_confirmed_path.format(side_it["in"], side_it["out"]))
        types = ["total", "sub", "del", "ins"]
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
        mean_position_gap = np.mean(position_gaps)
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
        #Write bridging reads
        side_headers_dict = {}
        all_headers = set()
        for side in sorted(repeat_edges[rep]):
            side_headers_dict[side] = {}
            part_list = _read_partitioning_file(partitioning.format(
                                                    side_it[side], side))
            for part in part_list:
                read_id, status, edge, top_sc, total_sc, header = part
                all_headers.add(header)
                if status == "Partitioned" and edge != "NA":
                    side_headers_dict[side][header] = (side, int(edge))
        bridging_reads = {}
        side_edges = []
        for side in sorted(repeat_edges[rep]):
            side_edges.append([])
            for edge in sorted(repeat_edges[rep][side]):
                side_edges[-1].append((side, edge))
        for edge_pair in sorted(product(*side_edges)):
            bridging_reads[edge_pair] = 0
        for header in all_headers:
            if (header in side_headers_dict["in"] and 
                header in side_headers_dict["out"]):
                in_edge = side_headers_dict["in"][header]
                out_edge = side_headers_dict["out"][header]
                bridging_reads[(in_edge, out_edge)] += 1
        for edge_pair in sorted(bridging_reads):
            pair_label = "|".join(["".join(map(str, x)) for x in edge_pair])
            f.write("{0}{1:21}\t{2}\n".format(pair_label, " Bridging Reads:",
                                           bridging_reads[edge_pair]))
        f.write("\n\n")
        #Write combos which are sets of bridging reads
        all_combos = _get_combos(side_edges[0], side_edges[1])
        combo_support = [0 for _ in all_combos]
        for i, combo in enumerate(all_combos):
            for edge_pair in combo:
                if edge_pair in bridging_reads:
                    combo_support[i] += bridging_reads[edge_pair]
        for i, combo in enumerate(all_combos):
            f.write("{0} {1}\n".format("Combo", i))
            coms = ["|".join(["".join(map(str, x)) for x in y]) for y in combo]
            combo_edges = " + ".join(coms)
            f.write("{0:12}\t{1}\n".format("Resolution:", combo_edges))
            f.write("{0:12}\t{1}\n\n".format("Support:", combo_support[i]))
        #Bridging conditions 
        bridged = False
        bridged_edges = None
        combo_inds = zip(combo_support, range(len(combo_support)))
        sorted_combos = sorted(combo_inds, reverse=True)
        if (sorted_combos[0][0] >= MIN_BRIDGE_COUNT and 
            sorted_combos[0][0] >= sorted_combos[1][0] * MIN_BRIDGE_FACTOR):
            bridged = True
            bridged_edges = all_combos[sorted_combos[0][1]]
        best_combo = sorted_combos[0][1]
        best_support = sorted_combos[0][0]
        best_against = 0
        for support, ind in sorted_combos[1:]:
            best_against += support
        second_combo = sorted_combos[1][1]
        second_support = sorted_combos[1][0]
        if bridged:
            f.write("BRIDGED\n")
            f.write("Bridging Combo: {0}\n".format(best_combo))
            br_ct_str = "{0} (min_bridge_count)".format(MIN_BRIDGE_COUNT)
            br_diff_str = "{0} * {1} (Combo {2} * min_bridge_factor)".format(
                second_support, MIN_BRIDGE_FACTOR, second_combo)
            f.write("Support = {0}\t> {1}\n{2:12}\t> {3}\n".format(
                best_support, br_ct_str, "", br_diff_str))
            f.write("Resolution:\n")
            for edge_pair in sorted(bridged_edges):
                f.write("{0[0]} {0[1]:2}  {1:3} {2[0]} {2[1]}\n".format(
                                             edge_pair[0], 
                                             "->", 
                                             edge_pair[1]))
            f.write("\n\n")
        else:
            f.write("UNBRIDGED\n")
            f.write("Best combo {0}\n".format(best_combo))
            f.write("{0:20}\t{1}\n".format("min_bridge_count", 
                                           MIN_BRIDGE_COUNT))
            f.write("{0:20}\t{1}\n\n\n".format("min_bridge_factor", 
                                               MIN_BRIDGE_FACTOR))
        summ_vals.extend([bridged, best_support, best_against])
        #If not bridged, find in/gap/out lengths and divergence rates
        if not bridged:
            #Write median in, out, and gap lengths
            side_lens = {s:0 for s in repeat_edges[rep]}
            for side in sorted(repeat_edges[rep]):
                trg_limits = []
                for edge_id in sorted(repeat_edges[rep][side]):
                    curr_cons_path = cons_align_path.format(side_it[side],
                                                             side, edge_id)
                    if os.path.isfile(curr_cons_path):
                        cons_align = _read_alignment(
                                curr_cons_path, 
                                template, 
                                CONS_ALN_RATE)
                        if cons_align and cons_align[0]:
                            if side == "in":
                                trg_limits.append(cons_align[0][0].trg_end)
                            elif side == "out":
                                trg_limits.append(template_len - 
                                                cons_align[0][0].trg_start)
                if trg_limits:
                    side_lens[side] = _get_median(trg_limits)
            gap_len = template_len - (side_lens["in"] + side_lens["out"])
            f.write("{0:30}\t{1}\n".format("Median in Sequence Length:", 
                                           side_lens["in"]))
            f.write("{0:30}\t{1}\n".format("Median out Sequence Length:", 
                                           side_lens["out"]))
            f.write("{0:30}\t{1}\n\n".format("Median Gap/Overlap Length:", 
                                             gap_len))
            
            #Write mean in and out divergence rates
            div_rates = {s:[] for s in repeat_edges[rep]}
            for side in sorted(repeat_edges[rep]):
                side_pairs = sorted(combinations(repeat_edges[rep][side], 2))
                for edge_one, edge_two in side_pairs:
                    cons_cons_file = cons_vs_cons_path.format(
                                        side_it[side], side, edge_one, 
                                        side_it[side], side, edge_two)
                    if os.path.isfile(cons_cons_file):
                        cons_vs_cons = _read_alignment(
                                cons_cons_file, 
                                consensuses.format(side, edge_two), 
                                CONS_ALN_RATE)
                        if cons_vs_cons and cons_vs_cons[0]:
                            edge_rate = _calculate_divergence(
                                            cons_vs_cons[0][0].qry_seq, 
                                            cons_vs_cons[0][0].trg_seq)
                            div_rates[side].append(edge_rate)
            mean_in_div = 0.0
            if div_rates["in"]:
                mean_in_div = np.mean(div_rates["in"])
            mean_out_div = 0.0
            if div_rates["out"]:
                mean_out_div = np.mean(div_rates["out"])
            weighted_mean_div = 0.0
            if side_lens["in"] + side_lens["out"] != 0:
                weighted_mean_div = (mean_in_div*side_lens["in"] + 
                                     mean_out_div*side_lens["out"]) / float(
                                     side_lens["in"] + side_lens["out"])
            f.write("{0:30}\t{1}\n".format("Mean in Divergence Rate:", 
                                            mean_in_div))
            f.write("{0:30}\t{1}\n".format("Mean out Divergence Rate:", 
                                            mean_out_div))
            f.write("{0:30}\t{1}\n\n".format("Weighted Mean Divergence Rate:", 
                                          weighted_mean_div))
            res_str = "No resolution so no resolved file for repeat {0}\n\n"
            f.write(res_str.format(rep))
            for i, edge in enumerate(sorted(repeat_edges[rep]["in"])):
                header = "Repeat_{0}_unbridged_copy_{1}".format(rep, i)
                resolved_repeats[header] = ""
                #seq_dict = {header:""}
                #fp.write_fasta_dict(seq_dict, resolved_seq_file.format(i))
            summ_vals.extend([""])
        #If bridged, find overlap and construct repeat copy sequences
        else:
            #Find end of repeat as min/max of in/out cons_vs_cons alignments
            edge_limits = {}
            for side in sorted(repeat_edges[rep]):
                side_pairs = sorted(combinations(repeat_edges[rep][side], 2))
                for edge_one, edge_two in side_pairs:
                    cons_cons_file = cons_vs_cons_path.format(
                                        side_it[side], side, edge_one, 
                                        side_it[side], side, edge_two)
                    if os.path.isfile(cons_cons_file):
                        cons_vs_cons = _read_alignment(
                            cons_cons_file, 
                            consensuses.format(side, edge_two), 
                            CONS_ALN_RATE)
                    if cons_vs_cons and cons_vs_cons[0]:
                        one_start = cons_vs_cons[0][0].qry_start
                        one_end = cons_vs_cons[0][0].qry_end
                        two_start = cons_vs_cons[0][0].trg_start
                        two_end = cons_vs_cons[0][0].trg_end
                        if side == "in":
                            if (side, edge_one) not in edge_limits:
                                edge_limits[(side, edge_one)] = one_start
                            elif one_start < edge_limits[(side, edge_one)]:
                                edge_limits[(side, edge_one)] = one_start
                            if (side, edge_two) not in edge_limits:
                                edge_limits[(side, edge_two)] = two_start
                            elif two_start < edge_limits[(side, edge_two)]:
                                edge_limits[(side, edge_two)] = two_start
                        elif side == "out":
                            if (side, edge_one) not in edge_limits:
                                edge_limits[(side, edge_one)] = one_end
                            elif one_end > edge_limits[(side, edge_one)]:
                                edge_limits[(side, edge_one)] = one_end
                            if (side, edge_two) not in edge_limits:
                                edge_limits[(side, edge_two)] = two_end
                            elif two_end > edge_limits[(side, edge_two)]:
                                edge_limits[(side, edge_two)] = two_end
            #For each edge_pair, find starting and ending indices of
            #in, out, and template sequences to construct sequences
            summ_resolution = []
            for i, edge_pair in enumerate(sorted(bridged_edges)):
                f.write("Repeat Copy {0}\n".format(i))
                f.write("{0[0]} {0[1]:2}  {1:3} {2[0]} {2[1]}\n".format(
                                             edge_pair[0], 
                                             "->", 
                                             edge_pair[1]))
                in_start = None
                out_end = None
                out_align = None
                in_align = None
                for side, edge_id in edge_pair:
                    if side == "in" and (side, edge_id) in edge_limits:
                        in_start = edge_limits[(side, edge_id)]
                    elif side == "out" and (side, edge_id) in edge_limits:
                        out_end = edge_limits[(side, edge_id)]
                    if os.path.isfile(cons_align_path.format(side_it[side], 
                                                             side, 
                                                             edge_id)):
                        cons_align = _read_alignment(
                                cons_align_path.format(side_it[side], 
                                                       side, 
                                                       edge_id), 
                                template, 
                                CONS_ALN_RATE)
                        if cons_align and cons_align[0]:
                            if side == "in":
                                in_align = cons_align[0][0]
                            elif side == "out":
                                out_align = cons_align[0][0]
                if not in_align:
                    in_end = 0
                    temp_start = 0
                    if in_start is None:
                        in_start = 0
                else:
                    in_end = in_align.qry_end
                    temp_start = in_align.trg_end
                    if in_start is None:
                        in_start = in_align.qry_start
                if not out_align:
                    temp_end = 0
                    out_start = 0
                    if out_end is None:
                        out_end = 0
                    out_qry_seq = ""
                    out_trg_seq = ""
                    out_trg_end = 0
                    out_qry_end = 0
                else:
                    temp_end = out_align.trg_start
                    out_start = out_align.qry_start
                    if out_end is None:
                        out_end = out_align.qry_end
                    out_qry_seq = out_align.qry_seq
                    out_trg_seq = out_align.trg_seq
                    out_trg_end = out_align.trg_end
                    out_qry_end = out_align.qry_end
                f.write("Alignment Indices:\n")
                f.write("{0:10}\t{1:5} - {2:5}\n".format("in", 
                                                         in_start, in_end))
                #f.write("{0:10}\t{1:5} - {2:5}\n".format("Template", 
                                                          #temp_start, 
                                                          #temp_end))
                f.write("{0:10}\t{1:5} - {2:5}\n".format("out", 
                                                         out_start, out_end))
                #Report gap/overlap length
                gap_len = temp_end - temp_start
                if gap_len >= 0:
                    f.write("{0}\t{1}\n".format("Gap between edges:", gap_len))
                else:
                    f.write("{0}\t{1}\n\n".format("Overlap between edges:", 
                                                  -gap_len))
                    #in sequence used to represent overlapping segment
                    #print check of overlapping segment
                    new_temp_end = temp_start
                    new_out_start = None
                    out_qry_aln, out_aln_qry = _index_mapping(out_qry_seq)
                    out_trg_aln, out_aln_trg = _index_mapping(out_trg_seq)
                    
                    in_edge = edge_pair[0][1]
                    out_edge = edge_pair[1][1]
                    if temp_start >= out_trg_end:
                        new_out_start = out_qry_end
                    else:
                        if temp_start < len(out_trg_aln):
                            out_aln_ind = out_trg_aln[temp_start]
                            if out_aln_ind < len(out_aln_qry):
                                new_out_start = (out_start + 
                                                 out_aln_qry[out_aln_ind])
                    """_check_overlap(
                            consensuses.format("in", in_edge), 
                            template,
                            consensuses.format("out", out_edge), 
                            -gap_len, in_start, in_end, temp_start, temp_end, 
                            out_start, out_end,
                            new_out_start, in_align.qry_seq, in_align.trg_seq, 
                            out_align.qry_seq, out_align.trg_seq, out_trg_aln, 
                            out_aln_trg, out_qry_aln, out_aln_qry, 
                            out_align.trg_end, out_align.qry_end, 
                            in_align, out_align)
                    """                    
                    temp_end = new_temp_end
                    if new_out_start:
                        out_start = new_out_start
                    f.write("Adjusted Alignment Indices:\n")
                    f.write("{0:10}\t{1:5} - {2:5}\n".format("in", 
                                                        in_start, in_end))
                    if temp_start != new_temp_end:
                        f.write("{0:10}\t{1:5} - {2:5}\n".format("Template", 
                                                        temp_start, 
                                                        new_temp_end))
                    f.write("{0:10}\t{1:5} - {2:5}\n\n\n".format("out", 
                                                        new_out_start, 
                                                        out_end))
                    
                in_edge = edge_pair[0][1]
                out_edge = edge_pair[1][1]
                header = "_".join(["Repeat_{0}".format(rep), 
                                        "bridged_copy_{0}".format(i), 
                                        "in_{0}_{1}_{2}".format(in_edge, 
                                                                in_start, 
                                                                in_end), 
                                        "template_{0}_{1}".format(temp_start, 
                                                                  temp_end), 
                                        "out_{0}_{1}_{2}".format(out_edge, 
                                                                 out_start, 
                                                                 out_end)])
                copy_seq = _construct_repeat_copy(
                        consensuses.format("in", in_edge), 
                        template,
                        consensuses.format("out", out_edge), 
                        in_start, in_end, 
                        temp_start, temp_end, 
                        out_start, out_end)
                resolved_repeats[header] = copy_seq
                if copy_seq:
                    seq_dict = {header:copy_seq}
                    fp.write_fasta_dict(seq_dict, 
                                        resolved_seq_file.format(rep, i))
                in_str = "".join(["in", str(in_edge)])
                out_str = "".join(["out", str(out_edge)])
                summ_resolution.append("|".join([in_str, out_str]))
            summ_vals.extend(["+".join(summ_resolution)])
    return bridged, resolved_repeats, summ_vals

def int_stats_postscript(rep, repeat_edges, integrated_stats, 
                         resolved_rep_path, res_vs_res):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    
    divs = []
    with open(integrated_stats, "a") as f:
        res_inds = range(len(repeat_edges[rep]["in"]))
        f.write("Resolved Repeat Sequence Alignments\n")
        for res_one, res_two in sorted(combinations(res_inds, 2)):
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            qry_seq = ""
            trg_seq = ""
            if os.path.isfile(res_vs_res.format(rep, res_one, res_two) and
                resolved_rep_path.format(rep, res_two)):
                res_align = _read_alignment(res_vs_res.format(rep, res_one, 
                                                              res_two), 
                                            resolved_rep_path.format(rep, 
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
            f.write("Copy {0}|Copy {1}\n".format(res_one, res_two))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Copy ", res_one, ":", 
                    qry_start, qry_end, qry_len))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Copy ", res_two, ":",   
                    trg_start, trg_end, trg_len))
            div_rate = _calculate_divergence(qry_seq, trg_seq)
            divs.append(div_rate)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
    return np.mean(divs)
        
def _get_partitioning_info(part_list, edges):
    edge_reads = {edge:0 for edge in edges}
    tied_reads = 0
    unassigned_reads = 0
    for part in part_list:
        read_id, status, edge, top_sc, total_sc, header = part
        if status == "Partitioned" and edge != "NA":
            edge_reads[int(edge)] += 1
        elif status == "Tied":
            tied_reads += 1
        elif status == "None":
            unassigned_reads += 1
        else:
            exception_str = "Unknown status {0} in partitioning file"
            raise Exception(exception_str.format(status))
    return edge_reads, tied_reads, unassigned_reads

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
    reads_dict = fp.read_fasta_dict(reads_file)
    read_lengths = sorted([len(x) for x in reads_dict.values()], reverse=True)
    summed_len = 0
    n50 = 0
    for l in read_lengths:
        summed_len += l
        if summed_len >= sum(read_lengths)/2:
            n50 = l
            break
    return n50

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

def _integrate_confirmed_pos(in_confirmed_path, out_confirmed_path):
    in_conf, in_rej, in_pos = _read_confirmed_positions(in_confirmed_path)
    out_conf, out_rej, out_pos = _read_confirmed_positions(out_confirmed_path)
    
    integrated_confirmed = {"total":[], "sub":[], "ins":[], "del":[]}
    integrated_rejected = {"total":[], "sub":[], "ins":[], "del":[]}
    
    for pos in sorted(in_pos["total"]):
        for pos_type in in_conf:
            if pos in in_conf[pos_type] or pos in out_conf[pos_type]:
                integrated_confirmed[pos_type].append(pos)
            elif pos in in_rej[pos_type] or pos in out_rej[pos_type]:
                integrated_rejected[pos_type].append(pos)
    return integrated_confirmed, integrated_rejected, in_pos

def _get_combos(in_list, out_list):
    all_combos = []
    for combo in _combo_helper(in_list, out_list):
        all_combos.append(combo)
    return all_combos

def _combo_helper(in_list, out_list):
    if not in_list or not out_list:
        yield []
        return
    else:
        in1 = in_list[0]
        for j in range(len(out_list)):
            combo = (in1, out_list[j])
            for rest in _combo_helper(in_list[1:], 
                                      out_list[:j] + out_list[j + 1:]):
                 yield [combo] + rest

def _get_aln_end(aln_start, aln_seq):
    return aln_start+len(aln_seq.replace("-",""))

def _check_overlap(in_file, temp_file, out_file, overlap, in_start, in_end, 
                   temp_start, temp_end, out_start, out_end, new_out_start, 
                   in_qry, in_trg, out_qry, out_trg, out_trg_aln, out_aln_trg, 
                   out_qry_aln, out_aln_qry, out_trg_end, out_qry_end, 
                   in_align, out_align):
    in_dict = fp.read_fasta_dict(in_file)
    in_seq = in_dict.values()[0]
    temp_dict = fp.read_fasta_dict(temp_file)
    temp_seq = temp_dict.values()[0]
    out_dict = fp.read_fasta_dict(out_file)
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


def _construct_repeat_copy(in_file, temp_file, out_file, in_start, in_end, 
                           temp_start, temp_end, out_start, out_end):
    if (not os.path.isfile(in_file) or 
        not os.path.isfile(temp_file) or 
        not os.path.isfile(out_file)):
        return ""
    in_dict = fp.read_fasta_dict(in_file)
    in_seq = in_dict.values()[0]
    temp_dict = fp.read_fasta_dict(temp_file)
    temp_seq = temp_dict.values()[0]
    out_dict = fp.read_fasta_dict(out_file)
    out_seq = out_dict.values()[0]
    seq = ''.join([in_seq[in_start:in_end], 
                   temp_seq[temp_start:temp_end], 
                   out_seq[out_start:out_end]])
    return seq

def init_summary(summary_file):
    with open(summary_file, "w") as f:
        summ_header_labels = ["Repeat", "Template", "Cov", "# Confirmed Pos", 
                              "Max Pos Gap", "Bridged?", "Support", "Against", 
                              "Avg Div", "Resolution"]
        spaced_header = map("{:13}".format, summ_header_labels)
        f.write("\t".join(spaced_header))
        f.write("\n")

def update_summary(rep, template_len, avg_cov, summ_vals, avg_div, 
                   summary_file):
    (confirmed_pos, max_pos_gap, bridged, 
     support, against, resolution) = tuple(summ_vals)
    summ_out = [rep, template_len, avg_cov, confirmed_pos, max_pos_gap, 
                bridged, support, against, avg_div, resolution]
    summ_out[2] = "{:.4}".format(summ_out[2])
    summ_out[5] = str(summ_out[5])
    summ_out[8] = "{:.4}".format(summ_out[8])
    spaced_summ = map("{:13}".format, map(str, summ_out))
    with open(summary_file, "a") as f:
        f.write("\t".join(spaced_summ))
        f.write("\n")
