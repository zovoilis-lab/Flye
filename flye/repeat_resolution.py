# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 03:50:31 2017

@author: jeffrey_yuan
"""

import logging
from collections import defaultdict
from itertools import izip, combinations, product
import multiprocessing
import signal
import copy
import os
import numpy as np

from flye.alignment import shift_gaps, SynchronizedSamReader
import flye.config as config
import flye.fasta_parser as fp
import flye.divergence as div

logger = logging.getLogger()

#Process Repeats functions
class ProcessingException(Exception):
    pass
       
def process_repeats(reads, repeats_dump, graph_edges, work_dir, 
                    template_file, extended_file, 
                    repeat_reads_file, partitioning_file, repeat_label, 
                    orient_labels, min_mult, max_mult, extend_len):
    """Generates repeat dirs and files given reads, repeats_dump and
    graph_edges files. Only returns repeats between min_mult and max_mult"""
    
    #Reads input files
    repeats_dict = _read_repeats_dump(repeats_dump)
    reads_dict = {}
    for read_file in reads:
        reads_dict.update(fp.read_fasta_dict(read_file))
    orig_graph = fp.read_fasta_dict(graph_edges)
    graph_dict = {int(h.split('_')[1]):orig_graph[h] for h in orig_graph}
    
    if not repeats_dict:
        logger.debug("Empty repeats_dump file: {0}".format(repeats_dump))
    if not reads_dict:
        raise ProcessingException("No reads found from {0}".format(reads))
    if not graph_dict:
        raise ProcessingException("No edges found from {0}".format(graph_edges))
        
    
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
        if (repeats_dict[rep][0] < min_mult or
                 repeats_dict[rep][0] > max_mult or
                 repeats_dict[-rep][0] < min_mult or
                 repeats_dict[-rep][0] > max_mult):
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
            template_path = os.path.join(orient_path, template_file)
            extended_path = os.path.join(orient_path, extended_file)
            repeat_reads_path = os.path.join(orient_path, repeat_reads_file)
            partitioning_path = os.path.join(orient_path, partitioning_file)
            
            in_label = "in"
            out_label = "out"
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
            """Headers will be in the form -h or +h,
            edge_dict is in the form >[Input,Output]_edge##_h,
            rev_comp of read will be written if the header is -h"""
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
                                                       total_score, header[1:]))
                        all_edge_headers[curr_rep][header[1:]] = read_id
                        read_id += 1
                
                extend_in_header = "Extended_Template_Input_{0}".format(edge_id)
                if edge_id > 0:
                    edge_seq = graph_dict[edge_id]
                elif edge_id < 0:
                    edge_seq = fp.reverse_complement(graph_dict[-edge_id])
                extended_seq = edge_seq[-extend_len:]
                extended_dicts[(in_label, edge_id)][extend_in_header] = (
                                        extended_seq + template_seq)
                
            
            for edge_id in outputs_dict:         
                repeat_edges[curr_rep]["out"].append(edge_id)
                extended_dicts[("out", edge_id)] = {}
                
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
                                                       total_score, header[1:]))
                        
                        all_edge_headers[curr_rep][header[1:]] = read_id
                        read_id += 1
                    else:
                        curr_read_id = all_edge_headers[curr_rep][header[1:]]
                    
                    if header[1:] not in out_headers:
                        status_label = "Partitioned"
                        edge_label = str(edge_id)
                        score = 1
                        total_score = 0
                        partitioning[out_label].append((curr_read_id, status_label, 
                                                        edge_label, score, 
                                                        total_score, header[1:]))
                        out_headers.add(header[1:])   
                
                extend_out_header = "Extended_Template_Output_{0}".format(edge_id)
                if edge_id > 0:
                    edge_seq = graph_dict[edge_id]
                elif edge_id < 0:
                    edge_seq = fp.reverse_complement(graph_dict[-edge_id])
                extended_seq = edge_seq[:extend_len]
                extended_dicts[(out_label, edge_id)][extend_out_header] = (
                                        template_seq + extended_seq)
            
            #Need to reiterate over in_headers to add in_headers to out-partitioning
            #while avoiding double-adding ones in both
            for edge_id in inputs_dict:
                headers = inputs_dict[edge_id]
                for header in headers:
                    if header[1:] not in out_headers:
                        curr_read_id = all_edge_headers[curr_rep][header[1:]]
                        status_label = "None"
                        edge_label = "NA"
                        score = 0
                        total_score = 0
                        partitioning[out_label].append((curr_read_id, status_label, 
                                                        edge_label, score, 
                                                        total_score, header[1:]))
                
            
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
            
            
            fp.write_fasta_dict(template_dict, template_path)
            for edge in extended_dicts:
                extended_edge_path = extended_path.format(edge[0], edge[1])
                fp.write_fasta_dict(extended_dicts[edge], extended_edge_path)
            fp.write_fasta_dict(repeat_reads_dict, repeat_reads_path)
            for side in [in_label, out_label]:
                _write_partitioning_file(partitioning[side], 
                                         partitioning_path.format(side))
            
            if not template_dict:
                raise ProcessingException("No template {0} found".format(curr_rep))
            for edge in extended_dicts:
                if not template_dict:
                    raise ProcessingException(
                        "No extended template {0} {1} {2} found".format(
                            curr_rep, edge[0], edge[1]))
            if not repeat_reads_dict:
                raise ProcessingException("No repeat reads {0} found".format(curr_rep))
            for side in [in_label, out_label]:
                if not partitioning[side]:
                    raise ProcessingException(
                        "Empty partitioning file {0}".format(
                            partitioning_path.format(side)))
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
    fp.write_fasta_dict(edge_reads, out_file)

#Partition Reads Functions

class PartitioningException(Exception):
    pass

def partition_reads(edges, it, side, position_path, cons_align_path, 
                    template, read_align_path, consensuses, 
                    confirmed_pos_path, part_file, 
                    headers_to_id, min_aln_rate, buffer_count,
                    pos_test, num_test):
    pos = div.read_positions(position_path)
    """"""
    test_pos = pos[:num_test]
    with open(pos_test, 'w') as f:
        f.write("Test positions ({0})\n".format(num_test))
        f.write("{0}\n\n".format(test_pos))
        
    cons_aligns = {}
    for edge_id in edges:
        cons_aligns[edge_id] = _read_alignment(cons_align_path.format(it, side, edge_id), 
                                               template, 
                                               min_aln_rate)
        if (not cons_aligns or 
                not cons_aligns[edge_id] or 
                not cons_aligns[edge_id][0]):
            raise PartitioningException("No alignment found for edge {0}".format(edge_id))
            
    confirmed_pos, rejected_pos, consensus_pos = _evaluate_positions(pos, cons_aligns, 
                                                                     pos_test, test_pos, 
                                                                     side)
    _write_confirmed_positions(confirmed_pos, rejected_pos, pos, confirmed_pos_path)
    read_aligns = {}
    for edge_id in edges:
        read_aligns[edge_id] = _read_alignment(read_align_path.format(it, side, edge_id), 
                                               consensuses[(it, side, edge_id)], 
                                               min_aln_rate)
    _classify_reads(read_aligns, consensus_pos, part_file, 
                    headers_to_id, buffer_count, pos_test, num_test)
    
def _read_alignment(alignment, target_path, min_aln_rate):
    alignments = []
    aln_reader = SynchronizedSamReader(alignment,
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

    def __init__(self, edge_id, qry_seq, trg_seq, qry_start, trg_start, trg_end):
        self.edge_id = edge_id
        self.qry_seq = shift_gaps(trg_seq, qry_seq)
        self.trg_seq = shift_gaps(self.qry_seq, trg_seq)
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

def _evaluate_positions(positions, cons_aligns, pos_test, test_pos, side):
    #Includes insertions!
    pos = set(positions)
    confirmed_pos = []
    rejected_pos = []
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
                
        
        if trg_ind in pos:
            if ((side == "in" and trg_ind < min_end_edge) or
                (side == "out" and trg_ind >= max_start_edge)):
                """"""
                if trg_ind in test_pos:
                    with open(pos_test, "a") as f:
                        f.write("Pos {0}:\n".format(trg_ind))
                ins_confirmed = False
                qry_confirmed = False
                qry_nuc = ""
                trg_nuc = ""
                for edge_id in alns:
                    aln = alns[edge_id]
                    """"""
                    if trg_ind in test_pos:
                        with open(pos_test, "a") as f:
                            f.write("Edge {0} in_alignment:\t{1}\n".format(edge_id, aln.in_alignment))
                    if aln.in_alignment:
                        #Directly add positions only to consensuses where insertions occur
                        #Add the position prior to curr_qry_ind to account for insertion
                        if aln.curr_ins_nuc:
                            ins_confirmed = True
                            consensus_pos[edge_id].append(aln.curr_qry_ind-1)
                            """"""
                            if trg_ind in test_pos:
                                with open(pos_test, "a") as f:
                                    f.write("ins_confirmed, cons_pos:\t{0}\n".format(aln.curr_qry_ind-1))
                            
                        if qry_nuc and qry_nuc != aln.curr_qry_nuc:
                            qry_confirmed = True
                        else:
                            qry_nuc = aln.curr_qry_nuc
                        if trg_nuc and trg_nuc != aln.curr_trg_nuc:
                            raise PartitioningException("Inconsistent trg_nuc, {0} {1} {2} {3}".format(edge_id, trg_ind, trg_nuc, aln.curr_trg_nuc))
                        else:
                            trg_nuc = aln.curr_trg_nuc
                        """"""
                        if trg_ind in test_pos:
                            with open(pos_test, "a") as f:
                                f.write("aln_ind:\t{0}\n".format(aln.curr_aln_ind))
                                f.write("qry_ind:\t{0}\n".format(aln.curr_qry_ind))
                                f.write("qry_nuc:\t{0}\n".format(aln.curr_qry_nuc))
                                f.write("trg_nuc:\t{0}\n".format(aln.curr_trg_nuc))
                                f.write("ins_nuc:\t{0}\n".format(aln.curr_ins_nuc))
                                f.write("around qry_seq:\t{0}\n".format(aln.qry_seq[aln.curr_aln_ind-5:aln.curr_aln_ind+5]))
                                f.write("around trg_seq:\t{0}\n".format(aln.trg_seq[aln.curr_aln_ind-5:aln.curr_aln_ind+5]))
                                f.write("current trg_nuc:\t{0}\n".format(trg_nuc))
                                f.write("current ins_confirmed:\t{0}\n".format(ins_confirmed))
                                f.write("current qry_confirmed:\t{0}\n\n".format(qry_confirmed))
                if ins_confirmed or qry_confirmed:
                    confirmed_pos.append(trg_ind)
                    #Add positions to consensuses for only substitutions/deletions
                    if qry_confirmed:
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
                                        f.write("Edge {0}, cons_pos:\t{1}\n".format(edge_id, aln.curr_qry_ind))
                else:
                    rejected_pos.append(trg_ind)
        
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
        f.write(">Confirmed_{0}_positions\n".format(len(confirmed)))
        f.write(",".join(map(str, sorted(confirmed))))
        f.write("\n")
        f.write(">Rejected_{0}_positions\n".format(len(rejected)))
        f.write(",".join(map(str, sorted(rejected))))
        f.write("\n")
        f.write(">{0}_tentative_positions\n".format(len(pos)))
        f.write(",".join(map(str, sorted(pos))))
        f.write("\n")

def _read_confirmed_positions(confirmed_file):
    confirmed = []
    rejected = []
    pos = []
    with open(confirmed_file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i == 1:
                confirmed = map(int, line.split(","))
            if i == 3:
                rejected = map(int, line.split(","))
            if i == 5:
                pos = map(int, line.split(","))
    return confirmed, rejected, pos
                

def _classify_reads(read_aligns, consensus_pos, part_file, 
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
                            f.write("\nEdge {0}, cons_pos {1}, read {2}\n".format(edge_id, pos, headers_to_id[read_header]))
                            f.write("pos_minus_start:\t{0}\n".format(pos_minus_start))
                            f.write("qry_seq:\t{0}\n".format(aln.qry_seq[aln_ind]))
                            f.write("trg_seq:\t{0}\n".format(aln.trg_seq[aln_ind]))
                            f.write("around qry_seq:\t{0}\n".format(aln.qry_seq[aln_ind-5:aln_ind+5]))
                            f.write("around trg_seq:\t{0}\n".format(aln.trg_seq[aln_ind-5:aln_ind+5]))
    
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
                if edge_score > top_score:
                    top_edge = edge_id
                    top_score = edge_score
                    total_score += edge_score
                    tie_bool = False
                elif edge_score == top_score:
                    tie_bool = True
                    total_score += edge_score
                else:
                    total_score += edge_score
            
            if tie_bool:
                status_label = "Tied"
                edge_label = "NA"
                partitioning.append((read_id, status_label, edge_label, 
                                    top_score, total_score, read_header))
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
    
    _write_partitioning_file(partitioning, part_file)

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

def init_side_stats(rep, side, repeat_edges, args, buffer_count, max_iter, 
                     min_edge_cov, position_path, partitioning, 
                     prev_parts, template_len, stats_file):
    pos = div.read_positions(position_path)
    #Count partitioned reads
    edge_below_cov = False
    part_list = _read_partitioning_file(partitioning)
    edge_reads, tied_reads, unassigned_reads = _get_partitioning_info(
                                                        part_list, 
                                                        repeat_edges[rep][side])
    #Check break condition for iteration loop
    for edge in repeat_edges[rep][side]:
        if edge_reads[edge] < min_edge_cov:
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
        f.write("{0:20}\t{1}\n".format("Repeat:", rep))
        f.write("{0:20}\t{1}\n".format("Side:", side))
        f.write("{0:20}\t".format("Edges:"))
        f.write(", ".join(map(str,sorted(repeat_edges[rep][side]))))
        f.write("\n")
        f.write("{0:20}\t{1}\n\n".format("Template Length:", template_len))
        f.write("Initial option values\n")
        f.write("{0:20}\t{1}\n".format("min_overlap:", args.min_overlap))
        f.write("{0:20}\t{1}\n".format("sub_thresh:", args.sub_thresh))
        f.write("{0:20}\t{1}\n".format("del_thresh:", args.del_thresh))
        f.write("{0:20}\t{1}\n".format("ins_thresh:", args.ins_thresh))
        f.write("{0:20}\t{1}\n".format("extend_len:", args.extend_len))
        f.write("{0:20}\t{1}\n".format("buffer_count:", buffer_count))
        f.write("{0:20}\t{1}\n".format("max_iter:", max_iter))
        f.write("{0:20}\t{1}\n".format("min_edge_cov:", min_edge_cov))
        f.write("\n")
        f.write("{0:20}\t{1}\n".format("Tentative positions:", len(pos)))
        f.write("{0:20}\t{1}\n".format("Total reads:", sum(edge_reads.values())))
        for edge in sorted(repeat_edges[rep][side]):
            f.write("{0}{1}{2:13}\t{3}\n".format("Edge ", edge, 
                                                 " starting_reads:", 
                                                 edge_reads[edge]))
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")
    
    return edge_below_cov

def update_side_stats(edges, it, side, cons_align_path, template, 
                      min_aln_rate, confirmed_pos_path, partitioning, 
                      min_edge_cov, prev_parts, stats_file):
    #Write stats for each iteration
    #Iter,Rep Lens,Confirmed/Rejected Pos,Partitioned Reads
    stats_out = [str(it)]
    for edge_id in sorted(edges):
        cons_align = _read_alignment(cons_align_path.format(it, side, edge_id), 
                                     template, 
                                     min_aln_rate)
        rep_len = 0
        if side == "in":
            rep_len = cons_align[0][0].qry_len - cons_align[0][0].qry_start
        elif side == "out":
            rep_len = cons_align[0][0].qry_end
        stats_out.extend([str(rep_len)])
    confirmed, rejected, pos = _read_confirmed_positions(confirmed_pos_path)
    stats_out.extend([str(len(confirmed)), str(len(rejected))])
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
        if edge_reads[edge] < min_edge_cov:
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

def finalize_side_stats(edges, it, side, cons_align_path, template, min_aln_rate, 
                 cons_vs_cons_path, consensuses, confirmed_pos_path, 
                 partitioning, max_iter, edge_below_cov, dup_part, 
                 stats_file):
    with open(stats_file, "a") as f:
        f.write("\n\n")
        f.write("{0:26}\t{1}\n".format("Final Iter:", it))
        f.write("Iteration terminated because:\n")
        if it == max_iter:
            f.write("Max iter reached\n")
        if edge_below_cov:
            f.write("Edge coverage fell below min_edge_cov\n")
        if dup_part:
            f.write("Partitioning was identical to a previous iteration\n")
        f.write("\n")
        #Write out alignment indices for edges vs template
        limit_ind = None
        limit_label = ""
        if side == "in":
            limit_label = "Min Template End"
        elif side == "out":
            limit_label = "Max Template Start"
        for edge_id in sorted(edges):
            cons_align = _read_alignment(cons_align_path.format(it, side, edge_id), 
                                     template, 
                                     min_aln_rate)
            trg_start = cons_align[0][0].trg_start
            trg_end = cons_align[0][0].trg_end
            if limit_ind is None or ((side == "in" and trg_end < limit_ind) or
                                     (side == "out" and trg_start >= limit_ind)):
                if side == "in":
                    limit_ind = trg_end
                elif side == "out":
                    limit_ind = trg_start
            f.write("Edge {0}|Template Alignment\n".format(edge_id))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_id, ":", 
                    cons_align[0][0].qry_start, 
                    cons_align[0][0].qry_end, 
                    cons_align[0][0].qry_len))
            f.write("{0:26}\t{1:5}-{2:5} of {3:5}\n".format("Template:",  
                    cons_align[0][0].trg_start, 
                    cons_align[0][0].trg_end, 
                    cons_align[0][0].trg_len))
        f.write("(Largest position considered)\n")
        f.write("{0:26}\t{1}\n\n".format(limit_label, limit_ind))
        #Write out alignment indices for edges vs edges
        edge_pairs = sorted(combinations(edges, 2))
        for edge_one, edge_two in edge_pairs:
            cons_vs_cons = _read_alignment(cons_vs_cons_path.format(it, side, edge_one, 
                                                                    it, side, edge_two), 
                                           consensuses[(it, side, edge_two)], 
                                           min_aln_rate)
            f.write("Edge {0}|Edge {1} Alignment\n".format(edge_one, edge_two))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_one, ":", 
                    cons_vs_cons[0][0].qry_start, 
                    cons_vs_cons[0][0].qry_end, 
                    cons_vs_cons[0][0].qry_len))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_two, ":", 
                    cons_vs_cons[0][0].trg_start, 
                    cons_vs_cons[0][0].trg_end, 
                    cons_vs_cons[0][0].trg_len))
            div_rate = _calculate_divergence(cons_vs_cons[0][0].qry_seq, 
                                             cons_vs_cons[0][0].trg_seq)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
        #Write overall position stats
        confirmed, rejected, pos = _read_confirmed_positions(confirmed_pos_path)
        remaining = len(pos) - (len(confirmed) + len(rejected))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Confirmed Positions:", 
                                              len(confirmed), len(pos), 
                                              len(confirmed)/float(len(pos))))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Rejected Positions:", 
                                              len(rejected), len(pos), 
                                              len(rejected)/float(len(pos))))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Remaining Positions:", 
                                              remaining, len(pos), 
                                              remaining/float(len(pos))))
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
                   template_info, stats_file):
    #Count edge reads
    side_reads = {}
    total_reads = 0
    for side in sorted(repeat_edges[rep]):
        part_list = _read_partitioning_file(partitioning.format(zero_it, side))
        total_reads = len(part_list)
        side_reads[side], tied_reads, unassigned_reads = _get_partitioning_info(
                                                        part_list, 
                                                        repeat_edges[rep][side])
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
    with open(stats_file, 'w') as f:
        f.write("{0:16}\t{1}\n".format("Repeat:", rep))
        f.write("{0:16}\t{1}\n".format("Template Length:", template_info))
        f.write("{0:16}\t{1}\n\n".format("#Total Reads:", total_reads))
        edge_headers = ["Side", "Edge", "#Reads"]
        spaced_edge_header = map("{:5}".format, edge_headers)
        f.write("\t".join(spaced_edge_header))
        f.write("\n")
        for side in sorted(repeat_edges[rep]):
            for edge_id in sorted(repeat_edges[rep][side]):
                edge_values = [side, edge_id, side_reads[side][edge_id]]
                spaced_values = map("{:5}".format, edge_values)
                f.write("\t".join(spaced_values))
                f.write("\n")
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")
        

def update_int_stats(rep, repeat_edges, side_it, cons_align_path, template, 
                     min_aln_rate, confirmed_pos_path, partitioning, 
                     stats_file):
    stats_out = []
    #Add side iters
    for side in sorted(repeat_edges[rep]):
        stats_out.extend([str(side_it[side])])
    #Find mean in, out, and gap lengths
    medians = {s:0 for s in repeat_edges[rep]}
    template_len = 0
    for side in sorted(repeat_edges[rep]):
        trg_limits = []
        for edge_id in sorted(repeat_edges[rep][side]):
            cons_align = _read_alignment(cons_align_path.format(side_it[side], 
                                                                side, 
                                                                edge_id), 
                                         template, 
                                         min_aln_rate)
            template_len = cons_align[0][0].trg_len
            if side == "in":
                trg_limits.append(cons_align[0][0].trg_end)
            elif side == "out":
                trg_limits.append(template_len - cons_align[0][0].trg_start)
        medians[side] = _get_median(trg_limits)
    gap_len = template_len - (medians["in"] + medians["out"])
    stats_out.extend([str(medians["in"]), str(gap_len), str(medians["out"])])
    #Add confirmed and rejected reads
    in_confirmed_path = confirmed_pos_path.format(side_it["in"], "in")
    out_confirmed_path = confirmed_pos_path.format(side_it["out"], "out")
    int_confirmed, int_rejected, pos = _integrate_confirmed_pos(in_confirmed_path,
                                                                out_confirmed_path)
    stats_out.extend([str(len(int_confirmed)), str(len(int_rejected))])
    #Get bridging reads for each pair of in/out edges
    side_headers_dict = {}
    all_headers = set()
    for side in sorted(repeat_edges[rep]):
        side_headers_dict[side] = {}
        part_list = _read_partitioning_file(partitioning.format(side_it[side], side))
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
        stats_out.extend("{0}".format(edge_pair))
        stats_out.extend(str(bridging_reads[edge_pair]))
    spaced_header = map("{:8}".format, stats_out)
    #Write to file
    with open(stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")
        
def finalize_int_stats(rep, repeat_edges, side_it, cons_align_path, template, min_aln_rate, 
                 cons_vs_cons_path, consensuses, confirmed_pos_path, 
                 partitioning, max_iter, edge_below_cov, dup_part, 
                 min_bridge_count, min_bridge_diff, stats_file):
    with open(stats_file, "a") as f:
        f.write("\n\n")
        for side in sorted(repeat_edges[rep]):
            f.write("{0}{1}{2}\t{3}\n".format("Final ", side, " Iter:", 
                                              side_it[side]))
        f.write("\n")
        #Overall confirmed and rejected positions
        template_len #NEED THIS
        in_confirmed_path = confirmed_pos_path.format(side_it["in"], "in")
        out_confirmed_path = confirmed_pos_path.format(side_it["out"], "out")
        int_confirmed, int_rejected, pos = _integrate_confirmed_pos(in_confirmed_path,
                                                                    out_confirmed_path)
        remaining = len(pos) - (len(int_confirmed) + len(int_rejected))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Confirmed Positions:", 
                                              len(int_confirmed), len(pos), 
                                              len(int_confirmed)/float(len(pos))))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Rejected Positions:", 
                                              len(int_rejected), len(pos), 
                                              len(int_rejected)/float(len(pos))))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Remaining Positions:", 
                                              remaining, len(pos), 
                                              remaining/float(len(pos))))
        f.write("\n")
        #Basic stats for confirmed positions
        av_div = 0.0
        if template_len != 0:
            av_div = len(int_confirmed) / float(template_len)
        position_gaps = [0 for x in range(len(int_confirmed)+1)]
        curr_pos = 0
        for i,p in enumerate(pos_list):
            position_gaps[i] = p-curr_pos
            curr_pos = p
        position_gaps[-1] = template_len-curr_pos
        mean_position_gap = np.mean(position_gaps)
        max_position_gap = max(position_gaps)
        f.write("{0:26}\t{1}\n".format("Template Length:", template_len))
        f.write("{0:26}\t{1}\n".format("Confirmed Positions:", len(int_confirmed)))
        f.write("{0:26}\t{1:.4f}\n".format("Confirmed Pos Average Divergence:", 
                                              av_div))
        f.write("{0:26}\t{1:.2f}\n".format("Mean Confirmed Pos Gap:", 
                                              mean_position_gap))
        f.write("{0:26}\t{1}\n".format("Max Confirmed Pos Gap:", 
                                              max_position_gap))
        f.write("\n")
        #Write bridging reads
        side_headers_dict = {}
        all_headers = set()
        for side in sorted(repeat_edges[rep]):
            side_headers_dict[side] = {}
            part_list = _read_partitioning_file(partitioning.format(side_it[side], side))
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
        f.write("\n")
        #Write combos which are sets of bridging reads
        all_combos = _get_combos(side_edges[0], side_edges[1])
        combo_support = [0 for _ in all_combos]
        for i, combo in enumerate(all_combos):
            for edge_pair in combo:
                if edge_pair in bridging_reads:
                    combo_support[i] += bridging_reads[edge_pair]
        for i, combo in enumerate(all_combos):
            f.write("{0} {1}\n".format("Combo", i))
            combo_edges = " + ".join(["|".join(["".join(map(str, x)) for x in y]) for y in combo])
            f.write("{0:26}\t{1}\n".format("Resolution:", combo_edges))
            f.write("{0:26}\t{1}\n\n".format("Support:", combo_support[i]))
        #Bridging conditions 
        bridged = False
        combo_inds = zip(combo_support, range(len(combo_support)))
        sorted_combos = sorted(combo_inds, reverse=True)
        if (sorted_combos[0][0] >= min_bridge_count and 
            sorted_combos[0][0] - sorted_combos[1][0] >= min_bridge_diff):
            bridged = True
        if bridged:
            f.write("BRIDGED\n")
            f.write("Combo {0}\n".format(sorted_combos[0][1]))
            br_ct_str = "{0} (min_bridge_count)".format(min_bridge_count)
            br_diff_str = "{0} + {1} (Combo {2} + min_bridge_diff)".format(
                sorted_combos[1][0], min_bridge_diff, sorted_combos[1][1])
            f.write("Support = {0} > {1} > {2}\n".format(
                sorted_combos[0][0], br_ct_str, br_diff_str))
            f.write("Resolution:\n")
            for edge_pair in all_combos[sorted_combos[0][1]]:
                f.write("{0[0]} {0[1]} {1} {2[0]} {2[1]}\n".format(edge_pair[0], 
                                             u"\u2192", 
                                             edge_pair[1]))
            f.write("\n")
        else:
            f.write("UNBRIDGED\n")
            f.write("Best combo {0}\n".format(sorted_combos[0][1]))
            f.write("{0}\t{1}\n".format("min_bridge_count", min_bridge_count))
            f.write("{0}\t{1}\n".format("min_bridge_diff", min_bridge_diff))
        #Write final in gap out, overlap, repeat sequences, statistics
        #Write where repeats overlap each other for endpoint, divergence rates
        

def _get_combos(in_list, out_list):
    if not in_list or not out_list:
        return []
    else:
        in1 = in_list[0]
        for j in range(len(out_list)):
            combo = (in1, out_list[j])
            for rest in _get_combos(in_list[1:], out_list[:j]+out_list[j+1:]):
                return [combo] + rest


    """finalize_side_stats(edges, it, side, cons_align_path, template, min_aln_length, 
                 cons_vs_cons_path, consensuses, confirmed_pos_path, 
                 partitioning, max_iter, edge_below_cov, dup_part, 
                 stats_file):
    with open(stats_file, "a") as f:
        f.write("\n\n")
        f.write("{0:26}\t{1}\n".format("Last Iter:", it))
        f.write("Iteration terminated because:\n")
        if it == max_iter:
            f.write("Max iter reached\n")
        if edge_below_cov:
            f.write("Edge coverage fell below min_edge_cov\n")
        if dup_part:
            f.write("Partitioning was identical to a previous iteration\n")
        f.write("\n")
        #Write out alignment indices for edges vs template
        limit_ind = None
        limit_label = ""
        if side == "in":
            limit_label = "Min Template End"
        elif side == "out":
            limit_label = "Max Template Start"
        for edge_id in sorted(edges):
            cons_align = _read_alignment(cons_align_path.format(it, side, edge_id), 
                                     template, 
                                     min_aln_length)
            trg_start = cons_align[0][0].trg_start
            trg_end = cons_align[0][0].trg_end
            if limit_ind is None or ((side == "in" and trg_end < limit_ind) or
                                     (side == "out" and trg_start >= limit_ind)):
                if side == "in":
                    limit_ind = trg_end
                elif side == "out":
                    limit_ind = trg_start
            f.write("Edge {0}|Template Alignment\n".format(edge_id))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_id, ":", 
                    cons_align[0][0].qry_start, 
                    cons_align[0][0].qry_end, 
                    cons_align[0][0].qry_len))
            f.write("{0:26}\t{1:5}-{2:5} of {3:5}\n".format("Template:",  
                    cons_align[0][0].trg_start, 
                    cons_align[0][0].trg_end, 
                    cons_align[0][0].trg_len))
        f.write("(Largest position considered)\n")
        f.write("{0:26}\t{1}\n\n".format(limit_label, limit_ind))
        #Write out alignment indices for edges vs edges
        edge_pairs = sorted(combinations(edges, 2))
        for edge_one, edge_two in edge_pairs:
            cons_vs_cons = _read_alignment(cons_vs_cons_path.format(it, side, edge_one, 
                                                                    it, side, edge_two), 
                                           consensuses[(it, side, edge_two)], 
                                           min_aln_length)
            f.write("Edge {0}|Edge {1} Alignment\n".format(edge_one, edge_two))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_one, ":", 
                    cons_vs_cons[0][0].qry_start, 
                    cons_vs_cons[0][0].qry_end, 
                    cons_vs_cons[0][0].qry_len))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_two, ":", 
                    cons_vs_cons[0][0].trg_start, 
                    cons_vs_cons[0][0].trg_end, 
                    cons_vs_cons[0][0].trg_len))
            div_rate = _calculate_divergence(cons_vs_cons[0][0].qry_seq, 
                                             cons_vs_cons[0][0].trg_seq)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
        #Write overall position stats
        confirmed, rejected, pos = _read_confirmed_positions(confirmed_pos_path)
        remaining = len(pos) - (len(confirmed) + len(rejected))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Confirmed Positions:", 
                                              len(confirmed), len(pos), 
                                              len(confirmed)/float(len(pos))))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Rejected Positions:", 
                                              len(rejected), len(pos), 
                                              len(rejected)/float(len(pos))))
        f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format("Total Remaining Positions:", 
                                              remaining, len(pos), 
                                              remaining/float(len(pos))))
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
        f.write("\n")    """
    
        
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
            raise Exception("No alignment conditions fit, {0} {1}".format(q, t))
    if curr_del != 0:
        del_count += 1
        curr_del = 0
    if curr_ins != 0:
        ins_count += 1
        curr_ins = 0
    
    indel_sim_rate = match_count / float(match_count + 
                                         mis_count + 
                                         del_count + 
                                         ins_count)
    return 1 - indel_sim_rate

def _get_median(lst):
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
    
    in_conf_set = set(in_conf)
    in_rej_set = set(in_rej)
    out_conf_set = set(out_conf)
    out_rej_set = set(out_rej)
    integrated_confirmed = []
    integrated_rejected = []
    
    for pos in sorted(in_pos):
        if pos in in_conf_set or pos in out_conf_set:
            integrated_confirmed.append(pos)
        elif pos in in_rej_set or pos in out_rej_set:
            integrated_rejected.append(pos)
    return integrated_confirmed, integrated_rejected, in_pos


def main_two():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    
    parser.add_argument('all_reads_file',action='store')
    parser.add_argument('repeats_file',action='store')
    parser.add_argument('graph_edges_file',action='store')
    parser.add_argument('output_stem',action='store')
    
    args = parser.parse_args()
    
    repeats_dict = readRepeatsFile(args.repeats_file)
    start_headers,reads = readAllLongFASTA(args.all_reads_file)
    headers = [h.split(' ')[0] for h in start_headers]
    headers_dict = {h:(i,r) for (i,(h,r)) in enumerate(zip(headers,reads))}
    graph_headers,graph_seqs = readAllLongFASTA(args.graph_edges_file)
    graph_dict = getGraphDict(graph_headers,graph_seqs)
    
    out_all_stem = args.output_stem+'.R%d.all.reads.fasta'
    out_start_stem = args.output_stem+'.R%d.C%d.start.reads.fasta'
    out_end_stem = args.output_stem+'.R%d.C%d.end.reads.fasta'
    out_template_stem = args.output_stem+'.R%d.template.fasta'
    
    for r in repeats_dict:
        count,all_reads_list,inputs_dict,outputs_dict = repeats_dict[r]
        writeReadsFromHeaderDict(out_all_stem%(r),headers_dict,all_reads_list)
        for i in inputs_dict:
            input_list = inputs_dict[i]
            writeReadsFromHeaderDict(out_start_stem%(r,i),headers_dict,input_list)
        for i in outputs_dict:
            output_list = outputs_dict[i]
            writeReadsFromHeaderDict(out_end_stem%(r,i),headers_dict,output_list)
        writeTemplateFromGraphDict(out_template_stem%(r),r,graph_dict)
    
    read_translation_file = args.output_stem+'.headers.to.read.num.fasta'
    with open(read_translation_file,'w') as rtf:
        for i,h in enumerate(start_headers):
            rtf.write('%s\n%d\n' % (h,i))
        
    print 'Total Time:\t%.3f seconds' % (time.time()-start_time)

def readRepeatsFile(repeats_file):
    repeats_dict = {}
    curr_repeat = -1
    count = 0
    all_reads_list = []
    all_bool = False
    inputs_dict = {}
    curr_input = -1
    curr_input_list = []
    input_bool = False
    outputs_dict = {}
    curr_output = -1
    curr_output_list = []
    output_bool = False
    with open(repeats_file,'r') as rf:
        for i,line in enumerate(rf):
            line = line.strip()
            if line:
                if line[0] == '#':    
                    all_bool = False
                    if input_bool:
                        inputs_dict[curr_input] = curr_input_list[:]
                        curr_input_list = []
                        curr_input = -1
                        input_bool = False
                    if output_bool:
                        outputs_dict[curr_output] = curr_output_list[:]
                        curr_output_list = []
                        curr_output = -1
                        output_bool = False
                    
                    parts = line.split('\t')
                    header = parts[0]
                    count = int(parts[1])
                    head_parts = header.split(' ')
                    if head_parts[0] == '#Repeat':
                        if curr_repeat != -1:
                            repeats_dict[curr_repeat] = (count,all_reads_list,inputs_dict,outputs_dict)
                            all_reads_list = []
                            inputs_dict = {}
                            outputs_dict = {}
                        
                        curr_repeat = int(head_parts[1])
                    elif head_parts[0] == '#All':
                        all_bool = True
                    elif head_parts[0] == '#Input':                        
                        curr_input = int(head_parts[1])
                        input_bool = True
                    elif head_parts[0] == '#Output':
                        curr_output = int(head_parts[1])
                        output_bool = True
                else:
                    if all_bool:
                        all_reads_list.append(line)
                    elif input_bool:
                        curr_input_list.append(line)
                    elif output_bool:
                        curr_output_list.append(line)
        if output_bool:
            outputs_dict[curr_output] = curr_output_list[:]
        if input_bool:
            inputs_dict[curr_input] = curr_input_list[:]
        if curr_repeat != -1:
            repeats_dict[curr_repeat] = (count,all_reads_list,inputs_dict,outputs_dict)
    return repeats_dict

def getGraphDict(headers,seqs):
    '''Assumes that headers are in the form >linear_###. Outputs a dict
    of ###:seq for each seq'''
    graph_dict = {}
    for h,s in zip(headers,seqs):
        num = int(h.split('_')[1])
        graph_dict[num] = s
    return graph_dict
                
def writeReadsFromHeaderDict(out_file,headers_dict,headers):
    '''Headers will be in the form -h or +h, header_dict is in the form >h,
    the read that will be written will be the revComp if the header is -h'''
    with open(out_file,'w') as of:
        rev = False
        for h in headers:
            if h[0] == '-':
                rev = True
            elif h[0] == '+':
                rev = False
            else:
                print 'Header not in the right form:\t%s' % h
            
            fixed_h = '>'+h[1:]
            
            i,r = headers_dict[fixed_h]
            
            if not rev:
                of.write('%s%s|%d|\n%s\n' % ('>',h,i,r))
            else:
                of.write('%s%s|%d|\n%s\n' % ('>',h,i,rev_comp(r)))

def writeTemplateFromGraphDict(out_file,num,graph_dict):
    with open(out_file,'w') as of:
        rev = False
        if num < 0:
            rev = True
            num = -num
        if not rev:
            of.write('>Template_%d\n%s\n' % (num,graph_dict[num]))
        else:
            of.write('>Template_%d\n%s\n' % (-num,rev_comp(graph_dict[num])))

#readAllLongFASTA accounts for reads spanning multiple lines and
#also '' for reads
def readAllLongFASTA(file_name):
    '''Given a FASTA file of reads where the header starts with a '>' and the
    corresponding read can span multiple succeeding lines, return a tuple of
    headers,seqs - read all lines, even if one is blank'''
    seq_file = open(file_name,'r')
    headers = []
    seqs = []
    curr_seq = ''
    for i,line in enumerate(seq_file):
        line = line.strip()
        if line and line[0] == '>':
            headers.append(line)
            if i > 0:
                seqs.append(curr_seq.upper())
            curr_seq = ''
        else:
            curr_seq += line
    seqs.append(curr_seq.upper())
    
    seq_file.close()
    return headers,seqs

def rev_comp(seq):
    '''Given a DNA sequence, returns the reverse complement of it'''
    trans = {'A':'T','C':'G','G':'C','T':'A','-':'-'}
    out = ''
    for s in seq[::-1]:
        out += trans[s]
    return out

#if __name__ == '__main__':
#    main()

class Profile:
    __slots__ = ("insertions", "matches", "nucl")

    def __init__(self):
        self.insertions = defaultdict(str)
        self.matches = defaultdict(int)
        self.nucl = "-"

def _thread_worker(aln_reader, contigs_info, platform, results_queue,
                   error_queue):
    try:
        aln_reader.init_reading()

        while not aln_reader.is_eof():
            ctg_id, ctg_aln = aln_reader.get_chunk()
            if ctg_id is None:
                break

            profile, aln_errors = _contig_profile(ctg_aln, platform,
                                                  contigs_info[ctg_id].length)
            #sequence = _flatten_profile(profile)
            results_queue.put((ctg_id, profile, aln_errors))

    except Exception as e:
        error_queue.put(e)

def _contig_profile(alignment, platform, genome_len):
    """
    Computes alignment profile
    """
    max_aln_err = config.vals["err_modes"][platform]["max_aln_error"]
    aln_errors = []
    profile = [Profile() for _ in xrange(genome_len)]
    for aln in alignment:
        #if aln.err_rate > max_aln_err: continue
        aln_errors.append(aln.err_rate)

        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)
        #qry_seq = aln.qry_seq
        #trg_seq = aln.trg_seq

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in izip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            prof_elem = profile[trg_pos]
            if trg_nuc == "-":
                prof_elem.insertions[aln.qry_id] += qry_nuc
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    return profile, aln_errors

def _call_positions(profile, sub_thresh, del_thresh):
    positions = {}
    
    for index,elem in enumerate(profile):
        pos_matches = elem.matches
        pos_insertions = elem.insertions
        pos_nucl = elem.nucl
        
        coverage = sum(pos_matches.values())
        
        mismatches = {key:pos_matches[key] for key in pos_matches
                                           if (key != pos_nucl and key != "-")}
    
        max_val = 0
        if len(mismatches):
            max_mismatch = max(mismatches, key=mismatches.get)
            max_val = mismatches[max_mismatch]
                
        del_key = "-"
        del_val = 0
        if del_key in pos_matches:
            del_val = pos_matches[del_key]            
        
        if coverage:
            if max_val / float(coverage) >= sub_thresh:
                positions[index] = (pos_nucl,max_mismatch)
            elif del_val / float(coverage) >= del_thresh:
                positions[index] = (pos_nucl,del_key)
        
    return positions

def find_divergence(alignment_path, contigs_path, contigs_info, min_aln_length,
                    platform, num_proc, sub_thresh, del_thresh):
    """
    Main function: takes in an alignment and finds the divergent positions
    """
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_fasta_dict(contigs_path),
                                       min_aln_length)
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    error_queue = manager.Queue()

    #making sure the main process catches SIGINT
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    threads = []
    for _ in xrange(num_proc):
        threads.append(multiprocessing.Process(target=_thread_worker,
                                               args=(aln_reader, contigs_info,
                                                     platform, results_queue,
                                                     error_queue)))
    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()

    if not error_queue.empty():
        raise error_queue.get()

    out_positions = {}
    total_positions = 0
    total_aln_errors = []
    #with open("debug_positions_check.txt",'w') as dpc:
    while not results_queue.empty():
        ctg_id, ctg_profile, aln_errors = results_queue.get()
        #dpc.write('Index\tMatches\tInsertions\tBase\n')
        #for index,elem in enumerate(ctg_profile):
        #    dpc.write("{0}\t{1}\t{2}\t{3}\n".format(index, elem.matches, elem.insertions, elem.nucl))
            
        ctg_positions = _call_positions(ctg_profile, sub_thresh, del_thresh)
        total_positions += len(ctg_positions)
        total_aln_errors.extend(aln_errors)
        if len(ctg_positions) > 0:
            out_positions[ctg_id] = ctg_positions
    
    mean_aln_error = float(sum(total_aln_errors)) / (len(total_aln_errors) + 1)
    logger.debug("Total positions: {0}".format(total_positions))
    logger.debug("Alignment error rate: {0}".format(mean_aln_error))

    return out_positions

def write_positions(positions, positions_file):
    """
    Writes dictionary with positions to file
    """
    with open(positions_file, "w") as f:
        for header in sorted(positions):
            ctg_pos = positions[header]
            f.write(">{0}_num:{1}\n".format(header, len(ctg_pos)))
            
            for elem in sorted(ctg_pos):
                var1, var2 = ctg_pos[elem]
                f.write("{0}:{1}|{2},".format(elem, var1, var2))
            f.write("\n")

def write_sub_n_template(contigs_path, positions, template_file):
    contigs = fp.read_fasta_dict(contigs_path)
    for ctg_header, pos_header in izip(contigs, positions):
        ctg_seq = contigs[ctg_header]
        ctg_pos = positions[pos_header]
        sub_n_seq = ''.join([s if i not in ctg_pos else 'N' for i,s in enumerate(ctg_seq)])
    
    with open(template_file, "w") as f:
        for ctg_header, pos_header in izip(sorted(contigs), sorted(positions)):
            ctg_pos = positions[pos_header]
            f.write(">{0}_sub_n\n{1}\n".format(ctg_header,sub_n_seq))
            

class PositionIOError(Exception):
    pass

def read_positions(positions_file):
    """
    Reads positions file into dictionary
    """
    header = None
    pos_parts = None
    pos_id = None
    pos_bases = None
    ctg_positions = {}
    all_positions = {}
    try:
        with open(positions_file, "r") as f:
            for line_id, line in enumerate(f):
                line = line.strip()
                if line.startswith(">"):
                    header = line[1:].split(" ")[0]
                    all_positions[header] = {}
                    ctg_positions = {}
                elif line:
                    pos_parts = line.split(",")
                    for p in pos_parts:
                        if p.strip():
                            pos_id = p.split(":")
                            pos_bases = pos_id[1].split("|")
                            ctg_positions[int(pos_id[0])] = pos_bases
                    all_positions[header] = copy.copy(ctg_positions)
                print len(all_positions)
    except IOError as e:
        raise PositionIOError(e)
    print all_positions.keys()
    return all_positions

def _pos_thread_worker(aln_reader, contigs_info, platform, results_queue, 
                       error_queue, positions):
    try:
        aln_reader.init_reading()

        while not aln_reader.is_eof():
            ctg_id, ctg_aln = aln_reader.get_chunk()
            if ctg_id is None:
                break

            qry_pos_bases, aln_errors, total_qry_len, trg_len, num_qry = _contig_pos_bases(ctg_aln, platform, 
                                                          positions[ctg_id],
                                                  contigs_info[ctg_id].length)
            results_queue.put((ctg_id, qry_pos_bases, aln_errors, total_qry_len, trg_len, num_qry))

    except Exception as e:
        error_queue.put(e)

def _contig_pos_bases(alignment, platform, ctg_positions, genome_len):
    """
    Computes alignment profile
    """
    max_aln_err = config.vals["err_modes"][platform]["max_aln_error"]
    aln_errors = []
    total_qry_len = 0
    trg_len = 0
    #pos_bases is a dictionary from position -> read_index -> read_base
    pos_bases = {pos:{} for pos in ctg_positions}
    for i,aln in enumerate(alignment):
        #if aln.err_rate > max_aln_err: continue
        aln_errors.append(aln.err_rate)
        total_qry_len += aln.qry_len
        trg_len = aln.trg_len

        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)
        #qry_seq = aln.qry_seq
        #trg_seq = aln.trg_seq

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in izip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            if trg_pos in ctg_positions:   
                pos_bases[trg_pos][i] = qry_nuc
            
            trg_pos += 1

    return pos_bases, aln_errors, total_qry_len, trg_len, len(alignment)

def _filter_with_doubles(positions, qry_pos_bases, thresh, ctg_id, mean_qry_len):
    max_dist = mean_qry_len/2
    pairs = _all_pairs_within_dist(positions, max_dist)
    print ctg_id
    print "Num positions", len(positions)
    print "Num pairs", len(pairs)
    #for pos in sorted(qry_pos_bases):
    #    print len(qry_pos_bases[pos])
    #    for read_ind in sorted(qry_pos_bases[pos]):
    #        print read_ind
        
    #empty_pairs is a list of indices corresponding to the position pairs
    #in pair_scores and pairs
    pairs, empty_pairs = _score_pairings(pairs, qry_pos_bases)
    
    graph_file = "adjacency_graph_{0}.txt".format(ctg_id)
    _write_graph(pairs, graph_file, positions, qry_pos_bases)
    
    
    
    filtered_pairs, empty_pairs = _filter_pairs(pairs, thresh, 
                                                empty_pairs)
    """
    inconsistent, nondivergent = _check_consistency(filtered_pairs)
    
    _print_failed(ctg_id, pairs, pair_scores, filtered_pairs, empty_pairs, inconsistent, 
                  nondivergent)
    
    new_positions = _filter_failed(positions, pairs, empty_pairs, inconsistent, 
                          nondivergent)
    
    count = 0
    while (inconsistent or nondivergent) and count < 10:
        print count
        print "inconsistent", len(inconsistent), inconsistent
        print "nondivergent", len(nondivergent), nondivergent
        print "all positions", len(new_positions)
        pairs = _make_position_sets(sorted(new_positions.keys()), within)
        pair_scores, empty_pairs = _score_pairings(pairs, qry_pos_bases)
        #print pair_scores
        filtered_pairs, empty_pairs = _filter_pairs(pair_scores, thresh, 
                                                    empty_pairs)
        
        inconsistent, nondivergent = _check_consistency(filtered_pairs)
        
        _print_failed(ctg_id, pairs, pair_scores, filtered_pairs, empty_pairs, inconsistent, 
                      nondivergent)
        
        new_positions = _filter_failed(new_positions, pairs, empty_pairs, inconsistent, 
                              nondivergent)
        count += 1
    return new_positions"""
    return []
    

def _make_position_sets(pos_list, within):
    """Creates a list of all adjacent positions of within size"""
    pairs = []
    for i in range(len(pos_list)-within):
        curr_pos = pos_list[i]
        
        for j in range(within):
            pair_pos = pos_list[i+j+1]
            pairs.append((curr_pos,pair_pos))
    return pairs

def _all_pairs_within_dist(pos, dist):
    """Creates a dictionary of pairs of positions that are less than dist apart"""
    pairs = {}
    for p1 in sorted(pos):
        pairs[p1] = {}
        for p2 in sorted(pos):
            if abs(p2 - p1) < dist:
                p1_base1, p1_base2 = pos[p1]
                p2_base1, p2_base2 = pos[p2]
                pairs[p1][p2] = ({(p1_base1, p2_base1) : 0, 
                                  (p1_base2, p2_base2) : 0},
                                 {(p1_base1, p2_base2) : 0,
                                  (p1_base2, p2_base1) : 0})
    return pairs

def _score_pairings(pairs, qry_pos_bases):
    empty_pairs = set()
    for p1 in sorted(pairs):
        for p2 in sorted(pairs[p1]):
            if p1 != p2:
                for read_ind in qry_pos_bases[p1]:
                    if read_ind in qry_pos_bases[p2]:
                        b1 = qry_pos_bases[p1][read_ind]
                        b2 = qry_pos_bases[p2][read_ind]
                        if (b1, b2) in pairs[p1][p2][0]:
                            pairs[p1][p2][0][(b1, b2)] += 1
                        elif (b1, b2) in pairs[p1][p2][1]:
                            pairs[p1][p2][1][(b1, b2)] += 1
            if sum(pairs[p1][p2][0].values() + pairs[p1][p2][1].values()) == 0:
                empty_pairs.add((p1, p2))
    return pairs, empty_pairs

def _write_graph(pairs, graph_file, pos, qry_pos_bases):
    with open(graph_file, "w") as f:
        f.write("P1\\P2\t")
        for p1 in sorted(pairs):
            f.write("{0} {1}|{2} {3} {4}\t".format(p1, pos[p1][0], pos[p1][1], 
                                               len(qry_pos_bases[p1]), qry_pos_bases[p1]))
        f.write("\n")
        
        for p1 in sorted(pairs):
            f.write("{0} {1}|{2} {3} {4}\t".format(p1, pos[p1][0], pos[p1][1], 
                                               len(qry_pos_bases[p1]), qry_pos_bases[p1]))
            for p2 in sorted(pairs):
                if p2 in pairs[p1]:
                    #f.write("{0}|{1}  ".format(p1, p2))
                    for b1, b2 in sorted(pairs[p1][p2][0]):
                        count = pairs[p1][p2][0][(b1, b2)]
                        if count > 0:
                            f.write("{0},{1}:{2} ".format(b1, b2, count))
                    f.write("vs ")
                    for b1, b2 in sorted(pairs[p1][p2][1]):
                        count = pairs[p1][p2][1][(b1, b2)]
                        if count > 0:
                            f.write("{0},{1}:{2} ".format(b1, b2, count))
                    f.write("\t")
                else:
                    f.write("\t")
            f.write("\n")
                

def _filter_pairs(pairs, thresh, empty_pairs):
    filtered_pairs = {}
    for p1 in pairs:
        for p2 in pairs[p1]:
            for b1, b2 in pairs[p1][p2]:
                pair_scores = pairs[p1][p2][(b1, b2)]
                for read_ind in i, pos_pair in enumerate(pair_scores):
                    total = sum(pos_pair.values())
                    for base_pair in pos_pair:
                        if pos_pair[base_pair] / float(total) >= thresh:
                            filtered_pairs[i][base_pair] = pos_pair[base_pair]
                    if not filtered_pairs[i]:
                        empty_pairs.add(i)
    return filtered_pairs, empty_pairs

def _check_consistency(pair_scores):
    min_pairs = 2
    inconsistent = []
    nondivergent = []
    #prev_bases = None means that any base will work, otherwise it will be a
    #set of strings (base/indel)
    prev_bases = None
    for i in range(len(pair_scores)):
        curr_pair = pair_scores[i]
        curr_cands = []
        if prev_bases is None:
            for pos_pair in curr_pair:
                curr_cands.append(pos_pair)
        else:
            for pos_pair in curr_pair:
                if pos_pair[0] in prev_bases:
                    curr_cands.append(pos_pair)
        
        #Inconsistent if not at least min_pairs are found
        if len(curr_cands) < min_pairs:
            inconsistent.append(i)
            prev_bases = None
        else:
            prev_bases = set()
            for first, second in curr_cands:
                prev_bases.add(second)
        
        if prev_bases is not None and len(prev_bases) < min_pairs:
            nondivergent.append(i)
            
    return inconsistent, nondivergent
        
def _print_failed(ctg_id, pairs, pair_scores, filtered_scores, empty_pairs, inconsistent,
                  nondivergent):
    with open("debug_position_pairs_output_{0}.txt".format(ctg_id),'w') as f:
        f.write("Empty Pairs:\nIndex\tPos_1 , Pos_2\n")
        for i in sorted(empty_pairs):
            pos_1, pos_2 = pairs[i]
            f.write("{0}\t{1} , {2}\n".format(i, pos_1, pos_2))
        
        f.write("\nInconsistent Pairs:\nIndex\tPos_1 , Pos_2\tCount\tTotal_Score\n")
        for i in inconsistent:
            pos_1, pos_2 = pairs[i]
            bases = filtered_scores[i]
            total_score = sum(bases.values())
            f.write("{0}\t{1} , {2}\t{3}\t{4}\n".format(i, pos_1, pos_2, 
                    len(bases), total_score))
        
        f.write("\nNondivergent Pairs:\nIndex\tPos_1 , Pos_2\tCount\tTotal_Score\n")
        for i in nondivergent:
            pos_1, pos_2 = pairs[i]
            bases = filtered_scores[i]
            total_score = sum(bases.values())
            f.write("{0}\t{1} , {2}\t{3}\t{4}\n".format(i, pos_1, pos_2, 
                    len(bases), total_score))
        
        f.write("\nAll Filtered Pair Scores:\n")
        f.write("Index\tPos_1 , Pos_2\tBase_1 , Base_2\tScore\n")
        for i, bases in enumerate(filtered_scores):
            pos_1, pos_2 = pairs[i]
            for base_pair in bases:
                base_1, base_2 = base_pair
                score = bases[base_pair]
                f.write("{0}\t{1} , {2}\t{3} , {4}\t{5}\n".format(i, pos_1, 
                        pos_2, base_1, base_2, score))
            f.write("\n")
        
        f.write("\nAll Pair Scores:\n")
        f.write("Index\tPos_1 , Pos_2\tBase_1 , Base_2\tScore\n")
        for i, bases in enumerate(pair_scores):
            pos_1, pos_2 = pairs[i]
            for base_pair in bases:
                base_1, base_2 = base_pair
                score = bases[base_pair]
                f.write("{0}\t{1} , {2}\t{3} , {4}\t{5}\n".format(i, pos_1, 
                        pos_2, base_1, base_2, score))
            f.write("\n")

def _filter_failed(positions, pairs, empty_pairs, inconsistent, nondivergent):
    new_positions = copy.deepcopy(positions)
    for i in nondivergent:
        pos_1, pos_2 = pairs[i]
        del new_positions[pos_2]
    
    #for i in inconsistent:
    #    pos_1, pos_2 = pairs[i]
    #    if pos_1 in new_positions:
    #        del new_positions[pos_1]
    
    return new_positions
        

def _match_positions_to_contigs(positions, contigs_info):
    ctg_key_to_pos_key = {}
    for pos_header in positions:
        ctg_key = "_".join(pos_header.split("_")[:2])#+["sub","n"])
        if ctg_key in contigs_info:
            ctg_key_to_pos_key[ctg_key] = pos_header
        #Raise error here if the pos_header doesn't match
            
    matched_positions = {}
    for ctg_id in contigs_info:
        if ctg_id in ctg_key_to_pos_key:
            matched_positions[ctg_id] = positions[ctg_key_to_pos_key[ctg_id]]
        else:
            #Raise potential error here if ctg_id is not found
            matched_positions[ctg_id] = {}
    return matched_positions
    
def confirm_positions(alignment_path, contigs_path, contigs_info, min_aln_length,
                    platform, num_proc, positions_file, confirm_thresh):
    """
    Uses frequency of base occurrences at pairs of positions in reads
    to confirm or reject positions. Can be used iteratively.
    """
    orig_positions = read_positions(positions_file)
    positions = _match_positions_to_contigs(orig_positions,contigs_info)
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_fasta_dict(contigs_path),
                                       min_aln_length)
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    error_queue = manager.Queue()

    #making sure the main process catches SIGINT
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    threads = []
    for _ in xrange(num_proc):
        threads.append(multiprocessing.Process(target=_pos_thread_worker,
                                               args=(aln_reader, contigs_info, 
                                                     platform, results_queue, 
                                                     error_queue, positions)))
    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()

    if not error_queue.empty():
        raise error_queue.get()

    out_positions = {}
    total_positions = 0
    total_aln_errors = []
    #with open("debug_positions_check.txt",'w') as dpc:
    while not results_queue.empty():
        ctg_id, ctg_pos_bases, aln_errors, total_qry_len, trg_len, num_qry = results_queue.get()
        #dpc.write('Index\tMatches\tInsertions\tBase\n')
        #for index,elem in enumerate(ctg_profile):
        #    dpc.write("{0}\t{1}\t{2}\t{3}\n".format(index, elem.matches, elem.insertions, elem.nucl))
        
        cov = 0.0
        if trg_len > 0:
            cov = total_qry_len/float(trg_len)
        mean_qry_len = 0
        if num_qry > 0:
            mean_qry_len = total_qry_len/float(num_qry)
        print ctg_id, total_qry_len, trg_len, cov
        print "Mean read length", mean_qry_len
        ctg_positions = _filter_with_doubles(positions[ctg_id], ctg_pos_bases, 
                                      confirm_thresh, ctg_id, mean_qry_len)
        total_positions += len(ctg_positions)
        total_aln_errors.extend(aln_errors)
        if len(ctg_positions) > 0:
            out_positions[ctg_id] = ctg_positions
    
    mean_aln_error = float(sum(total_aln_errors)) / (len(total_aln_errors) + 1)
    logger.debug("Total positions: {0}".format(total_positions))
    logger.debug("Alignment error rate: {0}".format(mean_aln_error))

    return out_positions