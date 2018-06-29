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

def find_coverage(frequency_file):
    header, freqs = div.read_frequency_path(frequency_file)
    cov_ind = header.index("Cov")
    all_covs = [f[cov_ind] for f in freqs]
    print min(all_covs), np.mean(all_covs), max(all_covs)
    return np.mean(all_covs)

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
    max_end_edge = max([_get_aln_end(alns[e].trg_start, alns[e].trg_seq) for e in alns])
    #end indices for conservatively defining confirmed positions
    min_end_edge = min([_get_aln_end(alns[e].trg_start, alns[e].trg_seq) for e in alns])
    max_start_edge = max([alns[e].trg_start for e in alns])
    
    for trg_ind in range(min_start_edge, max_end_edge):
        for edge_id in alns:
            aln = alns[edge_id]
            if aln.trg_start == trg_ind:
                aln.curr_aln_ind = 0
                aln.curr_qry_ind = aln.qry_start
                aln.in_alignment = True
            
            if aln.trg_start > trg_ind or _get_aln_end(aln.trg_start, aln.trg_seq) <= trg_ind:
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
                if pos >= aln.trg_start and pos < _get_aln_end(aln.trg_start, aln.trg_seq):
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
        f.write("{0:25}\t{1}\n".format("Repeat:", rep))
        f.write("{0:25}\t'{1}'\n".format("Side:", side))
        f.write("{0:25}\t".format("Edges:"))
        f.write(", ".join(map(str,sorted(repeat_edges[rep][side]))))
        f.write("\n")
        f.write("{0:25}\t{1}\n\n".format("Template Length:", template_len))
        f.write("Initial Option Values\n")
        f.write("{0:25}\t{1}\n".format("min_overlap:", args.min_overlap))
        f.write("{0:25}\t{1}\n".format("sub_thresh:", args.sub_thresh))
        f.write("{0:25}\t{1}\n".format("del_thresh:", args.del_thresh))
        f.write("{0:25}\t{1}\n".format("ins_thresh:", args.ins_thresh))
        f.write("{0:25}\t{1}\n".format("extend_len:", args.extend_len))
        f.write("{0:25}\t{1}\n".format("buffer_count:", buffer_count))
        f.write("{0:25}\t{1}\n".format("max_iter:", max_iter))
        f.write("{0:25}\t{1}\n".format("min_edge_cov:", min_edge_cov))
        f.write("\n")
        f.write("The following numbers are calculated based on moving ")
        f.write("into the repeat from the '{0}' direction\n\n".format(side))
        f.write("{0:25}\t{1}\n".format("Divergent Positions:", len(pos)))
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
            rep_len = _get_aln_end(cons_align[0][0].qry_start, cons_align[0][0].qry_seq)
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
        f.write("{0:26}\t{1}\n\n".format("Final Iter:", it))
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
            trg_end = _get_aln_end(cons_align[0][0].trg_start, cons_align[0][0].trg_seq)
            if limit_ind is None or ((side == "in" and trg_end < limit_ind) or
                                     (side == "out" and trg_start >= limit_ind)):
                if side == "in":
                    limit_ind = trg_end
                elif side == "out":
                    limit_ind = trg_start
            f.write("Edge {0}|Template Alignment\n".format(edge_id))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_id, ":", 
                    cons_align[0][0].qry_start, 
                    _get_aln_end(cons_align[0][0].qry_start, cons_align[0][0].qry_seq), 
                    cons_align[0][0].qry_len))
            f.write("{0:26}\t{1:5}-{2:5} of {3:5}\n".format("Template:",  
                    cons_align[0][0].trg_start, 
                    _get_aln_end(cons_align[0][0].trg_start, cons_align[0][0].trg_seq), 
                    cons_align[0][0].trg_len))
        f.write("\n")
        f.write("{0:26}\t{1}\n".format(limit_label, limit_ind))
        f.write("(End of positions considered)\n\n")
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
                    _get_aln_end(cons_align[0][0].qry_start, cons_align[0][0].qry_seq), 
                    cons_vs_cons[0][0].qry_len))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format("Edge ", edge_two, ":", 
                    cons_vs_cons[0][0].trg_start, 
                    _get_aln_end(cons_align[0][0].trg_start, cons_align[0][0].trg_seq), 
                    cons_vs_cons[0][0].trg_len))
            div_rate = _calculate_divergence(cons_vs_cons[0][0].qry_seq, 
                                             cons_vs_cons[0][0].trg_seq)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
        #Write overall position stats
        confirmed, rejected, pos = _read_confirmed_positions(confirmed_pos_path)
        remaining = len(pos) - (len(confirmed) + len(rejected))
        if side == "in":
            f.write("{0:26}\t{1}\n".format("Largest Confirmed Position:", 
                                           max(confirmed)))
        elif side == "out":
            f.write("{0:26}\t{1}\n".format("Smallest Confirmed Position:", 
                                           min(confirmed)))
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
                   all_reads_file, template_len, cov, int_stats_file):
    #Count edge reads
    side_reads = {}
    total_reads = 0
    all_side_reads = 0
    internal_reads = 0
    for side in sorted(repeat_edges[rep]):
        part_list = _read_partitioning_file(partitioning.format(zero_it, side))
        total_reads = len(part_list)
        side_reads[side], tied_reads, unassigned_reads = _get_partitioning_info(
                                                        part_list, 
                                                        repeat_edges[rep][side])
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
        edge_headers = ["Side", "Edge", "# Reads"]
        spaced_edge_header = map("{:5}".format, edge_headers)
        f.write("\t".join(spaced_edge_header))
        f.write("\n")
        for side in sorted(repeat_edges[rep]):
            for edge_id in sorted(repeat_edges[rep][side]):
                edge_values = [side, edge_id, side_reads[side][edge_id]]
                spaced_values = map("{:6}".format, edge_values)
                f.write("\t".join(spaced_values))
                f.write("\n")
        f.write("{0:16}\t{1}\n".format("Internal", internal_reads))
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")
        

def update_int_stats(rep, repeat_edges, side_it, cons_align_path, template, 
                     template_len, min_aln_rate, confirmed_pos_path, 
                     int_confirmed_path, partitioning, int_stats_file):
    stats_out = []
    #Add side iters
    for side in sorted(repeat_edges[rep]):
        stats_out.extend([str(side_it[side])])
    #Find median in, out, and gap lengths
    medians = {s:0 for s in repeat_edges[rep]}
    for side in sorted(repeat_edges[rep]):
        trg_limits = []
        for edge_id in sorted(repeat_edges[rep][side]):
            cons_align = _read_alignment(cons_align_path.format(side_it[side], 
                                                                side, 
                                                                edge_id), 
                                         template, 
                                         min_aln_rate)
            if side == "in":
                trg_limits.append(_get_aln_end(cons_align[0][0].trg_start, cons_align[0][0].trg_seq))
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
    _write_confirmed_positions(int_confirmed, int_rejected, pos, 
                               int_confirmed_path.format(side_it["in"], 
                                                         side_it["out"]))
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
        #stats_out.extend(["{0}".format(edge_pair)])
        stats_out.extend([str(bridging_reads[edge_pair])])
    spaced_header = map("{:8}".format, stats_out)
    #Write to file
    with open(int_stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")
        
def finalize_int_stats(rep, repeat_edges, side_it, cons_align_path, template, 
                       template_len, min_aln_rate, cons_vs_cons_path, 
                       consensuses, int_confirmed_path, partitioning, 
                       min_bridge_count, min_bridge_diff, int_stats_file, 
                       resolved_seq_file):
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
        for i,p in enumerate(int_confirmed):
            position_gaps[i] = p-curr_pos
            curr_pos = p
        position_gaps[-1] = template_len-curr_pos
        mean_position_gap = np.mean(position_gaps)
        max_position_gap = max(position_gaps)
        f.write("{0:26}\t{1}\n".format("Template Length:", template_len))
        f.write("{0:26}\t{1}\n".format("# Confirmed Positions:", len(int_confirmed)))
        f.write("{0:26}\t{1:.4f}\n".format("Confirmed Pos Avg Divergence:", 
                                              av_div))
        f.write("{0:26}\t{1:.2f}\n".format("Mean Confirmed Pos Gap:", 
                                              mean_position_gap))
        f.write("{0:26}\t{1}\n".format("Max Confirmed Pos Gap:", 
                                              max_position_gap))
        f.write("\n\n")
        summ_vals.extend([len(int_confirmed), max_position_gap])
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
            combo_edges = " + ".join(["|".join(["".join(map(str, x)) for x in y]) for y in combo])
            f.write("{0:12}\t{1}\n".format("Resolution:", combo_edges))
            f.write("{0:12}\t{1}\n\n".format("Support:", combo_support[i]))
        #Bridging conditions 
        bridged = False
        bridged_edges = None
        combo_inds = zip(combo_support, range(len(combo_support)))
        sorted_combos = sorted(combo_inds, reverse=True)
        if (sorted_combos[0][0] >= min_bridge_count and 
            sorted_combos[0][0] - sorted_combos[1][0] >= min_bridge_diff):
            bridged = True
            bridged_edges = all_combos[sorted_combos[0][1]]
        best_combo = sorted_combos[0][1]
        best_support = sorted_combos[0][0]
        best_against = best_support
        for support, ind in sorted_combos[1:]:
            best_against -= support
        second_combo = sorted_combos[1][1]
        second_support = sorted_combos[1][0]
        if bridged:
            f.write("BRIDGED\n")
            f.write("Bridging Combo: {0}\n".format(best_combo))
            br_ct_str = "{0} (min_bridge_count)".format(min_bridge_count)
            br_diff_str = "{0} + {1} (Combo {2} + min_bridge_diff)".format(
                second_combo, min_bridge_diff, second_support)
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
            f.write("{0}\t{1}\n".format("min_bridge_count", min_bridge_count))
            f.write("{0}\t{1}\n\n\n".format("min_bridge_diff", min_bridge_diff))
        summ_vals.extend([bridged, best_support, best_against])
        #If not bridged, find in/gap/out lengths and divergence rates
        if not bridged:
            #Write median in, out, and gap lengths
            side_lens = {s:0 for s in repeat_edges[rep]}
            for side in sorted(repeat_edges[rep]):
                trg_limits = []
                for edge_id in sorted(repeat_edges[rep][side]):
                    cons_align = _read_alignment(
                            cons_align_path.format(side_it[side], 
                                                   side, 
                                                   edge_id), 
                            template, 
                            min_aln_rate)
                    if side == "in":
                        trg_limits.append(_get_aln_end(cons_align[0][0].trg_start, cons_align[0][0].trg_seq))
                    elif side == "out":
                        trg_limits.append(template_len - cons_align[0][0].trg_start)
                side_lens[side] = _get_median(trg_limits)
            gap_len = template_len - (side_lens["in"] + side_lens["out"])
            f.write("{0:30}\t{1}\n".format("Median in Sequence Length:", side_lens["in"]))
            f.write("{0:30}\t{1}\n".format("Median out Sequence Length:", side_lens["out"]))
            f.write("{0:30}\t{1}\n\n".format("Median Middle Gap/Overlap Length:", gap_len))
            
            #Write mean in and out divergence rates
            div_rates = {s:[] for s in repeat_edges[rep]}
            for side in sorted(repeat_edges[rep]):
                side_pairs = sorted(combinations(repeat_edges[rep][side], 2))
                for edge_one, edge_two in side_pairs:
                    cons_vs_cons = _read_alignment(
                            cons_vs_cons_path.format(side_it[side], side, edge_one, 
                                                     side_it[side], side, edge_two), 
                            consensuses[(side_it[side], side, edge_two)], 
                            min_aln_rate)
                    edge_rate = _calculate_divergence(cons_vs_cons[0][0].qry_seq, 
                                                      cons_vs_cons[0][0].trg_seq)
                    div_rates[side].append(edge_rate)
            mean_in_div = np.mean(div_rates["in"])
            mean_out_div = np.mean(div_rates["out"])
            weighted_mean_div = (mean_in_div*side_lens["in"] + 
                                 mean_out_div*side_lens["out"]) / float(
                                 side_lens["in"] + side_lens["out"])
            f.write("{0:30}\t{1}\n".format("Mean in Divergence Rate:", mean_in_div))
            f.write("{0:30}\t{1}\n".format("Mean out Divergence Rate:", mean_out_div))
            f.write("{0:30}\t{1}\n\n".format("Weighted Mean Divergence Rate:", 
                                          weighted_mean_div))
            res_str = "No resolution so empty resolved seqs for repeat {0}\n\n"
            f.write(res_str.format(rep))
            for i, edge in enumerate(sorted(repeat_edges[rep]["in"])):
                header = "Repeat_{0}_unbridged_copy_{1}".format(rep, i)
                resolved_repeats[header] = ""
                seq_dict = {header:""}
                fp.write_fasta_dict(seq_dict, resolved_seq_file.format(i))
            summ_vals.extend([""])
        #If bridged, find overlap and construct repeat copy sequences
        else:
            #Find end of repeat as min/max of in/out cons_vs_cons alignments
            edge_limits = {}
            for side in sorted(repeat_edges[rep]):
                side_pairs = sorted(combinations(repeat_edges[rep][side], 2))
                for edge_one, edge_two in side_pairs:
                    cons_vs_cons = _read_alignment(
                            cons_vs_cons_path.format(side_it[side], side, edge_one, 
                                                     side_it[side], side, edge_two), 
                            consensuses[(side_it[side], side, edge_two)], 
                            min_aln_rate)
                    one_start = cons_vs_cons[0][0].qry_start
                    one_end = _get_aln_end(cons_vs_cons[0][0].qry_start, cons_vs_cons[0][0].qry_seq)
                    two_start = cons_vs_cons[0][0].trg_start
                    two_end = _get_aln_end(cons_vs_cons[0][0].trg_start, cons_vs_cons[0][0].trg_seq)
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
                    if side == "in":
                        in_start = edge_limits[(side, edge_id)]
                    elif side == "out":
                        out_end = edge_limits[(side, edge_id)]
                    cons_align = _read_alignment(
                            cons_align_path.format(side_it[side], 
                                                   side, 
                                                   edge_id), 
                            template, 
                            min_aln_rate)
                    if side == "in":
                        in_align = cons_align[0][0]
                    elif side == "out":
                        out_align = cons_align[0][0]
                in_end = _get_aln_end(in_align.qry_start, in_align.qry_seq)
                temp_start = _get_aln_end(in_align.trg_start, in_align.trg_seq)
                temp_end = out_align.trg_start
                out_start = out_align.qry_start
                f.write("Alignment Indices:\n")
                f.write("{0:10}\t{1:5} - {2:5}\n".format("in", in_start, in_end))
                #f.write("{0:10}\t{1:5} - {2:5}\n".format("Template", temp_start, temp_end))
                f.write("{0:10}\t{1:5} - {2:5}\n".format("out", out_start, out_end))
                #Report gap/overlap length
                gap_len = temp_end - temp_start
                if gap_len >= 0:
                    f.write("{0}\t{1}\n".format("Gap between edges:", gap_len))
                else:
                    f.write("{0}\t{1}\n\n".format("Overlap between edges:", -gap_len))
                    #in sequence used to represent overlapping segment
                    #print check of overlapping segment
                    new_temp_end = temp_start
                    new_out_start = None
                    out_qry_aln, out_aln_qry = _index_mapping(out_align.qry_seq)
                    out_trg_aln, out_aln_trg = _index_mapping(out_align.trg_seq)
                    
                    in_edge = edge_pair[0][1]
                    out_edge = edge_pair[1][1]
                    if temp_start >= _get_aln_end(out_align.trg_start, out_align.trg_seq):
                        new_out_start = _get_aln_end(out_align.qry_start, out_align.qry_seq)
                    else:
                        out_aln_ind = out_trg_aln[temp_start]
                        new_out_start = out_start + out_aln_qry[out_aln_ind]
                    """_check_overlap(
                            consensuses[(side_it["in"], "in", in_edge)], 
                            template,
                            consensuses[(side_it["out"], "out", out_edge)], 
                            -gap_len, in_start, in_end, temp_start, temp_end, out_start, out_end,
                            new_out_start, in_align.qry_seq, in_align.trg_seq, out_align.qry_seq, out_align.trg_seq, out_trg_aln, out_aln_trg, out_qry_aln, out_aln_qry, _get_aln_end(out_align.trg_start, out_align.trg_seq), _get_aln_end(out_align.qry_start, out_align.qry_seq), in_align, out_align)
                    """
                    f.write("Adjusted Alignment Indices:\n")
                    f.write("{0:10}\t{1:5} - {2:5}\n".format("in", in_start, in_end))
                    if temp_start != new_temp_end:
                        f.write("{0:10}\t{1:5} - {2:5}\n".format("Template", temp_start, new_temp_end))
                    f.write("{0:10}\t{1:5} - {2:5}\n\n\n".format("out", new_out_start, out_end))
                    temp_end = new_temp_end
                    out_start = new_out_start
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
                        consensuses[(side_it["in"], "in", in_edge)], 
                        template,
                        consensuses[(side_it["out"], "out", out_edge)], 
                        in_start, in_end, 
                        temp_start, temp_end, 
                        out_start, out_end)
                resolved_repeats[header] = copy_seq
                seq_dict = {header:copy_seq}
                fp.write_fasta_dict(seq_dict, resolved_seq_file.format(i))
                in_str = "".join(["in", str(in_edge)])
                out_str = "".join(["out", str(out_edge)])
                summ_resolution.append("|".join([in_str, out_str]))
            summ_vals.extend(["+".join(summ_resolution)])
    return bridged, resolved_repeats, summ_vals

def int_stats_postscript(rep, repeat_edges, integrated_stats, min_aln_rate, 
                         resolved_rep_path, res_vs_res):
    divs = []
    with open(integrated_stats, "a") as f:
        res_inds = range(len(repeat_edges[rep]["in"]))
        f.write("Resolved Repeat Sequence Alignments\n")
        for res_one, res_two in sorted(combinations(res_inds, 2)):
            res_align = _read_alignment(res_vs_res.format(res_one, res_two), 
                                        resolved_rep_path.format(res_two), 
                                        min_aln_rate)
            f.write("Copy {0}|Copy {1}\n".format(res_one, res_two))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format("Copy ", res_one, ":", 
                    res_align[0][0].qry_start, 
                    _get_aln_end(res_align[0][0].qry_start, res_align[0][0].qry_seq), 
                    res_align[0][0].qry_len))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format("Copy ", res_two, ":",   
                    res_align[0][0].trg_start, 
                    _get_aln_end(res_align[0][0].trg_start, res_align[0][0].trg_seq), 
                    res_align[0][0].trg_len))
            div_rate = _calculate_divergence(res_align[0][0].qry_seq, 
                                             res_align[0][0].trg_seq)
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
            for rest in _combo_helper(in_list[1:], out_list[:j]+out_list[j+1:]):
                 yield [combo] + rest

def _get_aln_end(aln_start, aln_seq):
    return aln_start+len(aln_seq.replace("-",""))

def _check_overlap(in_file, temp_file, out_file, overlap, in_start, in_end, temp_start, 
                   temp_end, out_start, out_end, new_out_start, in_qry, in_trg, out_qry, out_trg, out_trg_aln, out_aln_trg, out_qry_aln, out_aln_qry, out_trg_end, out_qry_end, in_align, out_align):
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
    print 'Out old start, old end, new start, out template start', out_start, out_end, new_out_start, temp_end
    print "Out_trg_end", out_trg_end
    print "Out_qry_end", out_qry_end
    print "In align qry inds", in_align.qry_start, in_align.qry_end, in_align.qry_len
    print "In align trg inds", in_align.trg_start, in_align.trg_end, in_align.trg_len
    print "Out align qry inds", out_align.qry_start, out_align.qry_end, out_align.qry_len
    print "Out align trg inds", out_align.trg_start, out_align.trg_end, out_align.trg_len
    print
    print "Overlap:\t{0}".format(overlap)
    print "In seq(-30 to end):\t{0}".format(in_seq[in_end-30:in_end])
    print "Template seq(-30 to end):\t{0}".format(temp_seq[temp_start-30:temp_start])
    #print "Out seq:\t{0}".format(out_seq[out_start:out_end])
    #print "AR In seq:\t{0}".format(in_seq[in_start-10:in_end+10])
    #print "AR Template seq:\t{0}".format(temp_seq[temp_end:temp_start+10])
    #print "AR Out seq:\t{0}".format(out_seq[out_start:out_end+10])
    print "New out seq(-30 to new start):\t{0}".format(out_seq[new_out_start-30:new_out_start])
    print "New out seq(new_start to +30):\t{0}".format(out_seq[new_out_start:new_out_start+30])
    print


def _construct_repeat_copy(in_file, temp_file, out_file, in_start, in_end, 
                           temp_start, temp_end, out_start, out_end):
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

def update_summary(rep, template_len, avg_cov, summ_vals, avg_div, summary_file):
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