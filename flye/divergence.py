# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 03:50:31 2017

@author: jeffrey_yuan
"""

import logging
from collections import defaultdict
from itertools import izip
import multiprocessing
import signal
import copy
import numpy as np

from flye.alignment import shift_gaps, SynchronizedSamReader
import flye.config as config
import flye.fasta_parser as fp

logger = logging.getLogger()

class Profile:
    __slots__ = ("insertions", "matches", "nucl")

    def __init__(self):
        self.insertions = defaultdict(int)
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
    #max_aln_err = config.vals["err_modes"][platform]["max_aln_error"]
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
                prof_elem.insertions[qry_nuc] += 1
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    return profile, aln_errors

def _count_freqs(elem):
    matches = elem.matches
    insertions = elem.insertions
    nucl = elem.nucl
    del_key = "-"
    
    coverage = sum(matches.values())
    
    mat_ct = 0
    if nucl in matches:
        mat_ct = matches[nucl]
        
    subs = {key:matches[key] for key in matches 
                                if (key != nucl and key != del_key)}
    max_sub_ct = 0
    max_sub_base = ""
    if subs:
        max_sub_base = max(subs, key=subs.get)
        max_sub_ct = subs[max_sub_base]
    
    del_ct = 0
    if del_key in matches:
        del_ct = matches[del_key]
    
    max_ins_key = "^"
    max_ins_ct = 0
    if insertions:
        max_ins_base = max(insertions, key=insertions.get)
        max_ins_ct = insertions[max_ins_base]
        max_ins_key = "^{0}".format(max_ins_base)
    
    return {'cov':coverage, 'mat_base':nucl, 'mat_ct':mat_ct, 
                            'sub_base':max_sub_base, 'sub_ct':max_sub_ct, 
                            'del_base':del_key, 'del_ct':del_ct, 
                            'ins_base':max_ins_key, 'ins_ct':max_ins_ct}

def _call_position(ind, counts, pos, called, sub_thresh, del_thresh, ins_thresh):
    over_thresh = False
    if counts['cov']:
        if counts['sub_ct'] / float(counts['cov']) >= sub_thresh:
            called['sub'] += 1
            over_thresh = True
        if counts['del_ct'] / float(counts['cov']) >= del_thresh:
            called['del'] += 1
            over_thresh = True
        if counts['ins_ct'] / float(counts['cov']) >= ins_thresh:
            called['ins'] += 1
            over_thresh = True
        
        if over_thresh:
            called['total'] += 1
            pos[ind] = counts
        
    return pos, called

def find_divergence(alignment_path, contigs_path, contigs_info, 
                    frequency_path, positions_path, div_sum_path, 
                    min_aln_rate, platform, num_proc, 
                    sub_thresh, del_thresh, ins_thresh):
    """
    Main function: takes in an alignment and finds the divergent positions
    """
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_fasta_dict(contigs_path),
                                       min_aln_rate)
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

    out_positions = []
    total_aln_errors = []
    while not results_queue.empty():
        ctg_id, ctg_profile, aln_errors = results_queue.get()
        
        positions, total_called = _write_frequency_path(frequency_path, 
                                                        ctg_profile,
                                                        sub_thresh, 
                                                        del_thresh, 
                                                        ins_thresh)
        pos_header = "".join(["Called_{0}_positions_".format(len(positions)), 
                              "with_thresholds_sub_{0}".format(sub_thresh), 
                              "_del_{0}_ins{1}".format(del_thresh, ins_thresh)])
        write_positions(positions_path, pos_header, positions)
        
        window_len = 1000
        sum_header = "Tentative Divergent Position Summary"
        write_div_summary(div_sum_path, sum_header, positions, total_called, 
                          len(ctg_profile), window_len)
        
        total_aln_errors.extend(aln_errors)
        if len(positions) > 0:
            out_positions.extend(positions)
            
    mean_aln_error = float(sum(total_aln_errors)) / (len(total_aln_errors) + 1)
    logger.debug("Total positions: {0}".format(len(out_positions)))
    logger.debug("Alignment error rate: {0}".format(mean_aln_error))

def _write_frequency_path(frequency_path, ctg_profile, sub_thresh, 
                          del_thresh, ins_thresh):
    #The set of called positions
    positions = {}
    #The total number of positions called for each category
    total_called = {'total':0, 'sub':0, 'del':0, 'ins':0}
    with open(frequency_path, 'w') as f:
        f.write("Index\tCov\tMatch\tCount\tSub\tCount\tDel\tCount\tIns\tCount\n")
        for index, elem in enumerate(ctg_profile):
            counts = _count_freqs(elem)
            f.write("{0}\t{c[cov]}\t{c[mat_base]}\t{c[mat_ct]}\t".format(index,c=counts))
            f.write("{c[sub_base]}\t{c[sub_ct]}\t{c[del_base]}\t".format(c=counts))
            f.write("{c[del_ct]}\t{c[ins_base]}\t{c[ins_ct]}\n".format(c=counts))
            
            #Adds this element to positions if it passes any threshold
            #and updates total_called appropriately
            positions, total_called = _call_position(index, counts, positions, 
                        total_called, sub_thresh, del_thresh, ins_thresh)
    return positions, total_called
    
def write_positions(positions_path, pos_header, positions):
    with open(positions_path, 'w') as f:
        f.write(">{0}\n".format(pos_header))
        f.write(",".join(map(str, sorted(positions.keys()))))
        f.write("\n")

def write_div_summary(div_sum_path, sum_header, positions, total_called, 
                      seq_len, window_len):
    pos_list = sorted(positions.keys())
    av_div = 0.0
    if seq_len != 0:
        av_div = len(positions) / float(seq_len)
    
    position_gaps = [0 for x in range(len(positions)+1)]
    curr_pos = 0
    for i,p in enumerate(pos_list):
        position_gaps[i] = p-curr_pos
        curr_pos = p
    position_gaps[-1] = seq_len-curr_pos
    
    mean_position_gap = np.mean(position_gaps)
    max_position_gap = max(position_gaps)
    
    window_len = 1000
    position_counts = [0 for x in range(((seq_len-1)/window_len)+1)]
    window_divs = [0.0 for x in range(((seq_len-1)/window_len)+1)]
    curr_p_i = 0
    for i in range(len(window_divs)):
        start = i*window_len
        end = (i+1)*window_len-1
        if i == len(window_divs)-1:
            end = seq_len-1
        
        curr_window_len = end-start+1
        
        if curr_p_i < len(pos_list) and pos_list[curr_p_i] < start:
            raise PositionIOError('Problem with position indices')
        while curr_p_i < len(pos_list) and pos_list[curr_p_i] <= end:
            position_counts[i] += 1
            curr_p_i += 1
        
        window_divs[i] = 0.0
        if curr_window_len != 0:
            window_divs[i] = position_counts[i]/float(curr_window_len)
    
    mean_window_div = np.mean(window_divs)
    median_window_div = np.median(window_divs)
    min_window_div = min(window_divs)
    
    with open(div_sum_path, 'w') as f:
        f.write("{0}\n\n".format(sum_header))
        
        f.write("Sequence Length:\t{0}\n".format(seq_len))
        f.write("Total Positions:\t{0}\n".format(len(positions)))
        f.write("Average Divergence:\t{:.4f}\n\n".format(av_div))
        
        f.write("Total Substitution Positions:\t{t[sub]}\n".format(t=total_called))
        f.write("Total Deletion Positions:\t{t[del]}\n".format(t=total_called))
        f.write("Total Insertion Positions:\t{t[ins]}\n".format(t=total_called))
        f.write("Total Positions:\t{t[total]}\n".format(t=total_called))
        mixed_count = (total_called['sub'] + total_called['del'] + 
                    total_called['ins']) - total_called['total']
        f.write("Mixed Positions:\t{0}\n\n".format(mixed_count))
        
        f.write('Mean Position Gap:\t%.2f\n' % mean_position_gap)
        f.write('Max Position Gap:\t%d\n\n' % max_position_gap)
        
        f.write('Window Length\t%d\n' % window_len)
        f.write('Mean Window Divergence:\t%.5f\n' % mean_window_div)
        f.write('Median Window Divergence:\t%.5f\n' % median_window_div)
        f.write('Min Window Divergence:\t%.5f\n\n' % min_window_div)
    
    

class PositionIOError(Exception):
    pass

def read_positions(positions_file):
    """
    Reads positions file into list
    """
    header = None
    positions = []
    try:
        with open(positions_file, "r") as f:
            for line_id, line in enumerate(f):
                line = line.strip()
                if line.startswith(">"):
                    header = line[1:]
                elif line:
                    pos_parts = line.split(",")
                    for p in pos_parts:
                        if p.strip():
                            positions = map(int, pos_parts)
    except IOError as e:
        raise PositionIOError(e)
    return positions


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
    
    
    filtered_pairs, empty_pairs = _filter_pairs(pairs, thresh, 
                                                empty_pairs)
    
    graph_file = "adjacency_graph_filtered_{0}.txt".format(ctg_id)
    _write_graph(filtered_pairs, graph_file, positions, qry_pos_bases)
    
    
    
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
                pairs[p1][p2] = {"cis":{(p1_base1, p2_base1) : 0, 
                                  (p1_base2, p2_base2) : 0},
                                 "trans":{(p1_base1, p2_base2) : 0,
                                  (p1_base2, p2_base1) : 0}}
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
                        if (b1, b2) in pairs[p1][p2]["cis"]:
                            pairs[p1][p2]["cis"][(b1, b2)] += 1
                        elif (b1, b2) in pairs[p1][p2]["trans"]:
                            pairs[p1][p2]["trans"][(b1, b2)] += 1
            if sum(pairs[p1][p2]["cis"].values() + pairs[p1][p2]["trans"].values()) == 0:
                empty_pairs.add((p1, p2))
    return pairs, empty_pairs

def _write_graph(pairs, graph_file, pos, qry_pos_bases):
    with open(graph_file, "w") as f:
        f.write("P1\\P2\t")
        for p1 in sorted(pairs):
            #f.write("{0} {1}|{2} {3} {4}\t".format(p1, pos[p1][0], pos[p1][1], 
            #                                   len(qry_pos_bases[p1]), qry_pos_bases[p1]))
            f.write("{0} {1}|{2} {3}\t".format(p1, pos[p1][0], pos[p1][1], 
                                               len(qry_pos_bases[p1])))
        f.write("\n")
        
        for p1 in sorted(pairs):
            #f.write("{0} {1}|{2} {3} {4}\t".format(p1, pos[p1][0], pos[p1][1], 
            #                                   len(qry_pos_bases[p1]), qry_pos_bases[p1]))
            f.write("{0} {1}|{2} {3}\t".format(p1, pos[p1][0], pos[p1][1], 
                                               len(qry_pos_bases[p1]))) 
            for p2 in sorted(pairs):
                if p2 in pairs[p1]:
                    #f.write("{0}|{1}  ".format(p1, p2))
                    if "cis" in pairs[p1][p2] and "trans" in pairs[p1][p2]:
                        for b1, b2 in sorted(pairs[p1][p2]["cis"]):
                            count = pairs[p1][p2]["cis"][(b1, b2)]
                            if count > 0:
                                f.write("{0},{1}:{2} ".format(b1, b2, count))
                        f.write("vs ")
                        #if "trans" in pairs[p1][p2]:
                        for b1, b2 in sorted(pairs[p1][p2]["trans"]):
                            count = pairs[p1][p2]["trans"][(b1, b2)]
                            if count > 0:
                                f.write("{0},{1}:{2} ".format(b1, b2, count))
                    f.write("\t")
                else:
                    f.write("\t")
            f.write("\n")
                

def _filter_pairs(pairs, thresh, empty_pairs):
    filtered_pairs = {}
    for p1 in pairs:
        filtered_pairs[p1] = {}
        for p2 in pairs[p1]:
            filtered_pairs[p1][p2] = {}
            cis_count = sum(pairs[p1][p2]["cis"].values())
            trans_count = sum(pairs[p1][p2]["trans"].values())
            if (p1, p2) not in empty_pairs:
                cis_frac = 0.0
                trans_frac = 0.0
                if cis_count + trans_count != 0:
                    cis_frac = cis_count/float(cis_count + trans_count)
                    trans_frac = trans_count/float(cis_count + trans_count)
                if cis_frac >= thresh:
                    filtered_pairs[p1][p2]["cis"] = pairs[p1][p2]["cis"]
                if trans_frac >= thresh:
                    filtered_pairs[p1][p2]["trans"] = pairs[p1][p2]["trans"]
            if not filtered_pairs[p1][p2]:
                empty_pairs.add((p1, p2))
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
    
def confirm_positions(alignment_path, contigs_path, contigs_info, min_aln_rate,
                    platform, num_proc, positions_file, confirm_thresh):
    """
    Uses frequency of base occurrences at pairs of positions in reads
    to confirm or reject positions. Can be used iteratively.
    """
    orig_positions = read_positions(positions_file)
    positions = _match_positions_to_contigs(orig_positions,contigs_info)
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_fasta_dict(contigs_path),
                                       min_aln_rate)
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