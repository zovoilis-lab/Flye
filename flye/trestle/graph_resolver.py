#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Modifies repeat graph using the Tresle output
"""

import logging
from itertools import izip, chain
from collections import defaultdict

import flye.utils.fasta_parser as fp

logger = logging.getLogger()


class Connection:
    __slots__ = ("id", "path", "sequence")
    def __init__(self, id=None, path=[], sequence=""):
        self.id = id
        self.path = path
        self.sequence = sequence


class RepeatInfo:
    __slots__ = ("id", "repeat_path", "all_reads", "in_reads", "out_reads",
                 "sequences", "multiplicity")

    def __init__(self, id, repeat_path, all_reads, in_reads, out_reads,
                 sequences, multiplicity):
        self.id = id
        self.repeat_path = repeat_path
        self.all_reads = all_reads
        self.in_reads = in_reads
        self.out_reads = out_reads
        self.sequences = sequences
        self.multiplicity = multiplicity


def get_simple_repeats(repeat_graph, alignments, edge_seqs):
    next_path_id = 1
    path_ids = {}
    repeats_dict = {}
    MULT = 2

    for path in repeat_graph.get_unbranching_paths():
        if not path[0].repetitive or path[0].self_complement:
            continue

        is_simple = True
        inputs = set()
        for in_edge in path[0].node_left.in_edges:
            inputs.add(in_edge.edge_id)
            if in_edge.repetitive:
                is_simple = False

        outputs = set()
        for out_edge in path[-1].node_right.out_edges:
            outputs.add(out_edge.edge_id)
            if out_edge.repetitive:
                is_simple = False

        if not is_simple or len(inputs) != MULT or len(outputs) != MULT:
            continue

        if path[0].edge_id not in path_ids:
            path_ids[path[0].edge_id] = next_path_id
            path_ids[-path[-1].edge_id] = -next_path_id
            next_path_id += 1
        path_id = path_ids[path[0].edge_id]

        repeat_edge_ids = set(map(lambda e: e.edge_id, path))
        inner_reads = []
        input_reads = defaultdict(list)
        output_reads = defaultdict(list)
        for read_aln in alignments:
            repeat_read = False
            for edge_aln in read_aln:
                if edge_aln.edge_id in repeat_edge_ids:
                    repeat_read = True
            if not repeat_read:
                continue

            inner_reads.append(read_aln[0].overlap.cur_id)
            for prev_edge, next_edge in izip(read_aln[:-1], read_aln[1:]):
                if (prev_edge.edge_id in inputs and
                        next_edge.edge_id == path[0].edge_id):
                    input_reads[prev_edge.edge_id].append(prev_edge.overlap.cur_id)

                if (prev_edge.edge_id == path[-1].edge_id and
                        next_edge.edge_id in outputs):
                    output_reads[next_edge.edge_id].append(next_edge.overlap.cur_id)


        #TODO: temporary restriction
        #if len(path) != 1:
        #    continue

        #repeats_dict[path_id] = (MULT, inner_reads, input_reads,
        #                         output_reads, map(lambda e: e.edge_id, path))

        #add edges sequences:
        print path_id
        sequences = {}
        for edge in chain(input_reads, output_reads):
            seq_id = repeat_graph.edges[edge].edge_sequences[0].edge_seq_name
            seq = edge_seqs[seq_id[1:]]
            if seq_id[0] == "-":
                seq = fp.reverse_complement(seq)
            sequences[edge] = seq

        template_seq = ""
        for edge in path:
            seq_id = edge.edge_sequences[0].edge_seq_name
            seq = edge_seqs[seq_id[1:]]
            if seq_id[0] == "-":
                seq = fp.reverse_complement(seq)
            template_seq += seq
        sequences["template"] = template_seq

        for h, s in sequences.items():
            print h, s[:100]

        repeats_dict[path_id] = RepeatInfo(path_id, map(lambda e: e.edge_id, path),
                                           inner_reads, input_reads, output_reads,
                                           sequences, MULT)

    return repeats_dict


def dump_repeats(repeats_info, filename):
    with open(filename, "w") as f:
        for repeat_id, info in repeats_info.iteritems():
            f.write("#Repeat {0}\t{1}\n\n".format(repeat_id, info.multiplicity))

            f.write("#All reads\t{0}\n".format(len(info.all_reads)))
            for read in info.all_reads:
                f.write(read + "\n")
            f.write("\n")

            for in_edge in info.in_reads:
                f.write("#Input {0}\t{1}\n".format(in_edge, len(info.in_reads[in_edge])))
                for read in info.in_reads[in_edge]:
                    f.write(read + "\n")
                f.write("\n")

            for out_edge in info.out_reads:
                f.write("#Output {0}\t{1}\n".format(out_edge, len(info.out_reads[out_edge])))
                for read in info.out_reads[out_edge]:
                    f.write(read + "\n")
                f.write("\n")


def apply_changes(repeat_graph, trestle_results,
                  resolved_repeats_fasta):
    #repeat_graph.output_dot("before.dot")
    connections = _get_connections(trestle_results)
    for conn in connections:
        repeat_graph.separate_path(conn.path, conn.id,
                                   resolved_repeats_fasta[conn.sequence])

    #repeat_graph.output_dot("after.dot")


def _get_connections(trestle_results):
    connections = []
    resolved_repeats = set()
    with open(trestle_results, "r") as f:
        for line in f:
            if line.startswith("Repeat"): continue

            tokens = line.strip().split()
            repeat_id, bridged = int(tokens[0]), tokens[5]
            if bridged == "True" and abs(repeat_id) not in resolved_repeats:
                resolved_repeats.add(abs(repeat_id))

                path_1, path_2 = tokens[9].split(":")
                seq_1, seq_2 = tokens[10].split(":")
                in_1, out_1 = path_1.split(",")
                in_2, out_2 = path_2.split(",")
                connection_1 = [int(in_1), repeat_id, int(out_1)]
                connection_2 = [int(in_2), repeat_id, int(out_2)]

                logger.info("Repeat {0}: {1}, {2}"
                    .format(repeat_id, connection_1, connection_2))

                new_seq_id = "trestle_resolved_" + str(repeat_id) + "_copy_"
                connections.extend([Connection(new_seq_id + "1", connection_1, seq_1),
                                    Connection(new_seq_id + "2", connection_2, seq_2)])

    return connections
