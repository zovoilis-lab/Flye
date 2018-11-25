#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This module provides repeat gaph parsing/serializing functions,
(as output by repeat graph construction module)
as well as some basic operations
"""

class RgEdge:
    __slots__ = ("node_left", "node_right", "edge_id", "repetitive",
                 "self_complement", "resolved", "mean_coverage",
                 "edge_sequences")

    def __init__(self, node_left=None, node_right=None,
                 edge_id=None):
        self.node_left = node_left
        self.node_right = node_right
        self.edge_id = edge_id

        self.repetitive = False
        self.self_complement = False
        self.resolved = False
        self.mean_coverage = 0
        self.edge_sequences = []

class EdgeSequence:
    __slots__ = ("edge_seq_name", "edge_seq_len", "orig_seq_id", "orig_seq_len",
                 "orig_seq_start", "orig_seq_end")

    def __init__(self, edge_seq_name=None, edge_seq_len=0, orig_seq_id=None,
                 orig_seq_len=0, orig_seq_start=0, orig_seq_end=0):
        self.edge_seq_name = edge_seq_name
        self.edge_seq_len = edge_seq_len

        self.orig_seq_id = orig_seq_id
        self.orig_seq_len = orig_seq_len
        self.orig_seq_start = orig_seq_start
        self.orig_seq_end = orig_seq_end

class RgNode:
    __slots__ = ("in_edges", "out_edge")

    def __init__(self):
        self.in_edges = []
        self.out_edges = []

class RepeatGraph:
    __slots__ = ("nodes", "edges")

    def __init__(self, edges_fasta):
        self.nodes = []
        self.edges = {} #key = edge id
        self.edges_fasta = edges_fasta

    def add_node(self):
        self.nodes.append(RgNode())
        return self.nodes[-1]

    def add_edge(self, edge):
        self.edges[edge.edge_id] = edge
        edge.node_left.out_edges.append(edge)
        edge.node_right.in_edges.append(edge)

    def complement_edge(self, edge):
        return self.edges[-edge.edge_id]

    def load_from_file(self, filename):
        id_to_node = {}
        cur_edge = None
        with open(filename, "r") as f:
            for line in f:
                tokens = line.strip().split()
                if tokens[0] == "Edge":
                    (edge_id, left_node, right_node, repetitive,
                     self_complement, resolved, mean_coverage) = tokens[1:]
                    if left_node not in id_to_node:
                        id_to_node[left_node] = self.add_node()
                    if right_node not in id_to_node:
                        id_to_node[right_node] = self.add_node()

                    cur_edge = RgEdge(id_to_node[left_node],
                                      id_to_node[right_node],
                                      _to_signed_id(int(edge_id)))
                    cur_edge.repetitive = bool(int(repetitive))
                    cur_edge.self_complement = bool(int(self_complement))
                    cur_edge.resolved = bool(int(resolved))
                    cur_edge.mean_coverage = int(mean_coverage)
                    self.add_edge(cur_edge)

                elif tokens[0] == "Sequence":
                    (edge_seq_name, edge_seq_len, orig_seq_id,
                    orig_seq_len, orig_seq_start, orig_seq_end) = tokens[1:]
                    edge_seq = EdgeSequence(edge_seq_name, int(edge_seq_len), orig_seq_id,
                                            int(orig_seq_len), int(orig_seq_start),
                                            int(orig_seq_end))
                    cur_edge.edge_sequences.append(edge_seq)

                else:
                    raise Exception("Error parsing " + filename)

    def dump_to_file(self, filename):
        next_node_id = 0
        node_ids = {}
        for node in self.nodes:
            node_ids[node] = next_node_id
            next_node_id += 1

        with open(filename, "w") as f:
            for edge in self.edges.values():
                f.write("Edge\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n"
                    .format(_to_unsigned_id(edge.edge_id), node_ids[edge.node_left],
                            node_ids[edge.node_right], int(edge.repetitive),
                            int(edge.self_complement), int(edge.resolved),
                            int(edge.mean_coverage)))

                for seq in edge.edge_sequences:
                    f.write("\tSequence\t{0} {1} {2} {3} {4} {5}\n"
                        .format(seq.edge_seq_name, seq.edge_seq_len,
                                seq.orig_seq_id, seq.orig_seq_id,
                                seq.orig_seq_len, seq.orig_seq_start,
                                seq.orig_seq_end))

    def output_dot(self, filename):
        next_node_id = 0
        node_ids = {}
        for node in self.nodes:
            node_ids[node] = next_node_id
            next_node_id += 1

        with open(filename, "w") as f:
            f.write("digraph {\nnodesep = 0.5;\n"
                    "node [shape = circle, label = \"\", height = 0.3];")
            for edge in self.edges.values():
                f.write("{0} -> {1} [label = \"{2}\", color = \"{3}\"]\n"
                    .format(node_ids[edge.node_left], node_ids[edge.node_right],
                            edge.edge_id, "red" if edge.repetitive else "black"))
            f.write("}")


def _to_signed_id(unsigned_id):
    return -(unsigned_id + 1) / 2 if unsigned_id % 2 else unsigned_id / 2 + 1


def _to_unsigned_id(signed_id):
    unsigned_id = abs(signed_id) * 2 - 2;
    return unsigned_id + int(signed_id < 0)
