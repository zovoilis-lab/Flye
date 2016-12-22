#!/usr/bin/env python

from __future__ import print_function
import sys
from collections import defaultdict, namedtuple
from itertools import combinations

from Bio import SeqIO
import networkx as nx


OVERHANG = 1000


class Overlap:
    def __init__(self, cur_id, cur_start, cur_end,
                 ext_id, ext_start, ext_end):
        self.cur_id = cur_id
        self.cur_start = cur_start
        self.cur_end = cur_end
        self.ext_id = ext_id
        self.ext_start = ext_start
        self.ext_end = ext_end
        self.paired = None
        self.knot = None

    def cur_range(self):
        return self.cur_end - self.cur_start

    def cur_intersect(self, other):
        return (min(self.cur_end, other.cur_end) -
                max(self.cur_start, other.cur_start))

    def similar(self, other):
        THRESHOLD = 150
        return (abs(self.cur_start - other.cur_start) < THRESHOLD and
                abs(self.cur_end - other.cur_end) < THRESHOLD and
                abs(self.ext_start - other.ext_start) < THRESHOLD and
                abs(self.ext_end - other.ext_end) < THRESHOLD)


Knot = namedtuple("Knot", ["start_id", "end_id", "length"])
SeqBlock = namedtuple("SeqBlock", ["start", "end", "knot"])


def main():
    if len(sys.argv) != 4:
        print("Usage: build-graph.py contigs overlaps dot_out")
        return 1

    fasta_seqs = {}
    fasta_names = {}
    seq_id = 0
    for seq in SeqIO.parse(sys.argv[1], "fasta"):
        fasta_seqs[seq_id] = seq.seq
        fasta_names[seq_id] = "+" + seq.id
        fasta_seqs[seq_id + 1] = seq.seq.reverse_complement()
        fasta_names[seq_id + 1] = "-" + seq.id
        seq_id += 2

    overlaps = defaultdict(list)
    with open(sys.argv[2], "r") as f:
        for line in f:
            tokens = line.split()
            cur_id, cur_start, cur_end = map(int, tokens[0:3])
            ext_id, ext_start, ext_end = map(int, tokens[4:7])
            ovlp = Overlap(cur_id, cur_start, cur_end,
                           ext_id, ext_start, ext_end)

            exists = False
            for other_ovlp in overlaps[cur_id]:
                if ovlp.similar(other_ovlp):
                    exists = True
                    break
            #if exists:
            #    continue

            paired_ovlp = Overlap(ext_id, ext_start, ext_end,
                                  cur_id, cur_start, cur_end)
            ovlp.paired, paired_ovlp.paired = paired_ovlp, ovlp

            MakeSet(ovlp)
            MakeSet(paired_ovlp)
            Union(ovlp, paired_ovlp)

            overlaps[ovlp.cur_id].append(ovlp)
            overlaps[ovlp.ext_id].append(paired_ovlp)

    #form homology clusters
    for seq in fasta_seqs:
        for o1, o2 in combinations(overlaps[seq], 2):
            if o1.cur_intersect(o2) > 0:
                Union(o1, o2)

    clusters = defaultdict(list)
    for seq in fasta_seqs:
        overlaps[seq].sort(key=lambda o: o.cur_start)
        for ovlp in overlaps[seq]:
            clusters[Find(ovlp)].append(ovlp)
            ovlp.knot = Find(ovlp)
    print(len(clusters))
    ##

    ##define knots
    graph = nx.MultiDiGraph()
    node_id = 0
    knots = {}

    for block in clusters:
        lengths = map(lambda o: o.cur_range(), clusters[block])
        #mean_len = sum(lengths) / len(lengths)
        mean_len = max(lengths)
        graph.add_edge(node_id, node_id + 1,
                       label="repeat, len = {0}".format(mean_len),
                       color="red")
        knots[block] = Knot(node_id, node_id + 1, mean_len)
        node_id += 2
    ###

    #traverse the sequences and build graph
    for seq in fasta_seqs:
        if not overlaps[seq]:
            continue

        #cluster overlapping overlaps
        for ovlp in overlaps[seq]:
            MakeSet(ovlp)
        for o1, o2 in combinations(overlaps[seq], 2):
            if o1.cur_intersect(o2) > 0:
                Union(o1, o2)
        seq_clusters = defaultdict(list)
        for ovlp in overlaps[seq]:
            seq_clusters[Find(ovlp)].append(ovlp)

        blocks = []
        for cl, cl_ovlps in seq_clusters.items():
            start = min(map(lambda o: o.cur_start, cl_ovlps))
            end = max(map(lambda o: o.cur_end, cl_ovlps))
            blocks.append(SeqBlock(start, end, cl_ovlps[0].knot))
        blocks.sort(key=lambda sb: sb.start)

        print(fasta_names[seq], len(blocks), len(overlaps[seq]))

        #sequence start
        if blocks[0].start > OVERHANG:
            graph.add_edge(node_id, knots[blocks[0].knot].start_id,
                           label="{0}[0:{1}]".format(fasta_names[seq],
                                                     blocks[0].start))
            node_id += 1

        #intermediate edges
        for o1, o2 in zip(blocks[:-1], blocks[1:]):
            graph.add_edge(knots[o1.knot].end_id, knots[o2.knot].start_id,
                           label="{0}[{1}:{2}]".format(fasta_names[seq],
                                                       o1.end, o2.start))

        #sequence end
        if len(fasta_seqs[seq]) - blocks[-1].end > OVERHANG:
            graph.add_edge(knots[blocks[-1].knot].end_id, node_id,
                           label="{0}[{1}:{2}]".format(fasta_names[seq],
                                     blocks[-1].end,
                                     len(fasta_seqs[seq])))
            node_id += 1

    nx.write_dot(graph, sys.argv[3])


#disjoint sets implementation
def MakeSet(x):
    x.parent = x
    x.rank = 0


def Union(x, y):
    xRoot = Find(x)
    yRoot = Find(y)
    if xRoot.rank > yRoot.rank:
        yRoot.parent = xRoot
    elif xRoot.rank < yRoot.rank:
        xRoot.parent = yRoot
    elif xRoot != yRoot: # Unless x and y are already in same set, merge them
        yRoot.parent = xRoot
        xRoot.rank = xRoot.rank + 1


def Find(x):
    if x.parent == x:
       return x
    else:
       x.parent = Find(x.parent)
       return x.parent
###############################


if __name__ == "__main__":
    main()
