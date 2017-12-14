#(c) 2017 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Simple scaffold sequence generator
"""

import flye.fasta_parser as fp
import flye.config as config

def generate_scaffolds(contigs_file, links_file, out_scaffolds):
    def rc(sign):
        return "+" if sign == "-" else "-"
    def unsigned(ctg):
        return ctg[1:]

    contigs_fasta = fp.read_fasta_dict(contigs_file)
    scaffolds_fasta = {}
    used_contigs = set()

    connections = {}
    with open(links_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            ctg_1, sign_1, ctg_2, sign_2 = line.split("\t")
            if ctg_1 in contigs_fasta and ctg_2 in contigs_fasta:
                connections[sign_1 + ctg_1] = sign_2 + ctg_2
                connections[rc(sign_2) + ctg_2] = rc(sign_1) + ctg_1

    scaffolds_fasta = {}
    for ctg in contigs_fasta:
        if ctg in used_contigs: continue

        used_contigs.add(ctg)
        scf = ["-" + ctg]
        #extending right
        while (scf[-1] in connections and
               unsigned(connections[scf[-1]]) not in used_contigs):
            scf.append(connections[scf[-1]])
            used_contigs.add(unsigned(scf[-1]))

        for i, ctg in enumerate(scf):
            scf[i] = rc(ctg[0]) + unsigned(ctg)
        scf = scf[::-1]

        #extending left
        while (scf[-1] in connections and
               unsigned(connections[scf[-1]]) not in used_contigs):
            scf.append(connections[scf[-1]])
            used_contigs.add(unsigned(scf[-1]))

        #generating sequence interleaved by Ns
        if len(scf) == 1:
            scaffolds_fasta[unsigned(ctg)] = contigs_fasta[unsigned(ctg)]
        else:
            scf_name = "scaffold_{0}".format(",".join(scf))
            scf_seq = []
            for scf_ctg in scf:
                if scf_ctg[0] == "+":
                    scf_seq.append(contigs_fasta[unsigned(scf_ctg)])
                else:
                    scf_seq.append(fp.reverse_complement(
                                    contigs_fasta[unsigned(scf_ctg)]))
            gap = "N" * config.vals["scaffold_gap"]
            scaffolds_fasta[scf_name] = gap.join(scf_seq)

    fp.write_fasta_dict(scaffolds_fasta, out_scaffolds)
