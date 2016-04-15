from __future__ import print_function
import sys
import os

import abruijn.alignment as aln
import abruijn.bubbles as bbl
import abruijn.polish as pol
import abruijn.fasta_parser as fp
import abruijn.assemble as asm


def run(reads_file, work_dir, num_proc):
    print("Running ABruijn", file=sys.stderr)
    aln.check_binaries()
    pol.check_binaries()
    asm.check_binaries()

    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    preassembly = os.path.join(work_dir, "read_edges.fasta")
    print("Assembling reads", file=sys.stderr)
    asm.assemble(reads_file, preassembly)
    alignment, genome_len = aln.get_alignment(preassembly, reads_file,
                                              NUM_PROC, work_dir)
    bubbles = bbl.get_bubbles(alignment, genome_len)
    print("Polishing draft assembly", file=sys.stderr)
    polished_seq = pol.polish(bubbles, NUM_PROC, work_dir)
    out_genome = os.path.join(work_dir, "contigs.fasta")
    fp.write_fasta_dict({"contig_1": polished_seq}, out_genome)
    print("Done! Your assembly is in file: " + out_genome, file=sys.stderr)


def main():
    NUM_PROC = 8
    if len(sys.argv) != 3:
        print("Usage: abruijn.py reads_file out_dir")
        return 1

    try:
        run(sys.argv[1], sys.argv[2], NUM_PROC)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException) as e:
        print("Error: ", e, file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
