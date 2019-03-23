import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-assembly")
    parser.add_argument("--out-dir")
    parser.add_argument("--min-overlap", type=int, default=500)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    out_dir = "{}/{}".format(os.getcwd(), args.out_dir)
    flye_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    flye_repeat_bin = "{}/bin/flye-repeat".format(flye_root)

    input_reads = "{}/input_reads.fasta".format(out_dir)
    config = "{}/flye/config/bin_cfg/asm_raw_reads.cfg".format(flye_root)

    cmd = [flye_repeat_bin]
    cmd.extend([args.input_assembly, input_reads])
    cmd.extend([out_dir, config])
    cmd.extend(["-v", str(args.min_overlap)])
    cmd.extend(["-t", str(args.threads)])
    cmd = " ".join(cmd)

    os.mkdir(args.out_dir)
    open(input_reads, "w").close()

    os.system(cmd)

    os.remove("{}/input_reads.fasta".format(out_dir))


if __name__ == "__main__":
    main()
