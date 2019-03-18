import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--flye-dir")
    parser.add_argument("--input-assembly")
    parser.add_argument("--out-dir")
    parser.add_argument("--min-overlap")
    args = parser.parse_args()

    os.mkdir(args.out_dir)
    open("{}/input_reads.fasta".format(args.out_dir), "w").close()

    repeat_module_binary = "{}/bin/flye-repeat".format(args.flye_dir)
    input_reads = "{}/input_reads.fasta".format(args.out_dir)
    config = "{}/flye/config/bin_cfg/asm_raw_reads.cfg".format(args.flye_dir)

    cmd = [repeat_module_binary]
    cmd.extend([args.input_assembly, input_reads])
    cmd.extend([args.out_dir, config])
    cmd.extend(["-v", args.min_overlap])

    cmd = " ".join(cmd)
    os.system(cmd)


if __name__ == "__main__":
    main()
