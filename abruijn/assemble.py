import subprocess

from abruijn.utils import which

ASSEMBLE_BIN = "abruijn-assemble"


class AssembleException(Exception):
    pass


def check_binaries():
    if not which(ASSEMBLE_BIN):
        raise AssembleException("assemble binary not found. Did you run 'make'?")


def assemble(reads_file, out_file):
    try:
        subprocess.check_call([ASSEMBLE_BIN, reads_file, out_file])
    except (subprocess.CalledProcessError, OSError) as e:
        raise AssembleException("Error in assemble binary: " + e)
