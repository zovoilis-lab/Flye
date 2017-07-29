#!/usr/bin/env python

#(c) 2017 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
A script for easy BLASR installation
"""

from __future__ import print_function
import sys, os, stat
import subprocess
import shutil
import argparse
import platform
try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve


HDF5_GIT = "https://github.com/mortenpi/hdf5/"
HDF5_REVISION = "8a275ab"

BLASR_GIT = "https://github.com/PacificBiosciences/blasr"
BLASR_REVISION = "eb2056b"

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
INSTALL_DIR = os.path.join(ROOT_DIR, "lib", "blasr_install")
BIN_DIR = os.path.join(ROOT_DIR, "bin")

HDF5_ROOT = os.path.join(INSTALL_DIR, "hdf5")
HDF5_PREFIX = os.path.join(HDF5_ROOT, "install")


def install_all():
    blasr_path = os.path.join(BIN_DIR, "blasr")
    try:
        if os.path.isdir(INSTALL_DIR):
            shutil.rmtree(INSTALL_DIR)
        os.mkdir(INSTALL_DIR)

        install_hdf()
        install_blasr()
        make_wrapper()

        wrapper_path = os.path.join(BIN_DIR, "blasr")
        print("\nSuccessfully installed to {0}".format(wrapper_path))

    except OSError as e:
        print("Error while installing - exiting", file=sys.stderr)
        print(e, file=sys.stderr)
        return False


def install_hdf():
    print("Installing HDF5", file=sys.stderr)
    os.chdir(INSTALL_DIR)
    subprocess.check_call(["git", "clone", HDF5_GIT])
    os.chdir("hdf5")
    subprocess.check_call(["git", "reset", "--hard", HDF5_REVISION])
    subprocess.check_call(["./configure", "--prefix={0}".format(HDF5_PREFIX),
                           "--enable-cxx"])
    subprocess.check_call(["make"])
    subprocess.check_call(["make", "install"])


def install_blasr():
    print("Installing BLASR", file=sys.stderr)

    hdf5_include = os.path.join(HDF5_PREFIX, "include")
    hdf5_lib = os.path.join(HDF5_PREFIX, "lib")

    os.chdir(INSTALL_DIR)
    subprocess.check_call(["git", "clone", BLASR_GIT])
    os.chdir("blasr")
    subprocess.check_call(["git", "pull", "--rebase", "origin", "master"])
    subprocess.check_call(["git", "reset", "--hard", BLASR_REVISION])
    subprocess.check_call(["make", "update-submodule"])
    subprocess.check_call(["./configure.py", "--shared", "--sub", "--no-pbbam",
                           "HDF5_INCLUDE={0}".format(hdf5_include),
                           "HDF5_LIB={0}".format(hdf5_lib)])
    subprocess.check_call(["make", "configure-submodule"])
    subprocess.check_call(["make", "build-submodule"])
    subprocess.check_call(["make", "blasr"])


def make_wrapper():
    blasr_bin = os.path.join(INSTALL_DIR, "blasr", "blasr")
    wrapper_path = os.path.join(BIN_DIR, "blasr")
    hdf5_lib = os.path.join(HDF5_PREFIX, "lib")

    if platform.system() == "Darwin":
        library_env = "DYLD_LIBRARY_PATH"
    else:
        library_env = "LD_LIBRARY_PATH"

    cur_env = ""
    try:
        cur_env = os.environ[library_env]
    except KeyError:
        pass

    ld_library_path = (hdf5_lib + ":" +
                      os.path.join(INSTALL_DIR, "blasr", "libcpp", "alignment") + ":" +
                      os.path.join(INSTALL_DIR, "blasr", "libcpp", "hdf") + ":" +
                      os.path.join(INSTALL_DIR, "blasr", "libcpp", "pbdata") + ":" +
                      cur_env)

    script_text = ("#!/bin/bash\n" +
                   "export {0}={1}\n".format(library_env, ld_library_path) +
                   "if [ ! -f {0} ]; then\n".format(blasr_bin) +
                   "\techo \"Something happened since the last BLASR installation. "+
                   "Please re-run install_blasr.py\"\n" +
                   "\texit 1\n"
                   "fi\n"
                   "{0} $@".format(blasr_bin))

    open(wrapper_path, "w").write(script_text)
    st = os.stat(wrapper_path)
    os.chmod(wrapper_path, st.st_mode | stat.S_IEXEC)


#Mimics UNIX "which" command
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def test_requirements():
    for tool in ["git", "make"]:
        if not which(tool):
            print("ERROR: building BLASR requires " + tool, file=sys.stderr)
            return False
    return True


def main():
    if not test_requirements():
        return 1

    if not install_all():
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
