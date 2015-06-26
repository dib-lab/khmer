#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

"""File handling/checking utilities for command-line scripts."""

from __future__ import print_function, unicode_literals

import os
import sys
import errno
from stat import S_ISBLK, S_ISFIFO


def check_input_files(file_path, force):
    """Check the status of the file.

    If the file is empty or doesn't exist AND if the file is NOT a
    fifo/block/named pipe then a warning is printed and sys.exit(1) is
    called
    """
    mode = None

    if file_path == '-':
        return
    try:
        mode = os.stat(file_path).st_mode
    except OSError:
        print("ERROR: Input file %s does not exist" %
              file_path, file=sys.stderr)

        if not force:
            print("Exiting", file=sys.stderr)
            sys.exit(1)
        else:
            return

    # block devices will be nonzero
    if S_ISBLK(mode) or S_ISFIFO(mode):
        return

    if not os.path.exists(file_path):
        print("ERROR: Input file %s does not exist; exiting" %
              file_path, file=sys.stderr)
        if not force:
            sys.exit(1)
    else:
        if os.stat(file_path).st_size == 0:
            print("ERROR: Input file %s is empty; exiting." %
                  file_path, file=sys.stderr)
            if not force:
                sys.exit(1)


def check_file_writable(file_path):
    """Return if file_path is writable, exit out if it's not."""
    try:
        file_obj = open(file_path, "a")
    except IOError as error:
        if error.errno == errno.EACCES:
            print("ERROR: File %s does not have write "
                  % file_path + "permission; exiting", file=sys.stderr)
            sys.exit(1)
        else:
            print("ERROR: " + error.strerror, file=sys.stderr)
    else:
        file_obj.close()
        return


def check_space(in_files, force, _testhook_free_space=None):
    """Check for available disk space.

    Estimate size of input files passed, then calculate disk space
    available and exit if disk space is insufficient.
    """
    # Get disk free space in Bytes assuming non superuser
    # and assuming all inFiles are in same disk
    in_file = in_files[0]

    dir_path = os.path.dirname(os.path.realpath(in_file))
    target = os.statvfs(dir_path)

    if _testhook_free_space is None:
        free_space = target.f_frsize * target.f_bavail
    else:
        free_space = _testhook_free_space

    # Check input file array, remove corrupt files
    valid_files = [f for f in in_files if os.path.isfile(f)]

    # Get input file size as worst case estimate of
    # output file size
    file_sizes = [os.stat(f).st_size for f in valid_files]
    total_size = sum(file_sizes)

    size_diff = total_size - free_space
    if size_diff > 0:
        print("ERROR: Not enough free space on disk "
              "for output files;\n"
              "       Need at least %.1f GB more."
              % (float(size_diff) / 1e9), file=sys.stderr)
        print("       Estimated output size: %.1f GB"
              % (float(total_size) / 1e9,), file=sys.stderr)
        print("       Free space: %.1f GB"
              % (float(free_space) / 1e9,), file=sys.stderr)
        if not force:
            sys.exit(1)


def check_space_for_hashtable(hash_size, force, _testhook_free_space=None):
    """Check we have enough size to write a hash table."""
    cwd = os.getcwd()
    dir_path = os.path.dirname(os.path.realpath(cwd))
    target = os.statvfs(dir_path)
    if _testhook_free_space is None:
        free_space = target.f_frsize * target.f_bavail
    else:
        free_space = _testhook_free_space  # allow us to test this code...

    size_diff = hash_size - free_space
    if size_diff > 0:
        print("ERROR: Not enough free space on disk "
              "for saved table files;"
              "       Need at least %s GB more."
              % (float(size_diff) / 1e9,), file=sys.stderr)
        print("       Table size: %.1f GB"
              % (float(hash_size) / 1e9,), file=sys.stderr)
        print("       Free space: %.1f GB"
              % (float(free_space) / 1e9,), file=sys.stderr)
        if not force:
            sys.exit(1)


def check_valid_file_exists(in_files):
    """Warn if input files are empty or missing.

    In a scenario where we expect multiple input files and
    are OK with some of them being empty or non-existent,
    this check warns to stderr if any input file is empty
    or non-existent.
    """
    for in_file in in_files:
        if os.path.exists(in_file):
            if os.stat(in_file).st_size > 0:
                return
            else:
                print('WARNING: Input file %s is empty' %
                      in_file, file=sys.stderr)
        else:
            print('WARNING: Input file %s not found' %
                  in_file, file=sys.stderr)
