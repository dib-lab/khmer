#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

import os
import sys


def check_file_status(file_path):
    """ Check status of file - return if file exists,
    is empty, or neither """
    if not os.path.exists(file_path):
        print >>sys.stderr, 'ERROR: Input file %s does not exist,\
         exiting' % file_path
        sys.exit(1)
    else:
        if os.stat(file_path).st_size == 0:
            print >>sys.stderr, 'ERROR: Input file %s is empty,\
                     exiting' % file_path
            sys.exit(1)


def check_space(in_files):
    """ Estimate size of input files passed, then calculate
    disk space available. Exit if insufficient disk space,
    """

    # Get disk free space in Bytes assuming non superuser
    # and assuming all inFiles are in same disk
    dir_path = os.path.dirname(in_files[0])
    target = os.statvfs(dir_path)
    free_space = target.f_frsize * target.f_bavail
    #<TODO>: If SU, use target.f_bfree

    # Check input file array, remove corrupt files
    valid_files = [f for f in in_files if os.path.isfile(f)]

    # Get input file size as worst case estimate of
    # output file size
    file_sizes = map(lambda f: os.stat(f).st_size, valid_files)
    total_size = reduce(lambda f1, f2: f1+f2, file_sizes)

    size_diff = total_size-free_space
    if size_diff > 0:
        print >>sys.stderr, 'ERROR: Not enough free space on disk, \
        need at least %s more,' % str(size_diff)
        sys.exit(1)


def check_space_for_hashtable(hash_size):
    """
    Check we have enough size to write a hash table
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    target = os.statvfs(dir_path)
    free_space = target.f_frsize * target.f_bavail

    size_diff = hash_size-free_space
    if size_diff > 0:
        print >>sys.stderr, 'ERROR: Not enough free space on disk, \
        need at least %s more,' % str(size_diff)
        sys.exit(1)


def check_valid_file_exists(in_files):
    """
    In a scenario where we expect multiple input files and
    are OK with some of them being empty or non-existent,
    this check ensures there is at least one valid file
    """
    for in_file in in_files:
        if os.path.exists(in_file):
            if os.stat(in_file).st_size > 0:
                return
            else:
                print >>sys.stderr, 'WARNING: Input file %s is empty,\
                     ' % in_file
        else:
            print >>sys.stderr, 'WARNING: Input file %s not found,\
                     ' % in_file
