#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: khmer-project@idyll.org
#

import os
import sys


def check_file_status(filePath):
    """ Check status of file - return if file exists,
    is empty, or neither """
    if not os.path.exists(filePath):
        print >>sys.stderr, 'ERROR: Input file %s does not exist,\
         exiting' % filePath
        sys.exit(1)
    else:
        if os.stat(filePath).st_size == 0:
            print >>sys.stderr, 'ERROR: Input file %s is empty,\
                     exiting' % filePath
            sys.exit(1)


def check_space(inFiles):
    """ Estimate size of inFiles passed, then calculate
    disk space available. Exit if insufficient disk space,
    """

    # Get disk free space in Bytes assuming non superuser
    # and assuming all inFiles are in same disk
    dirPath = os.path.dirname(inFiles[0])
    target = os.statvfs(dirPath)
    freeSpace = target.f_frsize * target.f_bavail
    #<TODO>: If SU, use target.f_bfree

    # Get input file size as worst case estimate of
    # output file size
    fileSizes = map(lambda f: os.stat(f).st_size, inFiles)
    totalSize = reduce(lambda f1, f2: f1+f2, fileSizes)

    sizeDiff = totalSize-freeSpace
    if sizeDiff > 0:
        print >>sys.stderr, 'ERROR: Not enough free space on disk, \
        need at least %s more,' % str(sizeDiff)
        sys.exit(1)


def check_space_for_hashtable(hashSize):
    """
    Check we have enough size to write a hash table
    """
    dirPath = os.path.dirname(os.path.realpath(__file__))
    target = os.statvfs(dirPath)
    freeSpace = target.f_frsize * target.f_bavail

    sizeDiff = hashSize-freeSpace
    if sizeDiff > 0:
        print >>sys.stderr, 'ERROR: Not enough free space on disk, \
        need at least %s more,' % str(sizeDiff)
        sys.exit(1)
