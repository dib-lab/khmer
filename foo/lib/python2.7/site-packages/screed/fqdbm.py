#!/usr/bin/env python

# Copyright (c) 2008-2010, Michigan State University

import sys
from __init__ import read_fastq_sequences
import DBConstants

# A shell interface to the screed FQDBM database writing function
if __name__ == "__main__":
    # Make sure the user entered the command line arguments correctly
    if len(sys.argv) != 2:
        sys.stderr.write("ERROR: USAGE IS: %s <dbfilename>\n" % sys.argv[0]);
        exit(1)

    filename = sys.argv[1]
    read_fastq_sequences(filename)

    print "Database saved in %s%s" % (sys.argv[1], DBConstants.fileExtension)
    exit(0)
