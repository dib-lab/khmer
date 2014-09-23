#!/usr/bin/env python

# Copyright (c) 2008-2010, Michigan State University

from screed import ToFasta
import sys, os

# Shell interface to the ToFasta screed conversion function
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage: %s <dbfilename> <outputfilename>" % sys.argv[0]
        exit(1)

    dbFile = sys.argv[1]
    outputFile = sys.argv[2]

    if not os.path.isfile(dbFile):
        print "No such file: %s" % dbFile
        exit(1)
    if os.path.isfile(outputFile):
        os.unlink(outputFile)
    
    ToFasta(dbFile, outputFile)
