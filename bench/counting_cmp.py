#! /usr/bin/env python

import sys

import khmer

cg = khmer.Countgraph(31, 1e9, 4)
cg.consume_fasta(sys.argv[1])
