#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys

import khmer
import screed

parser = argparse.ArgumentParser()
parser.add_argument('--controls', metavar='FILE', nargs='+',
                    help='list of countgraph files corresponding to control '
                    'sample(s)')
parser.add_argument('--case', metavar='FILE', help='countgraph of '
                    'sample from case/condition of interest')
parser.add_argument('-x', '--control_threshold', metavar='X', type=int,
                    default=1, help='k-mers with abund > X in any control '
                    'sample are uninteresting; default=1')
parser.add_argument('-y', '--case_lower_threshold', metavar='Y', type=int,
                    default=5, help='ignore k-mers from case with abund < Y; '
                    'default=5')
parser.add_argument('case_fastq')
args = parser.parse_args()

print('Loading case countgraph', args.case, '...', file=sys.stderr)
case = khmer.load_countgraph(args.case)
print('Case countgraph loaded, k={:d}'.format(case.ksize()), file=sys.stderr)

controls = list()
for ctlfile in args.controls:
    print('Loading control countgraph', ctlfile, '...', file=sys.stderr)
    countgraph = khmer.load_countgraph(ctlfile)
    assert countgraph.ksize() == case.ksize()
    controls.append(countgraph)
    print('Control countgraph loaded', file=sys.stderr)

print('Iterating over case reads', args.case_fastq, '...', sys.stderr)
for n, record in enumerate(screed.open(args.case_fastq)):
    if n+1 % 1e6 == 0:
        print('    processed', n+1, 'reads...', file=sys.stderr)
    for kmer in case.get_kmers(record.sequence):
        case_abund = case.get_count(kmer)
        if case_abund < args.case_lower_threshold:
            continue
        control_abunds = [c.get_count(kmer) for c in controls]
        ctl_thresh_pass = [a > args.control_threshold for a in control_abunds]
        if False in ctl_thresh_pass:
            continue

        ctl_abund_str = '\t'.join([str(a) for a in control_abunds])
        print(kmer, case_abund, ctl_abund_str)
