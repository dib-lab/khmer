#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys

import khmer
import screed

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--controls', metavar='FILE', nargs='+',
                        help='list of countgraph files corresponding to '
                        'control sample(s)')
    parser.add_argument('--case', metavar='FILE', help='countgraph of '
                        'sample from case/condition of interest')
    parser.add_argument('-x', '--control_threshold', metavar='X', type=int,
                        default=1, help='k-mers with abund > X in any control '
                        'sample are uninteresting; default=1')
    parser.add_argument('-y', '--case_lower_threshold', metavar='Y', type=int,
                        default=5, help='ignore k-mers from case with abund '
                        '< Y; default=5')
    parser.add_argument('--dry-run', action='store_true',
                        help='do not execute, just print arguments')
    parser.add_argument('case_fastq')
    return parser

def main(args):
    if args.dry_run:
        print(*sys.argv)
        exit(0)

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

    print('Iterating over case reads', args.case_fastq, '...', file=sys.stderr)
    intreadcount = 0
    intkmercount = 0
    lpaths = dict()
    for n, record in enumerate(screed.open(args.case_fastq)):
        if n > 0 and n % 1e6 == 0:
            print('    processed', n, 'reads...', file=sys.stderr)
        novel_kmers = dict()
        for i, kmer in enumerate(case.get_kmers(record.sequence)):
            case_abund = case.get(kmer)
            if case_abund < args.case_lower_threshold:
                continue
            control_abunds = [c.get(kmer) for c in controls]
            ctl_thresh_pass = [a > args.control_threshold for a in control_abunds]
            if True in ctl_thresh_pass:
                continue

            novel_kmers[i] = [kmer, case_abund] + control_abunds
            linear_path = case.assemble_linear_path(kmer)
            if linear_path not in lpaths:
                lpaths[linear_path] = []
            lpaths[linear_path].append(record.name)
            intkmercount += 1
        if len(novel_kmers) > 0:
            intreadcount += 1
            print('@', record.name, '\n', record.sequence, '\n+\n', record.quality, sep='')
            for i in sorted(novel_kmers):
                kmer, case_abund, ctl1, ctl2 = novel_kmers[i]
                abundstr = ' '.join([str(abund) for abund in [case_abund, ctl1, ctl2]])
                print(' ' * i, kmer, ' ' * 10, abundstr, '#', sep='')
                sys.stdout.flush()

    message = 'Found {} novel kmers in {} reads'.format(intkmercount, intreadcount)
    message += ', {} linear paths'.format(len(lpaths))
    print(message, file=sys.stderr)
    for i, linear_path in enumerate(lpaths):
        readnames = lpaths[linear_path]
        lpathname = '##>lpath{} {} reads {}'.format(
                        i+1, len(readnames),
                        ' '.join(readnames)
                    )
        print('##>', lpathname, '\n', '##', linear_path, sep='')

if __name__ == '__main__':
    main(get_parser().parse_args())
