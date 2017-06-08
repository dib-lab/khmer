#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys

from Bio import SeqIO
import pysam

def match(record, seq):
    refrseq = str(seq[record.pos:record.pos+record.rlen].seq)
    return refrseq == records.seq

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seqid', metavar='SEQ',
                        help='dump reads not mapped to SEQ')
    parser.add_argument('--out', type=argparse.FileType('w'),
                        help='output file; default is terminal (stdout)')
    parser.add_argument('--dry-run', action='store_true',
                        help='do not execute, just print arguments')
    parser.add_argument('fasta')
    parser.add_argument('bam')
    return parser

def main(args):
    if args.dry_run:
        print(*sys.argv)
        exit(0)

    with open(args.fasta, 'r') as seqfile:
        seqs = {record.id: record for record in SeqIO.parse(seqfile, 'fasta')}

    bam = pysam.AlignmentFile(args.bam, 'rb')
    for i, record in enumerate(bam):
        if (i+1) % 10000 == 0:
            print('...processed', i, 'records', file=sys.stderr)

        if args.seqid:
            if record.reference_id < 0 or bam.get_reference_name(record.reference_id) != args.seqid:
                continue

        if record.is_secondary or record.is_supplementary:
            continue

        matchcigar = '{:d}M'.format(record.rlen)
        if record.cigarstring != matchcigar:
            continue

        seq = seqs[bam.get_reference_name(record.tid)]
        refrseq = str(seq[record.pos:record.pos+record.rlen].seq)
        if refrseq.upper() != record.seq.upper():
            print('@', record.qname, '\n', record.seq, '\n+\n', record.qual, sep='', file=args.out)

if __name__ == '__main__':
    main(get_parser().parse_args())
