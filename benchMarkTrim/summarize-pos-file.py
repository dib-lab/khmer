#! /usr/bin/env python2
import sys
import argparse
import screed

def read_pos_file(filename):
    for line in open(filename):
        line = line.strip()
        try:
            read, posns = line.split(' ', 2)
            posns = map(int, posns.split(','))
        except ValueError:
            read = line
            posns = []
            continue
            
        yield read, posns

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('posfile')
    parser.add_argument('reads')
    parser.add_argument('--limitreads', default=None)
    args = parser.parse_args()

    print 'reading files...', args.posfile, args.reads
    posdict = dict(read_pos_file(args.posfile))

    limitnames = None
    if args.limitreads:
        limitnames = set([ readname for readname, _ in \
                           read_pos_file(args.limitreads) ])
    
    all_reads = 0
    sum_bp = 0

    print 'reading sequences...'
    for n, record in enumerate(screed.open(args.reads)):
        if n % 100000 == 0:
            print >>sys.stderr, '...', n

        if args.limitreads and record.name not in limitnames:
            continue

        all_reads += 1
        sum_bp += len(record.sequence)

    print 'done!'

    n_reads = 0
    n = 0
    m = 0
    skipped = 0
    for k, v in posdict.iteritems():
        if args.limitreads and k not in limitnames:
            skipped += 1
            continue

        n_reads += 1

        if not v:
            continue

        n += 1
        m += len(v)

    print 'XXX', all_reads, n_reads

    print 'posfile %s: %d mutated reads of %d; %d mutations total' % \
          (args.posfile, n, n_reads, m)
    print 'skipped:', skipped
    print '%d bp total' % (sum_bp,)
    print 'overall error rate: %f%%' % (100. * m / float(sum_bp))


if __name__ == '__main__':
    main()


