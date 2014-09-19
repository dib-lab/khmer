#! /usr/bin/env python
import argparse
import screed
import sys
import khmer


def output_single(read):
    if hasattr(read, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.accuracy)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min-coverage', type=int, default=None)
    parser.add_argument('-M', '--max-coverage', type=int, default=None)
    parser.add_argument('input_counting_table')
    parser.add_argument('input_readfile')
    parser.add_argument('output_readfile')
    args = parser.parse_args()

    print >>sys.stderr, 'min_coverage: %s' % args.min_coverage
    print >>sys.stderr, 'max_coverage: %s' % args.max_coverage

    if not (args.min_coverage or args.max_coverage):
        print >>sys.stderr, "neither min nor max coverage specified!? exiting!"
        sys.exit(1)

    if args.min_coverage and args.max_coverage and \
       args.max_coverage < args.min_coverage:
        print >>sys.stderr, "min_coverage > max_coverage!? exiting!"
        sys.exit(1)

    htable = khmer.load_counting_hash(args.input_counting_table)
    output_file = args.output_readfile
    output_fp = open(output_file, 'w')

    n_kept = 0
    n = 0
    for n, record in enumerate(screed.open(args.input_readfile)):
        if n % 100000 == 0:
            print >>sys.stderr, '...', n, n_kept

        seq = record.sequence.upper()
        if 'N' in seq:
            seq = seq.replace('N', 'G')

        try:
            med, _, _ = htable.get_median_count(seq)
        except ValueError:
            continue

        keep = True
        if args.min_coverage and med < args.min_coverage:
            keep = False

        if args.max_coverage and med > args.max_coverage:
            keep = False

        if keep:
            n_kept += 1

            output_fp.write(output_single(record))

    print >>sys.stderr, 'consumed %d reads; kept %d' % (n, n_kept)

if __name__ == '__main__':
    main()
