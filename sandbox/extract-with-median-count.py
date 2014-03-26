#! /usr/bin/env python
import screed
import sys
import argparse


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('min_count', type=int)
    parser.add_argument('max_count', type=int)
    parser.add_argument('filenames', nargs='+')

    args = parser.parse_args()
    min_count = args.min_count
    max_count = args.max_count

    for filename in args.filenames:
        for n, record in enumerate(screed.open(filename, parse_description=False)):
            if n % 10000 == 0:
                print >>sys.stderr, '...', filename, n
            kmed = record.name.split()[-1]
            assert kmed.startswith('kmed'), record.name
            kmed = kmed.split('=')[1]
            kmed = int(kmed)

            if kmed >= min_count and kmed < max_count:
                print '>%s\n%s' % (record.name, record.sequence)

if __name__ == '__main__':
    main()
