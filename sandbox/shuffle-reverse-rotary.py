from __future__ import print_function
#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
import sys
import screed
import os.path
import argparse

ROTARY_SIZE = 100


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(filenames, nargs='+')
    args = parser.parse_args()

    prefix = os.path.basename(args.filenames[0])

    fp_d = {}
    for n in range(0, ROTARY_SIZE):
        num = ROTARY_SIZE - n
        fp_d[n] = open(prefix + '.%03d' % num, 'w')

    total = 0
    for filename in args.filenames:
        for record in screed.open(filename):
            total += 1
            if total % 10000 == 0:
                print('...', total)
            loc = total % ROTARY_SIZE
            fp_d[loc].write('>%s\n%s\n' % (record.name, record.sequence))

    print('reverse-rotary shuffled %d sequences into %d files (%s.NNN)' % \
        (total, ROTARY_SIZE, prefix))


if __name__ == '__main__':
    main()
