#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import sys
import screed


def main():
    filename = sys.argv[1]
    prefix = sys.argv[2]
    size = int(float(sys.argv[3]))          # e.g. 1e9

    division = -1
    for n, record in enumerate(screed.open(filename)):
        if n % 100000 == 0:
            print('...', n)

        if n % size == 0:
            division += 1
            new_name = '%s.%04d.fa' % (prefix, division)
            print('opening', new_name)
            fp = open(new_name, 'w')

        fp.write('>%s\n%s\n' % (record['name'], record['sequence']))


if __name__ == '__main__':
    main()
