#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import screed
import sys

CUTOFF = 200


def main():
    n = 0
    prefix = sys.argv[1]
    for filename in sys.argv[2:]:
        for record in screed.open(filename):
            if len(record.sequence) >= CUTOFF:
                n += 1
                print('>%s.%s %s\n%s' % (prefix, n, record.name, record.sequence))


if __name__ == '__main__':
    main()
