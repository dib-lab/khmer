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


def main():
    for record in screed.open(sys.argv[1]):
        name = record['name']
        sequence = record['sequence']

        name = name.split()[0]

        print('>%s\n%s' % (name, sequence,))


if __name__ == '__main__':
    main()
