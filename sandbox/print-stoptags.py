#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import os

K = 32


def main():
    ht = khmer.new_hashbits(32, 1, 1)
    ht.load_stop_tags(sys.argv[1])
    ht.print_stop_tags(os.path.basename(sys.argv[1]) + '.txt')


if __name__ == '__main__':
    main()
