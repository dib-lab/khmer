#! /usr/bin/env python
import khmer, sys, os

K = 32

ht = khmer.new_hashbits(32, 1, 1)
ht.load_stop_tags(sys.argv[1])
ht.print_stop_tags(os.path.basename(sys.argv[1]) + '.txt')
