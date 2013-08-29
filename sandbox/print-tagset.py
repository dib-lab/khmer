#! /usr/bin/env python
import khmer
import sys
import os

K = 32

ht = khmer.new_hashbits(32, 1, 1)
ht.load_tagset(sys.argv[1])
print 'loaded!'
ht.print_tagset(os.path.basename(sys.argv[1]) + '.txt')
