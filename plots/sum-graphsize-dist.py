#! /usr/bin/env python
import sys

total = 0
totaltotal = int(sys.argv[2])
for line in open(sys.argv[1]):
    size, num, sumnum = line.split()
    size = int(size)
    num = int(num)
    total += num*size

    print size, num, sumnum, total, float(total) / float(totaltotal) * 100.
