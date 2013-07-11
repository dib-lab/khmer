#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print >>sys.stderr, "USAGE: plot_two_column_bar.py <input_file> [title]"
    sys.exit(1)

x = []
y = []

xlab = None
ylab = None

for line in open(sys.argv[1]):
    if line[0] == "#":
        lexemes = line[1:].strip().split()
        xlab = lexemes[0]
        ylab = lexemes[1]
    else:
        lexemes = line.strip().split()
        if len(lexemes) > 1:
            x.append(int(lexemes[0]))
            y.append(float(lexemes[1]))

print "Bins: %d" % max(y)
plt.hist(y, max(y))
if len(sys.argv) > 3:
    plt.title(sys.argv[2])

if xlab:
    plt.xlabel(ylab)
    plt.ylabel("count")

plt.savefig("%s_histo.png" % sys.argv[1])
