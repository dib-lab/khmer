#! /usr/bin/env python
import screed, sys, itertools

s1_file = sys.argv[1]
s2_file = sys.argv[2]

for r1, r2 in itertools.izip(screed.open(s1_file), screed.open(s2_file)):
    print '@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s' % (r1.name,
                                              r1.sequence,
                                              r1.accuracy,
                                              r2.name,
                                              r2.sequence,
                                              r2.accuracy)
