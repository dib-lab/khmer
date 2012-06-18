#! /usr/bin/env python
import screed, sys, itertools

s1_file = sys.argv[1]
s2_file = sys.argv[2]

for r1, r2 in itertools.izip(screed.open(s1_file), screed.open(s2_file)):
    name1 = r1.name
    if not name1.endswith('/1'):
        name1 += '/1'
    name2 = r2.name
    if not name2.endswith('/2'):
        name2 += '/2'

    print '@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s' % (name1,
                                              r1.sequence,
                                              r1.accuracy,
                                              name2,
                                              r2.sequence,
                                              r2.accuracy)
