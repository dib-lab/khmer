#! /usr/bin/env python
import sys
import screed
import os.path

for filename in sys.argv[1:]:
    assert 'R1' in filename
    filename2 = filename.replace('R1', 'R2')

    r1 = iter(screed.open(filename)).next()
    r2 = iter(screed.open(filename2)).next()

    assert r1.name == r2.name, (r1.name, r2.name)

    final = filename.replace('R1', '')
    print 'python /root/khmer/sandbox/interleave.py %s %s | gzip -9c > %s' % (filename, filename2, final)
