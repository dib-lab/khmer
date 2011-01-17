#! /usr/bin/env python
import sys
from screed.fasta import fasta_iter

n = 0
for filename in sys.argv[1:]:
    sys.stderr.write('... %s %d\n' % (filename, n))
    idx = filename.find('group')
    assert idx != -1, filename

    group_num = int(filename[idx + 5:].split('.')[0])

    for record in fasta_iter(open(filename + '/contigs.fa')):
        print n, group_num, len(record['sequence'])
        n += 1
