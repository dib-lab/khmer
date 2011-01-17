#! /usr/bin/env python
import screed
import sys

single_fp = open(sys.argv[2] + '.se', 'w')
paired_fp = open(sys.argv[2] + '.pe', 'w')

last_record = None
last_name = None
for record in screed.fasta.fasta_iter(open(sys.argv[1]), parse_description=False):
    name = record['name']
    sequence = record['sequence']

    name = name.split()[0]

    fp = single_fp
    if name.endswith('2') and last_record and name[:-1] == last_name[:-1]:
        fp = paired_fp
        print >>paired_fp, '>%s\n%s' % (last_name, last_record['sequence'])
    
    print >>fp, '>%s\n%s' % (name, sequence,)
    last_name = name
    last_record = record
