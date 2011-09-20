#! /usr/bin/env python
import screed
import sys

def is_pair(name1, name2):
    if name1.endswith('/1') and name2.endswith('/2'):
        s1 = name1.split('/')[0]
        s2 = name2.split('/')[0]
        if s1 == s2:
            assert(s1)
            return True
        
    return False

infile = sys.argv[1]
outfile = infile
if len(sys.argv) > 2:
    outfile = sys.argv[2]

single_fp = open(outfile + '.se', 'w')
paired_fp = open(outfile + '.pe', 'w')

last_record = None
last_name = None

print 'splitting pe/se sequences from %s to %s.{pe,se}' % (infile, outfile)
for n, record in enumerate(screed.open(sys.argv[1])):
    if n % 100000 == 0 and n > 0:
       print '...', n
    name = record['name'].split()[0]
    sequence = record['sequence']

    if last_record:
        if is_pair(last_name, name):
           print >>paired_fp, '>%s\n%s' % (last_name, last_record['sequence'])
           print >>paired_fp, '>%s\n%s' % (name, record['sequence'])
           name, record = None, None
        else:
           print >>single_fp, '>%s\n%s' % (last_name, last_record['sequence'])

    last_name = name
    last_record = record

if last_record:
    if is_pair(last_name, name):
        print >>paired_fp, '>%s\n%s' % (last_name, last_record['sequence'])
        print >>paired_fp, '>%s\n%s' % (name, record['sequence'])
        name, record = None, None
    else:
        print >>single_fp, '>%s\n%s' % (last_name, last_record['sequence'])
        name, record = None, None

if record:
   print >>single_fp, '>%s\n%s' % (name, record['sequence'])

single_fp.close()
paired_fp.close()

### check, at the end, to see if it worked!
paired_fp = open(outfile + '.pe')
if not paired_fp.read(1):
    raise Exception("no paired reads!? check file formats...")
