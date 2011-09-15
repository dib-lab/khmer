#! /usr/bin/env python
import screed
import sys

infile = sys.argv[1]
outfile = infile
if len(sys.argv) > 2:
    outfile = sys.argv[2]

single_fp = open(outfile + '.se', 'w')
paired_fp = open(outfile + '.pe', 'w')

last_record = None
last_name = None
for record in screed.open(sys.argv[1]):
    name = record['name'].split()[0]
    sequence = record['sequence']

    if last_record:
        if last_name.endswith('/1') and name.endswith('/2') and name[:-1] == last_name[:-1]:
           fp = paired_fp
           print >>paired_fp, '>%s\n%s' % (last_name, last_record['sequence'])
           print >>paired_fp, '>%s\n%s' % (name, record['sequence'])
           name, record = None, None
        else:
           print >>single_fp, '>%s\n%s' % (last_name, last_record['sequence'])

    last_name = name
    last_record = record

if last_record:
   if last_name.endswith('/1') and name.endswith('/2') and name[:-1] == last_name[:-1]:
      fp = paired_fp
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
