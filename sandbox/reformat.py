import sys
import screed

for record in screed.open(sys.argv[1], parse_description=False):
    name, descr = record.name.split(" ", 1)
    name = name.lstrip('>')
    name += "\t".join(descr.split())

    print '>%s\n%s' % (name, record.sequence,)
