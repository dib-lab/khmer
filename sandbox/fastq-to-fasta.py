import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed
from screed.fastq import fastq_iter

for n, record in enumerate(fastq_iter(open(sys.argv[1]))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n

    sequence = record['sequence']
    name = record['name']

    if 'N' in sequence:
        continue

    print ">" + name
    print sequence
