import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed
from screed.fastq import fastq_iter

import khmer

hasher = khmer.HashtableIntersect(17, *khmer.PRIMES_1b)

filename = '/scratch/gpgc/iowa/850.2.fq'

for n, record in enumerate(fastq_iter(open(filename))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n
    if n > 10**6:
        break
        
    sequence = record['sequence']
    if 'N' in sequence:
        continue

    hasher.consume(sequence)

for n, record in enumerate(fastq_iter(open(filename))):
    if n % 10000 == 0:
        print>>sys.stderr, '2...', n
    if n > 10**6:
        break
        
    sequence = record['sequence']
    if 'N' in sequence:
        continue

    if hasher.get_min_count(sequence) > 5:
        print record['name']
