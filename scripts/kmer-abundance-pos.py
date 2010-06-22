import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed
from screed.fastq import fastq_iter
import glob

import khmer

K = 17

hasher = khmer.HashtableIntersect(K, *khmer.PRIMES_1b)

files = glob.glob(sys.argv[1])

for filename in files:
    for n, record in enumerate(fastq_iter(open(filename))):
        if n % 10000 == 0:
            print>>sys.stderr, '...', n
        #if n > 10**5:
        #    break
        
        sequence = record['sequence']
        if 'N' in sequence:
            continue

        hasher.consume(sequence)

readAbund = [0]*98

for filename in files:
    for n, record in enumerate(fastq_iter(open(filename))):
        if n % 10000 == 0:
            print>>sys.stderr, '2...', n
        #if n > 10**5:
        #    break
        
        sequence = record['sequence']
        if 'N' in sequence:
            continue

        for i in range(len(sequence) - K + 1):
            kmer = sequence[i:(i+K)]
            readAbund[i] += hasher.get_min_count(kmer)

print readAbund
