import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed
from screed.fastq import fastq_iter
from screed.fasta import fasta_iter

K = 15

import khmer

kh = khmer.new_hashtable(K, 1999999973)

filename = '/scratch/gpgc/iowa/850.2.fq'
#filename = 'foo'

for n, record in enumerate(fastq_iter(open(filename))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n
    if n > 10**6:
        break

    sequence = record['sequence']
    if 'N' in sequence:
        continue

    kh.consume(sequence)

bins = [0] * 256
for n, record in enumerate(fastq_iter(open(filename))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n
    if n > 10**5:
        break

    sequence = record['sequence']
    if 'N' in sequence:
        continue

    for i in range(0, len(sequence) - K + 1):
        subseq = sequence[i:i+K]

        abundance = kh.get(subseq)
        bins[abundance] += 1

for n, count in enumerate(bins):
    print n, count
